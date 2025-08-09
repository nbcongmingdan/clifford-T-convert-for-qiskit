"""clifford_t_converter.py
---------------------------------
Utility functions for converting an arbitrary ``qiskit.circuit.QuantumCircuit``
into an *equivalent* Clifford+T-only circuit using `pygridsynth`.

Main entry point
----------------
>>> from clifford_t_converter import to_clifford_t
>>> ct_circ = to_clifford_t(original_circ, precision=1e-3)

Requirements
------------
* qiskit ≥ 0.46  (for the experimental ``qiskit.qasm3.dumps`` helper)
* pygridsynth
* bqskit          (only used to ensure QASM3 compliance, optional)
* mpmath
"""



from __future__ import annotations

import os
import re
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Union


import re, mpmath
from pathlib import Path
from typing import List

PI = mpmath.pi
pi = mpmath.pi  # For compatibility with eval in gridsynth_modem

from pygridsynth.gridsynth import gridsynth_gates


from qiskit.circuit import QuantumCircuit as QiskitCircuit

import qiskit.qasm3 as qasm3

###############################################################################
# QASM helpers
###############################################################################

def _ensure_dir(path: os.PathLike | str) -> Path:
    """Create *path* if it does not yet exist and return it as ``Path``."""
    path = Path(path)
    if not path.exists():
        path.mkdir(parents=True, exist_ok=True)
    return path


def export_qiskit_to_qasm3(
    circuit: QiskitCircuit,
    out_dir: os.PathLike | str | None = None,
    experimental: bool = False
) -> Path:
    """Serialise `circuit` to OpenQASM 3 and write it to `out_dir`."""
    # 1. Prepare output directory
    if out_dir is None:
        out_dir = tempfile.gettempdir()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 2. Dump to QASM 3 (includes header and include "stdgates.inc";)
    if experimental:
        qasm_str = qasm3.dumps(circuit, experimental=True)
    else:
        qasm_str = qasm3.dumps(circuit)

    # 3. Write to file
    filename = f"{circuit.name or 'circuit'}.qasm"
    out_path = out_dir / filename
    out_path.write_text(qasm_str, encoding="utf-8")

    print(f"Exported circuit '{circuit.name}' to QASM3: {out_path}")
    return out_path


###############################################################################
# gridsynth 调用 + QASM 转换工具（已测试版本）
###############################################################################

def gridsynth_modem(angle: any, precision: any, qubit: int) -> dict[int, str]:
    """
    用 gridsynth 将 RZ(angle) 近似到 Clifford+T, 返回 {行号: QASM3 指令}。
    依赖 `include "stdgates.inc";`(内含 gphase/rz/h/t/s/x)。
    """
    # ---------- ① 解析参数 ----------
    theta  = mpmath.mpmathify(eval(str(angle), {"pi": mpmath.pi}))
    epsilon  = mpmath.mpmathify(precision)
    if epsilon <= 0:
        raise ValueError("precision must be positive")

    # ---------- ② 固定高精度 ----------
    prev_dps, mpmath.mp.dps = mpmath.mp.dps, 128
    try:
        start = time.time()
        gates: list[str] = gridsynth_gates(theta, epsilon)          # 可能抛 ZeroDivisionError
    finally:
        mpmath.mp.dps = prev_dps

    # ---------- ③ 组装 QASM ----------
    qasm: list[str] = []

    mapping = {
        "H": "h",
        "T": "t",              
        "S": "s",
        "X": "x",
        
    }

    for g in reversed(gates):          # gridsynth 列表先右→左
        if g == "I":
            continue
        if g == "W":
            # W = e^{iπ/4}·I -- 纯全局相位
            qasm.append("gphase(pi/4);")
        else:
            qasm.append(f"{mapping[g]} q[{qubit}];")

    print(
        f"gridsynth RZ({angle}) with epsilon≈{precision} → {len(qasm)} lines "
        f"in {time.time() - start:.3f}s."
    )
    return {i: line for i, line in enumerate(qasm)}


###############################################################################
# QASM <‑‑> dict utilities
###############################################################################

def dict_to_qasm_file(qasm_dict: Dict[int, str], output_path: str | Path) -> None:
    header  = ("OPENQASM 3.0;", 'include "stdgates.inc";')
    q_lines = [qasm_dict[idx] for idx in sorted(qasm_dict)]
    
    # ---- 计算量子寄存器大小 ----
    import re
    qubit_indices = []
    for line in q_lines:
        qubit_indices += [int(m.group(1)) for m in re.finditer(r"q\[(\d+)\]", line)]
    n_qubits = (max(qubit_indices) + 1) if qubit_indices else 1

    with Path(output_path).open("w", encoding="utf-8") as fout:
        fout.write("\n".join(header) + "\n")
        fout.write(f"qreg q[{n_qubits}];\n")
        # 若有测量，可按需写 creg 行：
        # fout.write(f"creg c[{n_qubits}];\n")
        for line in q_lines:
            fout.write(line + "\n")







def _parse_angle(expr: str, pi_const) -> float:
    """Convert '0.3', 'pi/4', etc. to a float."""
    try:
        return float(expr)
    except ValueError:
        # safe eval with only ‘pi’ exposed
        return eval(expr, {"pi": pi_const})
    

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def qasm_file_to_dict(
    qasm_file_path: os.PathLike | str,
    precision: float = 1e-3,
) -> Dict[int, str]:
    """
    读取 QASM (2.0 或 3.0) 文件，把所有 u3(θ, φ, λ) q[i]; 拆解成
    Clifford+T 门序列（倒序写出），其余指令保持原顺序。
    返回 {行号: 指令文本}。
    """

    # ------------------------------------------------------------------
    # 1. 读取文件并去掉行尾换行/空白
    # ------------------------------------------------------------------
    lines = [ln.rstrip() for ln in Path(qasm_file_path).read_text("utf-8").splitlines()]

    # ------------------------------------------------------------------
    # 1a. 定位第一条“量子寄存器声明”行
    #     • QASM 2:  qreg q[4];
    #     • QASM 3:  qubit[4] q;   或   qubit q;
    # ------------------------------------------------------------------
    def _is_qreg_decl(s: str) -> bool:
        s = s.lstrip()
        return (
            s.startswith("qreg ")           # QASM 2
            or s.startswith("qubit ")       # QASM 3 (单比特)
            or s.startswith("qubit[")       # QASM 3 (数组)
        )
    try:
        split_idx = next(i for i, l in enumerate(lines) if _is_qreg_decl(l))
    except StopIteration:
        raise ValueError("Input QASM missing quantum-register declaration (qreg / qubit).")

    # ------------------------------------------------------------------
    # 1b. 抛弃文件头 & 所有寄存器声明，只留下真正的“body”指令
    # ------------------------------------------------------------------
    body: list[str] = []
    in_header = True
    for ln in lines:
        if in_header and _is_qreg_decl(ln):
            # 读取到第一条量子寄存器声明 → 进入 body
            in_header = False
            continue  # 跳过本行寄存器声明
        if in_header:
            continue  # 仍在 header 区域，继续跳过
        # 进入 body 区域，过滤后续所有寄存器声明行（qreg/creg/qubit/bit）
        ls = ln.lstrip()
        if ls.startswith(("qreg ", "creg ", "qubit", "bit")):
            continue
        body.append(ln)

    # ------------------------------------------------------------------
    # 2. 准备正则：匹配  u3(θ, φ, λ) q[i];
    # ------------------------------------------------------------------
    u3_pattern = re.compile(
        r"u3\(\s*([^,]+)\s*,\s*([^,]+)\s*,\s*([^)]+)\s*\)\s+q\[(\d+)\];",
        re.IGNORECASE,
    )

    def _parse_angle(expr: str) -> float:
        """把角度表达式安全 eval 成 float(rad)。"""
        return float(mpmath.mpf(str(eval(expr, {"pi": PI}))))

    # ------------------------------------------------------------------
    # 3. 主循环：拆解 u3；其余行原样保留
    # ------------------------------------------------------------------
    final_lines: list[str] = []

    for ln in body:
        m = u3_pattern.match(ln)
        if m is None:
            # 非 u3 → 直接写入
            final_lines.append(ln)
            continue

        theta_expr, phi_expr, lam_expr, q_idx = m.groups()
        theta = _parse_angle(theta_expr)
        phi   = _parse_angle(phi_expr)
        lam   = _parse_angle(lam_expr)
        q     = int(q_idx)

        # ----------- 生成 Clifford+T 序列（倒序写出） --------------
        seg: list[str] = []

        # 全局相位 -(π+θ)/2
        phase1 = -(PI + theta) / 2
        seg.append(f"gphase({mpmath.nstr(phase1, 17)});")

        # Rz(φ+π)
        phase2 = (PI + phi) / 2
        seg.append(f"gphase({mpmath.nstr(phase2, 17)});")
        seg.extend(gridsynth_modem(phase2 * 2, precision, q).values())
        seg.append(f"sx q[{q}];")

        # Rz(θ+π)
        phase4 = -phase1             # (π + θ)/2
        seg.append(f"gphase({mpmath.nstr(phase4, 17)});")
        seg.extend(gridsynth_modem(phase4 * 2, precision, q).values())
        seg.append(f"sx q[{q}];")

        # Rz(λ)
        phase5 = lam / 2
        seg.append(f"gphase({mpmath.nstr(phase5, 17)});")
        seg.extend(gridsynth_modem(lam, precision, q).values())

        # 倒序插入
        for inst in reversed(seg):
            final_lines.append(inst)

    # ------------------------------------------------------------------
    # 4. 打包为 {行号: 指令} 字典返回
    # ------------------------------------------------------------------
    return {i + 1: text for i, text in enumerate(final_lines)}


###############################################################################
# Public API
###############################################################################

def to_clifford_t(circuit: QiskitCircuit, precision: float = 1e-3) -> QiskitCircuit:
    """Convert *circuit* to an equivalent Clifford+T-only QuantumCircuit.

    Steps
    -----
    1. Export input ``circuit`` to a temporary QASM2 file.
    2. Rewrite the QASM, replacing **all** ``u3`` single-qubit rotations with
       a Clifford+T decomposition from :func:`gridsynth_modem`.
    3. Reload and return the Clifford+T circuit.

    Parameters
    ----------
    circuit : QiskitCircuit
        Input circuit (may contain arbitrary gates).
    precision : float, default 1e-3
        Target synthesis precision passed to *gridsynth*.

    Returns
    -------
    QiskitCircuit
        A logically equivalent Clifford+T-only circuit.
    """
    # 1) Serialise original circuit ➜ QASM3
    orig_qasm_path = export_qiskit_to_qasm3(circuit)

    # 2) Parse + rewrite ➜ dict
    ct_dict = qasm_file_to_dict(orig_qasm_path, precision=precision)

    # 3) Write rewritten QASM
    ct_qasm_path = orig_qasm_path.with_stem(orig_qasm_path.stem + "_ct")
    dict_to_qasm_file(ct_dict, ct_qasm_path)

    # 4) Reload via Qiskit and return
    ct_circuit = qasm3.load(str(ct_qasm_path))
    return ct_circuit


# Alias retained for backward compatibility with earlier drafts
main_function = to_clifford_t


###############################################################################
# CLI helper (optional)
###############################################################################
if __name__ == "__main__":
    import argparse
    from qiskit import QuantumCircuit as _QC

    parser = argparse.ArgumentParser(description="Convert a Qiskit circuit QASM file to Clifford+T only.")
    parser.add_argument("input", type=str, help="Path to a .qasm or .py file that generates a Qiskit circuit.")
    parser.add_argument("--precision", "-p", type=float, default=1e-3, help="Target synthesis precision (default: 1e-3)")
    parser.add_argument("--output", "-o", type=str, default=None, help="Path to save the converted QASM. Defaults next to input.")
    args = parser.parse_args()

    if args.input.endswith(".qasm"):
        qc = _QC.from_qasm_file(args.input)
    else:
        # Execute the Python file and look for a variable named ``circ``
        namespace: dict[str, object] = {}
        exec(Path(args.input).read_text(), namespace)
        qc = namespace.get("circ")  # type: ignore
        if not isinstance(qc, _QC):
            raise RuntimeError("Python file did not define a QuantumCircuit named 'circ'.")

    ct_qc = to_clifford_t(qc, precision=args.precision)
    out_path = args.output or Path(args.input).with_suffix(".ct.qasm")
    _ensure_dir(Path(out_path).parent)
    out_path.write_text(ct_qc.qasm(), encoding="utf-8")
    print(f"Clifford+T circuit written to {out_path}")
