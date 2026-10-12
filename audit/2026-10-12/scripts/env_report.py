"""Record the exact computational environment (run inside the WSL2 project folder).

    python env_report.py > env_report.json

Read-only: installs nothing, changes nothing. Captures interpreter, OS/WSL details, CPU/RAM,
thread settings, scientific and quantum package versions, BLAS/LAPACK build info, PySCF
runtime settings and `pip freeze`, so results can be tied to the environment that produced them.
"""
import importlib
import importlib.metadata as md
import io
import json
import os
import platform
import subprocess
import sys
from contextlib import redirect_stdout

DISTS = ["numpy", "scipy", "pyscf", "h5py", "matplotlib", "pytest", "geometric", "pyberny",
         "qiskit", "qiskit-aer", "qiskit-nature", "qiskit-algorithms", "openfermion",
         "openfermionpyscf", "pennylane", "cirq-core", "ffsim", "rdkit"]


def _read(path):
    try:
        with open(path) as f:
            return f.read()
    except OSError:
        return ""


def main():
    out = {"python": sys.version, "executable": sys.executable,
           "platform": platform.platform(), "machine": platform.machine()}
    proc = _read("/proc/version").strip()
    out["kernel"] = proc
    out["is_wsl"] = "microsoft" in proc.lower()
    out["cwd"] = os.getcwd()
    out["cwd_on_windows_drive(/mnt/*: slow I/O under WSL2)"] = os.getcwd().startswith("/mnt/")
    out["cpu_count"] = os.cpu_count()
    mem = {l.split(":")[0]: l.split(":")[1].strip() for l in _read("/proc/meminfo").splitlines() if ":" in l}
    out["mem_total"] = mem.get("MemTotal")
    out["swap_total"] = mem.get("SwapTotal")
    out["env"] = {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                                                  "PYSCF_TMPDIR", "PYSCF_MAX_MEMORY", "TMPDIR", "VIRTUAL_ENV", "CONDA_PREFIX")}
    versions = {}
    for d in DISTS:
        try:
            versions[d] = md.version(d)
        except md.PackageNotFoundError:
            versions[d] = None
    out["package_versions"] = versions
    try:
        import numpy
        buf = io.StringIO()
        with redirect_stdout(buf):
            numpy.show_config()
        out["numpy_show_config"] = buf.getvalue()
    except Exception as exc:
        out["numpy_show_config"] = repr(exc)
    try:
        pyscf = importlib.import_module("pyscf")
        from pyscf import lib
        out["pyscf"] = {"version": pyscf.__version__, "path": pyscf.__file__,
                        "TMPDIR": lib.param.TMPDIR, "MAX_MEMORY_MB": lib.param.MAX_MEMORY,
                        "num_threads": lib.num_threads()}
        try:
            importlib.import_module("pyscf.dft.libxc")
            out["pyscf"]["xc_library"] = "libxc"
        except Exception:
            out["pyscf"]["xc_library"] = "libxc unavailable"
    except Exception as exc:
        out["pyscf"] = repr(exc)
    try:
        out["pip_freeze"] = subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True,
                                           text=True, timeout=120).stdout.splitlines()
    except Exception as exc:
        out["pip_freeze"] = repr(exc)
    try:
        out["git_head"] = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        out["git_dirty_files"] = len(subprocess.run(["git", "status", "--porcelain"], capture_output=True,
                                                    text=True).stdout.splitlines())
    except Exception:
        out["git_head"] = None
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
