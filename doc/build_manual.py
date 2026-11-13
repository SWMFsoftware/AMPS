#!/usr/bin/env python3
"""Build the local manual without importing AMPS/SWMF make configuration.

The standard-library-only driver discovers compiler names on PATH and accepts
an optional executable path, never shell text. Compiler locations are discovered
on each machine, not read back from an earlier build receipt. It keeps products
below doc/, preserves failed logs, and
publishes AMPS.pdf only after successful compilation and reference checks.
It does not execute the scientific commands printed in the manual.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
BUILD = ROOT / "build-manual"


def run(argv, env, log_path, timeout=600):
    """Capture the complete tool transcript; a failed pass never replaces PDF."""
    print("Running:", " ".join(str(x) for x in argv), flush=True)
    with log_path.open("w", encoding="utf-8") as log:
        try:
            p = subprocess.run(
                [str(x) for x in argv], cwd=str(ROOT), env=env,
                stdout=log, stderr=subprocess.STDOUT, timeout=timeout)
        except subprocess.TimeoutExpired:
            raise RuntimeError("manual compiler timed out; log: " + str(log_path))
    if p.returncode:
        print(log_path.read_text(errors="replace")[-12000:])
        raise RuntimeError("compiler exit {}; log: {}".format(p.returncode, log_path))


def choose_engine(request):
    """Return (invocation, identity path), resolving tools on this machine.

    A normal PATH installation is invoked as pdflatex/tectonic, so the build
    command needs no machine-specific compiler prefix. The separately returned
    absolute path is only for provenance/hashing. Preserve a symlink's public
    name: TeX installations commonly use the invoked name to choose a format.
    Explicit paths and relative PATH entries need an absolute invocation,
    because compilation itself runs from the manual directory, not the caller's
    original working directory. These paths are computed, never hard-coded.
    """
    def installed(name):
        candidate = shutil.which(name)
        if not candidate:
            return None
        identity = Path(candidate).absolute()
        invocation = (name if Path(name).name == name and
                      Path(candidate).is_absolute() else str(identity))
        return invocation, identity

    if request != "auto":
        # Expand ~ only for an explicitly supplied path. The default search
        # never assumes a home directory, TeX Live version or installation root.
        request = os.path.expanduser(request)
        candidate = installed(request)
        if candidate:
            return candidate
        if Path(request).is_file():
            identity = Path(request).absolute()
            return str(identity), identity
        # The documented explicit name must find the same project-local
        # Tectonic as auto mode. A local installation is intentionally not
        # added to the user's PATH; this also makes --only-cached usable
        # without changing the shell environment. Explicit paths and requests
        # for other engines never silently select a different compiler.
        local = ROOT / ".build-tools" / "tectonic"
        if request == "tectonic" and local.is_file():
            return str(local), local
        raise RuntimeError("TeX engine not found: " + request)
    for name in ("pdflatex", "tectonic"):
        candidate = installed(name)
        if candidate:
            return candidate
    candidate = ROOT / ".build-tools" / "tectonic"
    if candidate.is_file():
        return str(candidate), candidate
    raise RuntimeError(
        "No pdflatex/tectonic found. Install a TeX distribution and put its "
        "compiler on PATH, or optionally supply LATEX_ENGINE=/path/to/compiler. "
        "See doc/README.md.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", default=os.environ.get("LATEX_ENGINE", "auto"),
                        help="auto, a compiler name on PATH, or an optional executable path; defaults to LATEX_ENGINE or auto")
    parser.add_argument("--only-cached", action="store_true",
                        help="Tectonic: disable TeX package network retrieval")
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("--remove-pdf", action="store_true")
    args = parser.parse_args()
    if BUILD.is_symlink():
        raise RuntimeError("refuse to use a symlink as manual build directory")
    if args.clean:
        if BUILD.exists():
            shutil.rmtree(str(BUILD))
        if args.remove_pdf:
            (ROOT / "AMPS.pdf").unlink(missing_ok=True)
        print("Removed generated manual products only.")
        return 0
    if args.remove_pdf:
        parser.error("--remove-pdf requires --clean")

    subprocess.run([sys.executable, "-B", str(ROOT / "check_manual.py")],
                   cwd=str(ROOT), check=True)
    compiler_command, engine = choose_engine(args.engine)
    BUILD.mkdir(exist_ok=True)
    env = os.environ.copy()
    # Tectonic bootstraps ordinary TeX packages. Never write a user-wide cache
    # or require a system package installation just to build this document.
    env["XDG_CACHE_HOME"] = str(ROOT / ".tex-cache")
    env["XDG_CONFIG_HOME"] = str(ROOT / ".tex-config")
    if "tectonic" in engine.name.lower():
        argv = [compiler_command, "--keep-logs", "--keep-intermediates",
                "--outdir", BUILD, "--makefile-rules", BUILD / "dependencies.mk"]
        if args.only_cached:
            argv.append("--only-cached")
        run(argv + ["AMPS.tex"], env, BUILD / "compiler.log")
    else:
        # Three passes settle the contents, tables, citations and cross refs.
        # Bibliography is local TeX, so no BibTeX/style download is necessary.
        for number in range(1, 4):
            run([compiler_command, "-interaction=nonstopmode", "-halt-on-error",
                 "-file-line-error", "-recorder",
                 "-output-directory=" + str(BUILD), "AMPS.tex"],
                env, BUILD / ("compiler-pass-{}.log".format(number)))
    pdf = BUILD / "AMPS.pdf"
    if not pdf.is_file() or not pdf.read_bytes().startswith(b"%PDF-"):
        raise RuntimeError("compiler did not produce a PDF")
    log = (BUILD / "AMPS.log").read_text(errors="replace")
    unresolved = [x for x in log.splitlines()
                  if ("undefined" in x and ("Reference" in x or "Citation" in x))
                  or "There were undefined references" in x
                  or "Rerun to get cross-references right" in x]
    if unresolved:
        raise RuntimeError("unresolved manual references: " + "; ".join(unresolved))
    content = pdf.read_bytes()
    staged = ROOT / "AMPS.pdf.partial"
    staged.write_bytes(content)
    staged.replace(ROOT / "AMPS.pdf")
    receipt = {
        "status": "PASS", "engine_command": compiler_command,
        "engine": str(engine),
        "engine_sha256": hashlib.sha256(engine.read_bytes()).hexdigest(),
        "pdf_sha256": hashlib.sha256(content).hexdigest(),
        "pdf_bytes": len(content),
        "overfull_boxes": log.count("Overfull"),
        "scope": "Documentation compilation only; no AMPS physics runs.",
        "compiler_log_directory": str(BUILD),
    }
    (BUILD / "build-receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print("PASS:", ROOT / "AMPS.pdf", "({} bytes)".format(len(content)))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, subprocess.SubprocessError) as error:
        print("ERROR:", error, file=sys.stderr)
        sys.exit(1)
