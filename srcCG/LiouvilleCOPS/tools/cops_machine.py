"""Machine-specific lines of a LiouvilleCOPS input file (user request 2026-11-03).

A line may end with ``@ <machine> [<machine> ...]`` (names separated by commas and/or blanks; shell patterns
``*``, ``?``, ``[...]``; case-insensitive), for example

    KernelDirectory = /home/vtenishe/SPICE/Kernels/ROSETTA/esa_2022  @ passion
    KernelDirectory = /home/vtenishev/SPICE/Kernels/ROSETTA/esa_2022 @ voyager*

The ``@`` must follow a blank and is the last ``@`` of the line before its ``#`` comment.

Names of the running machine: the host name (short and as returned) and ``$SLURM_CLUSTER_NAME`` (set inside Slurm
jobs, so a cluster's compute nodes match its cluster name). ``$COPS_MACHINE``, when set, replaces them (to choose a
machine's lines explicitly).

Rule, the same in ``src/Input.cpp``: a key with at least one line for this machine takes its value(s) only from
the lines for this machine, in the order of the file; otherwise from its lines without ``@``. Lines for other
machines are ignored. This holds for repeated keys too (``Kernel``): qualify all lines of a list to change it
per machine.
"""
import fnmatch
import os
import re
import socket

_QUALIFIER = re.compile(r"\s@\s*([^@]*)$")


def machine_names():
    """Names this process answers to (lower case, no duplicates, in order)."""
    forced = os.environ.get("COPS_MACHINE", "").strip()
    if forced:
        return [forced.lower()]
    host = socket.gethostname()
    names = [host.split(".")[0], host, os.environ.get("SLURM_CLUSTER_NAME", "")]
    out = []
    for n in names:
        n = n.strip().lower()
        if n and n not in out:
            out.append(n)
    return out


def split_qualifier(active):
    """(text without the qualifier, [patterns] or None) of a comment-free line."""
    m = _QUALIFIER.search(active)
    if m is None:
        return active, None
    patterns = [p for p in re.split(r"[,\s]+", m.group(1).strip()) if p]
    if not patterns:
        raise ValueError("'@' without a machine name: %s" % active.strip())
    return active[:m.start()].rstrip(), patterns


def matches(patterns, names=None):
    names = machine_names() if names is None else names
    return any(fnmatch.fnmatchcase(n, p.lower()) for p in patterns for n in names)


def resolve_text(text, names=None):
    """The input text with every line that does not apply on this machine commented out.

    Lines for this machine lose their qualifier; the result has no ``@`` lines, so every reader of the plain
    ``Key = Value`` format reads it unchanged. A note at the top names the machine. Returns ``(text, info)``,
    where ``info`` maps each machine-specific key to ``{"machines": [...], "applied": bool}``.
    """
    names = machine_names() if names is None else names
    lines = text.splitlines()
    parsed = []          # (index, key, qualifier patterns or None, matched, active text without qualifier)
    for i, line in enumerate(lines):
        active, mark, comment = line.partition("#")
        if not active.strip() or "=" not in active:
            continue
        body, pats = split_qualifier(active)
        key = body.split("=", 1)[0].strip()
        parsed.append((i, key, pats, pats is not None and matches(pats, names), body, mark + comment))
    specific = {key for _, key, pats, ok, _, _ in parsed if ok}
    info = {}
    for i, key, pats, ok, body, comment in parsed:
        if pats is not None:
            rec = info.setdefault(key, {"machines": [], "applied": False})
            rec["machines"] += [p for p in pats if p not in rec["machines"]]
            rec["applied"] |= ok
        if pats is None and key not in specific:
            continue
        if ok:
            lines[i] = body + ("  " + comment if comment else "") + "   # @ " + " ".join(pats)
        elif pats is None:
            lines[i] = "# Overridden on this machine: " + lines[i]
        else:
            lines[i] = "# Other machine: " + lines[i]
    if info:
        lines.insert(0, "# Machine-specific lines resolved for %s (tools/cops_machine.py)" % ", ".join(names))
    return "\n".join(lines) + ("\n" if text.endswith("\n") else ""), info


def missing_message(key, info, names=None):
    """Error text for a required key that is given only for other machines (None if that is not the case)."""
    rec = info.get(key)
    if rec is None or rec["applied"]:
        return None
    names = machine_names() if names is None else names
    return ("%s is given only for the machines %s; this machine is %s (set COPS_MACHINE to choose one)"
            % (key, ", ".join(rec["machines"]), ", ".join(names)))
