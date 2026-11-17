"""Own native subprocess scratch and preserve reviewed input path contracts.

This module contains no model initialization or numerical physics. Native
relative arguments historically used srcSEP as their working directory. Keep
that input base explicit while moving PT, mesh caches and incidental output
into unique, retained, output-owned directories. Never stage a writable link
back into the source tree or search for an alternative input silently.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
from typing import Any, Dict, List, Sequence, Tuple


class NativeExecutionError(RuntimeError):
    """An input cannot be isolated without inventing its path interpretation."""


def create_native_work_directory(output: Path, label: str) -> Path:
    """Allocate one retained directory per process, including discovery.

    A unique leaf prevents two campaigns from sharing fixed PIC names even
    when their caller accidentally reuses an output parent. Report destinations
    remain caller-owned; this function never removes an older run's evidence.
    """
    parent=output.resolve()/"native-work"
    parent.mkdir(parents=True,exist_ok=True)
    safe=re.sub(r"[^A-Za-z0-9_.-]","_",label)
    return Path(tempfile.mkdtemp(prefix=safe+"-",dir=str(parent)))


def _path(value: str, original_root: Path) -> Path:
    """Resolve only the declared old base; no existence-based fallback search."""
    item=Path(value).expanduser()
    return (item if item.is_absolute() else original_root/item).resolve()


def _retain(path: Path, work: Path, category: str, records: List[Dict[str,str]]) -> Path:
    """Copy exact input bytes so the native process cannot modify the original.

    The recognized manifest references are flat CSV inputs, not executable
    configuration trees. A referenced asset's own semantics remain with its
    schema; unknown nested formats are rejected by _manifest rather than
    guessing that their paths are relative to the copied filename.
    """
    try:data=path.read_bytes()
    except OSError as error:raise NativeExecutionError("cannot retain native input "+str(path)+": "+str(error)) from error
    sha=hashlib.sha256(data).hexdigest()
    target=work/"inputs"/(category+"-"+sha+"-"+path.name)
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_bytes(data)
    records.append({"original":str(path),"retained":str(target),"sha256":sha,"kind":category})
    return target


def _manifest(path: Path, work: Path, original_root: Path,
              records: List[Dict[str,str]]) -> Path:
    """Retain native manifests and canonicalize only declared leaf references.

    The flat CV/XM protocols pass filenames literally to C++ and therefore use
    the old process base, srcSEP, rather than the manifest's parent. CLI fault
    inputs contain no filenames. Other schemas need an explicit adapter before
    isolation can preserve nested path meaning; rejecting them is preferable
    to reading another file or changing physical inputs unnoticed.
    """
    raw=_retain(path,work,"manifest",records)
    try:lines=raw.read_text(encoding="utf-8").splitlines()
    except UnicodeError as error:raise NativeExecutionError("native test-input is not UTF-8: "+str(path)) from error
    if not lines:raise NativeExecutionError("native test-input manifest is empty: "+str(path))
    if lines[0]=="srcsep-cli-fault-v1":return raw
    if lines[0]!="srcsep-cv01-native-args-v1" and not lines[0].startswith("srcsep-controlled-native-args-v1:"):
        raise NativeExecutionError("native test-input schema has no declared nested-path isolation contract: "+lines[0])
    tokens=lines[1:]
    if not tokens or len(tokens)%2 or any(not token for token in tokens):
        raise NativeExecutionError("native argument manifest requires complete nonempty name/value pairs")
    references={"--initial-particles","--source-history-csv"}
    for index in range(0,len(tokens),2):
        if tokens[index] in references:
            retained=_retain(_path(tokens[index+1],original_root),work,"asset",records)
            tokens[index+1]=str(retained)
    # Multiple --test-input occurrences preserve the parser's last-wins
    # contract. Retain each distinct manifest separately so earlier provenance
    # never names a file subsequently overwritten by a later argument.
    normalized=work/"inputs"/("resolved-"+hashlib.sha256(raw.read_bytes()).hexdigest()+".txt")
    normalized.write_text("\n".join([lines[0]]+tokens)+"\n",encoding="utf-8")
    records.append({"original":str(raw),"retained":str(normalized),
                    "sha256":hashlib.sha256(normalized.read_bytes()).hexdigest(),"kind":"path-resolved-manifest"})
    return normalized


def isolate_native_arguments(arguments: Sequence[str], work: Path,
                             original_root: Path) -> Tuple[List[str],Dict[str,Any]]:
    """Resolve forwarded path arguments before changing the process directory.

    Python-owned executable/report/output paths are already absolute. Native
    forwarded paths retain the old srcSEP base, in both value syntaxes. The
    versioned --input and reduced-event parsers own their internal base-directory
    contracts; passing their absolute top-level path preserves that ownership.
    Only the flat --test-input protocols above require leaf-asset rewriting.
    """
    paths={"--test-input","--test-output-dir","--test-output","--test-json","--test-junit",
           "--input","-input","--reduced-shock-event","--initialization-output-dir"}
    result=[];records: List[Dict[str,str]]=[];index=0
    while index<len(arguments):
        token=arguments[index];name,separator,value=token.partition("=")
        if name not in paths:
            result.append(token);index+=1;continue
        if not separator:
            if index+1>=len(arguments):raise NativeExecutionError("native path option requires a value: "+name)
            value=arguments[index+1];index+=1
        if not value:raise NativeExecutionError("native path option requires a nonempty value: "+name)
        absolute=_path(value,original_root)
        resolved=_manifest(absolute,work,original_root,records) if name=="--test-input" else absolute
        result.extend([name,str(resolved)]);index+=1
    provenance={"schema":"srcsep-native-execution-v1","working_directory":str(work),
                "native_relative_path_base":str(original_root.resolve()),"input_artifacts":records,
                "requested_arguments":list(arguments),"arguments":result,
                "output_scope":"PT/mesh/cache/incidental files remain in this unique output-owned directory"}
    # Publish complete path provenance before spawning native code. A process
    # cannot consume partially prepared input and leave an apparently complete
    # provenance file after a failed write or interrupted preparation.
    destination=work/"native-execution.json"
    temporary=destination.with_suffix(".json.tmp")
    with temporary.open("w",encoding="utf-8") as stream:
        json.dump(provenance,stream,indent=2)
        stream.write("\n");stream.flush();os.fsync(stream.fileno())
    os.replace(temporary,destination)
    return result,provenance
