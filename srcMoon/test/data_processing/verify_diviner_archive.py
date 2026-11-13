#!/usr/bin/env python3
"""Verify selected local Diviner files against the official archive manifest.

Download provenance only, never replace local raw observations. The archive
MD5 file is read as evidence of official byte identity; SHA-256 still records
all local run identities. A product-label version is not assumed to be the
RDR calibration version. The instrument-team release statement is preserved
separately from the archive's product/collection labels.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import urllib.request

SOURCES={
 "derived.md5":"https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/urn-nasa-pds-lro_diviner_derived1.md5",
 "derived-bundle.xml":"https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/bundle.xml",
 "archive.md5":"https://pds-geosciences.wustl.edu/lro/lro-l-dlre-4-rdr-v1/lrodlr_1001_260921.md5",
 "bundle.xml":"https://pds-geosciences.wustl.edu/lro/lro-l-dlre-4-rdr-v1/bundle_diviner_calibrated.xml",
 "readme.txt":"https://pds-geosciences.wustl.edu/lro/lro-l-dlre-4-rdr-v1/readme.txt",
 "instrument-data.html":"https://www.diviner.ucla.edu/data"}

def digest(path,algorithm="sha256"):
    value=hashlib.new(algorithm)
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):value.update(block)
    return value.hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--processed-dir",type=Path,required=True)
    p.add_argument("--data-root",type=Path,default=Path("/data/vtenishe/moon_validation_data"))
    p.add_argument("--provenance-cache",type=Path,help="reuse hash-verified provenance from an earlier preserved result")
    a=p.parse_args();out=a.output_dir.resolve()
    if out.exists() and any(out.iterdir()):p.error("output directory must be new or empty")
    out.mkdir(parents=True,exist_ok=True)
    report={"schema":"srcMoon-D03-official-byte-verification-v1","status":"ERROR",
      "raw_modified":False,"sources":{},"files":[],"missing_official_entries":[],"mismatches":[]}
    try:
        cached=json.loads((a.provenance_cache/"result.json").read_text())["sources"] if a.provenance_cache else {}
        for name,url in SOURCES.items():
            target=out/name
            prior=cached.get(name)
            if prior and prior["url"]==url and digest(Path(prior["path"]))==prior["sha256"]:
                shutil.copy2(prior["path"],target)
                report["sources"][name]={"url":url,"path":str(target),"sha256":digest(target)}
                continue
            with urllib.request.urlopen(url,timeout=60) as response,target.open("wb") as stream:
                size=0
                while True:
                    block=response.read(1024*1024)
                    if not block:break
                    size+=len(block)
                    if size>256*1024*1024:raise ValueError("provenance document exceeds 256 MiB safety bound")
                    stream.write(block)
            report["sources"][name]={"url":url,"path":str(target),"sha256":digest(target)}
        official={}
        for manifest_name in ("archive.md5","derived.md5"):
          with (out/manifest_name).open() as stream:
            for line in stream:
                fields=line.strip().split(maxsplit=1)
                if len(fields)==2 and len(fields[0])==32:
                    name=fields[1].lstrip("*.").replace("\\","/").lower()
                    # Only a basename is used after verifying uniqueness;
                    # duplicate basenames with different hashes are ambiguous
                    # evidence and must not be silently matched.
                    key=Path(name).name
                    entry=(fields[0].lower(),name,SOURCES[manifest_name])
                    if key in official and (official[key] is None or official[key][0]!=entry[0]):official[key]=None
                    else:official[key]=entry
        package=a.data_root/"diviner";raw={}
        for filename in ("D03-GCP/d03_gcp_qa.json","D03-PRP/d03_prp_qa.json"):
            qa=json.loads((a.processed_dir/filename).read_text())
            if qa["status"]!="READY":raise ValueError("processor QA is not READY: "+filename)
            for name,h in qa["raw_sha256"].items():raw[str(package/name)]=h
        qa=json.loads((package/"processed/i17/d03_rdr_qa.json").read_text())
        for entry in qa["raw_files"]:raw[entry["path"]]=entry["sha256"]
        for name,h in raw.items():
            path=Path(name);expected=official.get(path.name.lower())
            if expected is None:
                report["missing_official_entries"].append(name);continue
            actual_sha=digest(path)
            if actual_sha!=h:raise ValueError("frozen local SHA-256 mismatch: "+name)
            actual_md5=digest(path,"md5")
            entry={"path":name,"sha256":actual_sha,"official_md5":expected[0],"local_md5":actual_md5,
                   "official_path":expected[1],"source_manifest_url":expected[2]}
            report["files"].append(entry)
            if actual_md5!=expected[0]:report["mismatches"].append(entry)
        # Missing manifest entries mean unavailable provenance, never agreement.
        # Altered bytes with an authoritative comparator are corruption ERROR.
        report["status"]="ERROR" if report["mismatches"] else "SKIPPED" if report["missing_official_entries"] else "PASS"
    except (OSError,ValueError,KeyError) as error:report["error"]=str(error)
    (out/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in report.items() if k!="files"},indent=2))
    return 2 if report["status"]=="ERROR" else 0

if __name__=="__main__":raise SystemExit(main())
