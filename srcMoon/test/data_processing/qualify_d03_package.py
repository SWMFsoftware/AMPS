#!/usr/bin/env python3
"""Reproduce scoped D03 package metadata from evaluated archive/processor QA.

No observation is edited. GCP/PRP science CSV bytes are copied unchanged;
existing I17 history/selection/raw files remain immutable. READY is restricted
to these declared input selections and never implies model/operator readiness.
Publication refuses existing metadata or processed destinations.
"""
import argparse
import json
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from integration_common import sha256


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--processed-dir",type=Path,required=True)
    p.add_argument("--archive-proof",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--data-root",type=Path,default=Path("/data/vtenishe/moon_validation_data"))
    p.add_argument("--publish",action="store_true",help="install only new generated metadata/processed directories; never touch raw")
    a=p.parse_args();out=a.output_dir.resolve();package=a.data_root.resolve()/"diviner"
    if out.exists() and any(out.iterdir()):p.error("output must be new or empty")
    out.mkdir(parents=True,exist_ok=True)
    proof=json.loads(a.archive_proof.read_text())
    if proof["status"]!="PASS" or proof["mismatches"] or proof["missing_official_entries"]:
        raise ValueError("all selected raw files must match authoritative archive checksums")
    sources=proof["sources"]
    for evidence in sources.values():
        if sha256(Path(evidence["path"]))!=evidence["sha256"]:raise ValueError("archive provenance altered")
    if "RDR Version 3.1" not in Path(sources["instrument-data.html"]["path"]).read_text():
        raise ValueError("preserved instrument-team page does not identify calibration version 3.1")
    raw=[]
    for entry in proof["files"]:
        path=Path(entry["path"])
        if sha256(path)!=entry["sha256"]:raise ValueError("raw changed since official verification: "+str(path))
        label=path if path.suffix.lower()==".xml" else path.with_suffix(".xml")
        tree=ET.parse(label)
        lid=tree.findtext(".//{*}Identification_Area/{*}logical_identifier")
        version=tree.findtext(".//{*}Identification_Area/{*}version_id")
        if not lid or not version:raise ValueError("native product identity/version absent: "+str(label))
        raw.append(dict(entry,product_lid=lid,product_label_version=version))
    collections={};processed={};processor_receipts={}
    for subdirectory,key,filename in (("D03-GCP","gcp","d03_gcp_qa.json"),
                                     ("D03-PRP","prp","d03_prp_qa.json")):
        receipt=a.processed_dir/subdirectory/filename
        qa=json.loads(receipt.read_text())
        if qa["status"]!="READY":raise ValueError("processor is not READY: "+str(receipt))
        collection=out/"processed"/key;collection.mkdir(parents=True)
        entries=qa["processed"]
        entries=[entries] if "path" in entries else list(entries.values())
        for product in entries:
            source=Path(product["path"])
            if sha256(source)!=product["sha256"]:raise ValueError("processed file altered: "+str(source))
            target=collection/source.name;shutil.copy2(source,target)
            processed[key+"/"+source.name]={"path":str(package/"processed"/key/source.name),"sha256":product["sha256"]}
        # This is a reproducible metadata relocation only. CSV science bytes,
        # selections, missing values, and their hashes are unchanged.
        def relocate(value):
            if isinstance(value,dict):return {k:(str(package/"processed"/key/Path(v).name) if k=="path" else relocate(v)) for k,v in value.items()}
            if isinstance(value,list):return [relocate(v) for v in value]
            return value
        (collection/filename).write_text(json.dumps(relocate(qa),indent=2,sort_keys=True)+"\n")
        collections[key]={"status":"READY","qa_path":str(package/"processed"/key/filename),
            "native_semantics":qa.get("native_product",qa.get("native_semantics")),"selection":qa["selection"],"counts":qa["counts"]}
        processor_receipts[key]={"path":str(receipt.resolve()),"sha256":sha256(receipt)}
    i17_path=package/"processed/i17/d03_rdr_qa.json";i17=json.loads(i17_path.read_text())
    if i17["status"]!="READY":raise ValueError("I17 processor is not READY")
    for record in (i17["processed"],i17["processor"]):
        if sha256(Path(record["path"]))!=record["sha256"]:raise ValueError("I17 processed/processor identity changed")
    if sha256(Path(i17["contract"]["snapshot_path"]))!=i17["contract"]["sha256"]:raise ValueError("I17 frozen selection changed")
    collections["i17"]={"status":"READY","qa_path":str(i17_path),"counts":i17["counts"],"selection_contract":i17["contract"],"calibration":i17["calibration"]}
    processed["i17/history"]=i17["processed"]
    provenance={"schema":"srcMoon-D03-package-provenance-v1",
      "datasets":{"D03":{"status":"READY","scope":"18 GCP TAB/XML pairs, 2 PRP TAB/XML pairs, frozen June/December I17 selection (87 RDR ZIP/XML pairs)",
        "rdr_calibration_version":"3.1","rdr_calibration_version_source":sources["instrument-data.html"],
        "note":"Calibration release version is distinct from each native PDS4 product-label version."}},
      "raw_files":raw,"processed_files":processed,"collections":collections,
      "archive_proof":{"path":str(a.archive_proof.resolve()),"sha256":sha256(a.archive_proof)},
      "processor_receipts":processor_receipts,
      "processing_script":{"path":str(Path(__file__).resolve()),"sha256":sha256(Path(__file__))},
      "command":sys.argv,
      "processing_scripts":{name:{"path":str(Path(__file__).with_name(name).resolve()),
          "sha256":sha256(Path(__file__).with_name(name))} for name in
          ("prepare_d03_gcp.py","prepare_d03_prp.py","prepare_d03_rdr_eclipse.py","verify_diviner_archive.py")},
      "manual_edits":False,"synthetic_data_used":False,"model_residuals_used":False,
      "limitations":["GCP tbol is bolometric brightness temperature, not a uniquely inferred kinetic surface temperature.",
        "PRP 2-cm annual-average temperature and annual-maximum surface temperature are distinct native quantities.",
        "GCP/PRP per-record uncertainties and RDR cross-channel calibration covariance are not supplied; none was invented.",
        "READY applies to declared data selections only; M1, observation-operator and untouched hold-out model gates remain mandatory."]}
    qa={"schema":"srcMoon-D03-package-QA-v1","overall_status":"READY","dataset_status":{"D03":"READY"},
      "checks":{"official_archive_bytes":True,"frozen_raw_sha256":True,"processor_qa_ready":True,
        "processed_sha256":True,"native_product_ids":True,"frozen_i17_selection":True},
      "raw_file_count":len(raw),"collections":collections,"failures":[],"observational_model_validation":False}
    (out/"provenance.json").write_text(json.dumps(provenance,indent=2,sort_keys=True)+"\n")
    (out/"qa_report.json").write_text(json.dumps(qa,indent=2,sort_keys=True)+"\n")
    (out/"README.md").write_text("# Diviner D03 selected-input package\n\n"
      "READY covers the selected 18 GCP pairs, two PRP pairs and frozen June/December I17 inputs only. "
      "All 214 selected raw files match the preserved official PDS checksum manifests. "
      "The instrument team's release statement identifies archived RDR calibration 3.1; native product-label versions remain separate.\n\n"
      "GCP tbol is bolometric brightness temperature. PRP annual-average temperature at 2 cm depth is distinct from annual-maximum surface temperature. "
      "Missing uncertainties/covariance remain missing. See provenance.json for exact files, native identities, units, selections, processing receipts and limitations.\n\n"
      "Raw files are immutable. Reproduce with AMPS/srcMoon/test/data_processing processors, verify_diviner_archive.py and qualify_d03_package.py. "
      "Processed GCP/PRP CSVs are copied unchanged; metadata paths are relocated by that script. "
      "The original I17 history, quality flags, selection and June/December split are preserved. No model residuals are read.\n\n"
      "Data readiness does not qualify the physical surface state, interpolation/observation operator, I17 thresholds, or a hold-out comparison. "
      "Formal M2 still requires same-campaign I13 and every numerical/physical gate.\n")
    if a.publish:
        destinations=[package/name for name in ("README.md","provenance.json","qa_report.json","processed/gcp","processed/prp")]
        if any(path.exists() for path in destinations):raise FileExistsError("refusing to overwrite existing D03 package metadata/processed products")
        for key in ("gcp","prp"):shutil.copytree(out/"processed"/key,package/"processed"/key)
        for name in ("README.md","provenance.json","qa_report.json"):shutil.copy2(out/name,package/name)
    print(json.dumps({"status":"READY","scope":"data-only","raw_file_count":len(raw),"published":a.publish,"candidate":str(out)},indent=2))
    return 0

if __name__=="__main__":raise SystemExit(main())
