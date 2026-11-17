#!/usr/bin/env python3
"""Extract explicitly reviewed native D02 rectangles without model selection.

This prepares regional spatial support, not an illumination time history or
an observational comparison. Native projected raster-cell area is explicitly
distinct from true lunar surface area. No region, epoch, logical-model cell,
missing native value, uncertainty or investigator approval is invented.
Canonical data are read-only; all outputs must be outside the data tree.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime
from functools import lru_cache
import json
import math
from pathlib import Path
import sys

TEST_ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(TEST_ROOT))
from m2_illumination import (Unavailable,checked_file,json_document,number,
                            spatial_support,digest)
from prepare_d02_reference import (EXPECTED,integer,scalar,read_i16,
                                   inverse_polar_stereographic)
from audit_d02_temporal_contract import label_value

SCHEMA="srcMoon-I15-spatial-extraction-contract-v1"


def _integer(value,label,minimum=1):
    value=number(value,minimum=minimum)
    if int(value)!=value:raise ValueError(label+" must be an integer")
    return int(value)


def _label(path):
    """Interpret only the frozen native spherical polar-stereographic ABI."""
    text=path.read_text(encoding="ascii")
    exact={"SAMPLE_TYPE":"LSB_INTEGER","MAP_PROJECTION_TYPE":"POLAR STEREOGRAPHIC",
           "COORDINATE_SYSTEM_NAME":"MEAN EARTH/POLAR AXIS OF DE421",
           "POSITIVE_LONGITUDE_DIRECTION":"EAST","KEYWORD_LATITUDE_TYPE":"PLANETOCENTRIC",
           "UNIT":"NONE"}
    for key,value in exact.items():
        if label_value(text,key)!=value:raise ValueError("unsupported native D02 "+key)
    if integer(text,"SAMPLE_BITS")!=16:raise ValueError("native D02 samples must be 16-bit")
    if not label_value(text,"MAP_SCALE").endswith("<m/pix>"):
        raise ValueError("native map scale lacks explicit metre/pixel units")
    for key in ("A_AXIS_RADIUS","B_AXIS_RADIUS","C_AXIS_RADIUS"):
        if not label_value(text,key).endswith("<km>"):
            raise ValueError("native projection radius lacks explicit km units")
    radius=scalar(text,"A_AXIS_RADIUS")
    if radius!=scalar(text,"B_AXIS_RADIUS") or radius!=scalar(text,"C_AXIS_RADIUS"):
        raise ValueError("unsupported nonspherical native projection")
    pole=scalar(text,"CENTER_LATITUDE")
    if pole not in (-90,90) or scalar(text,"CENTER_LONGITUDE")!=0:
        raise ValueError("unsupported native stereographic centre")
    result={"product_id":label_value(text,"PRODUCT_ID"),"lines":integer(text,"LINES"),
        "samples":integer(text,"LINE_SAMPLES"),"scale_m":scalar(text,"MAP_SCALE"),
        "radius_m":radius*1000,"north":pole>0,"line_offset":scalar(text,"LINE_PROJECTION_OFFSET"),
        "sample_offset":scalar(text,"SAMPLE_PROJECTION_OFFSET"),"multiplier":scalar(text,"SCALING_FACTOR"),
        "offset":scalar(text,"OFFSET"),"latitude_min":scalar(text,"MINIMUM_LATITUDE"),
        "latitude_max":scalar(text,"MAXIMUM_LATITUDE")}
    if result["lines"]<1 or result["samples"]<1 or not math.isfinite(result["scale_m"]) or result["scale_m"]<=0:
        raise ValueError("invalid native image dimensions or map scale")
    return result


def _reviewed_contract(path):
    if not path.is_file():raise Unavailable("reviewed native-region selection contract unavailable")
    value=json_document(path)
    if value.get("schema")!=SCHEMA:raise ValueError("unsupported I15 native-region extraction schema")
    review=value.get("review",{})
    if (value.get("status")!="APPROVED" or not isinstance(review,dict) or review.get("approved") is not True or
        value.get("frozen_before_model_evaluation") is not True):
        raise Unavailable("native region/selection is not explicitly reviewed and frozen")
    if not review.get("reviewer") or not review.get("approved_utc"):
        raise ValueError("native-region review requires investigator and approval time")
    if datetime.fromisoformat(review["approved_utc"].replace("Z","+00:00")).tzinfo is None:
        raise ValueError("native-region approval needs an explicit time zone")
    if value.get("area_measure")!="native-polar-stereographic-map-plane":
        raise Unavailable("native projected-area convention has not been explicitly approved")
    if value.get("selection_model_independent") is not True or not value.get("selection_basis"):
        raise Unavailable("native region selection lacks its model-independent basis")
    if not value.get("time_support_id") or not value.get("time_support_provenance"):
        raise Unavailable("exact-product time-support metadata are unavailable; no epoch is inferred")
    checked_file(path,value["time_support_provenance"])
    if not isinstance(value.get("domains"),list) or not value["domains"]:
        raise Unavailable("no explicitly selected native region")
    return value


def prepare(contract_path,data_root,output_dir,*,expected_hashes=None):
    """Reproduce the scoped table; optional hash map is unit-fixture-only.

    The CLI always uses the committed EXPECTED production hashes. A caller's
    temporary toy bytes can test the native decoder, but cannot qualify the
    canonical package. No AMPS executable/result is accepted by this API.
    """
    contract_path=Path(contract_path).resolve();data_root=Path(data_root).resolve();output=Path(output_dir).resolve()
    if output==data_root or data_root in output.parents:
        raise ValueError("regional artifacts must be outside the canonical validation-data tree")
    if output.exists() and any(output.iterdir()):raise ValueError("refusing to overwrite previous regional evidence")
    output.mkdir(parents=True,exist_ok=True)
    result={"id":"I15-spatial-preparation","status":"SKIPPED","scope":"native-region-processing-only",
        "formal_i15_status":"SKIPPED","observational_scoring":False,"model_residuals_used":False,
        "canonical_data_modified":False,"artifacts":{}}
    try:
        contract=_reviewed_contract(contract_path)
        expected=EXPECTED if expected_hashes is None else expected_hashes
        if expected_hashes is not None:
            result["scope"]="control-fixture-native-decoder-only-not-science"
        package=data_root/"lola";raw=package/"raw/illumination"
        metadata={}
        for name in ("README.md","provenance.json","qa_report.json"):
            path=package/name
            if not path.is_file():raise Unavailable("canonical D02 package metadata unavailable: "+name)
            metadata[name]={"path":str(path),"sha256":digest(path)}
        package_qa=json_document(package/"qa_report.json")
        if package_qa.get("dataset_status",{}).get("D02") not in ("READY","PARTIAL"):
            raise Unavailable("D02 native package is not integrity-qualified for spatial processing")
        mapping={}
        if contract.get("logical_mapping_file"):
            with checked_file(contract_path,contract["logical_mapping_file"]).open(newline="") as stream:
                for row in csv.DictReader(stream):
                    key=(row["native_product_id"],_integer(row["native_line"],"line"),_integer(row["native_sample"],"sample"))
                    if key in mapping:raise ValueError("duplicate native logical-model mapping")
                    mapping[key]=_integer(row["logical_surface_element"],"logical element",0)
        labels={};raw_records={};rows=[];native_seen=set();domains=[];area_scale=None;uninterpreted=[]
        for region in contract["domains"]:
            if not region.get("source_id"):raise Unavailable("native product pair needs an explicit source identity")
            psr_id=region.get("native_product_id");visibility_id=region.get("visibility_product_id")
            if not psr_id or not visibility_id:raise Unavailable("both exact native product IDs are required")
            if not psr_id.startswith("LPSR_") or not visibility_id.startswith("AVGVISIB_"):
                raise ValueError("wrong native PSR/visibility product family")
            for product in (psr_id,visibility_id):
                for suffix in (".LBL",".IMG"):
                    name=product+suffix;path=raw/name
                    if name not in expected:raise ValueError("product is not in the frozen native D02 registry: "+name)
                    if not path.is_file():raise Unavailable("required native D02 file unavailable: "+str(path))
                    if name not in raw_records:
                        actual=digest(path)
                        if actual!=expected[name]:raise ValueError("frozen native D02 hash mismatch: "+name)
                        raw_records[name]={"path":str(path),"sha256":actual}
                if product not in labels:
                    labels[product]=_label(raw/(product+".LBL"))
                    if labels[product]["product_id"]!=product:raise ValueError("native label/product identity mismatch")
                    if (raw/(product+".IMG")).stat().st_size!=2*labels[product]["lines"]*labels[product]["samples"]:
                        raise ValueError("native IMG/label byte dimensions disagree")
            psr,vis=labels[psr_id],labels[visibility_id]
            for key in ("lines","samples","scale_m","radius_m","north","line_offset","sample_offset"):
                if psr[key]!=vis[key]:raise ValueError("native PSR/visibility grids do not co-register: "+key)
            if area_scale is None:area_scale=psr["scale_m"]
            elif area_scale!=psr["scale_m"]:raise Unavailable("mixed native map scales require a separately reviewed area contract")
            bounds={key:_integer(region[key],key) for key in
                ("line_min","line_max","sample_min","sample_max","scoring_line_min","scoring_line_max","scoring_sample_min","scoring_sample_max")}
            if not (1<=bounds["line_min"]<bounds["scoring_line_min"]<=bounds["scoring_line_max"]<bounds["line_max"]<=psr["lines"] and
                    1<=bounds["sample_min"]<bounds["scoring_sample_min"]<=bounds["scoring_sample_max"]<bounds["sample_max"]<=psr["samples"]):
                raise ValueError("selected native bounds lack valid complete neighbour padding")
            domains.append(dict(bounds,native_product_id=psr_id))

            @lru_cache(maxsize=None)
            def sample(product,line,sample_index):
                label=labels[product]
                return read_i16(raw/(product+".IMG"),2*((line-1)*label["samples"]+sample_index-1))

            for line in range(bounds["line_min"],bounds["line_max"]+1):
                for column in range(bounds["sample_min"],bounds["sample_max"]+1):
                    key=(psr_id,line,column)
                    if key in native_seen:raise ValueError("overlapping native region selection")
                    native_seen.add(key)
                    latitude,longitude=inverse_polar_stereographic(line-1,column-1,psr["north"],psr["scale_m"],
                        psr["line_offset"],psr["sample_offset"],psr["radius_m"])
                    visibility_dn=sample(visibility_id,line,column);psr_dn=sample(psr_id,line,column)
                    visibility=visibility_dn*vis["multiplier"]+vis["offset"];classification=psr_dn*psr["multiplier"]+psr["offset"]
                    valid=(0<=visibility<=1 and classification in (0,1) and psr["latitude_min"]<=latitude<=psr["latitude_max"])
                    ambiguous=False
                    for di in (-1,0,1):
                        for dj in (-1,0,1):
                            if not (di or dj):continue
                            if not (1<=line+di<=psr["lines"] and 1<=column+dj<=psr["samples"]):
                                ambiguous=True;continue
                            neighbour=sample(psr_id,line+di,column+dj)*psr["multiplier"]+psr["offset"]
                            if neighbour not in (0,1):valid=False
                            if neighbour!=classification:ambiguous=True
                    site=f"{psr_id}:{line}:{column}"
                    if not valid:uninterpreted.append(site)
                    rows.append({"site_id":site,"reference_average_occurrence":format(visibility,".17g") if valid else "",
                        "reference_psr":int(classification) if valid else "","psr_boundary_ambiguous":int(ambiguous),
                        "pixel_area_m2":format(psr["scale_m"]**2,".17g"),
                        "psr_resolution_allowance_m2":format(psr["scale_m"]**2 if ambiguous else 0,".17g"),
                        "logical_surface_element":mapping.get(key,""),"time_support_id":contract["time_support_id"],
                        "native_product_id":psr_id,"native_line":line,"native_sample":column,"source_id":region["source_id"],
                        "latitude_deg":format(latitude,".17g"),"longitude_east_deg":format(longitude,".17g"),
                        "visibility_product_id":visibility_id,"native_visibility_dn":visibility_dn,"native_psr_dn":psr_dn,
                        "quality":"NATIVE_VALUE_UNINTERPRETED" if not valid else "GOOD",
                        "area_measure":"native-polar-stereographic-map-plane"})
        table=output/"site_reference.csv"
        with table.open("w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator="\n");writer.writeheader();writer.writerows(rows)
        support={"schema":"srcMoon-I15-native-regional-support-v1","status":"PARTIAL" if uninterpreted else "READY",
            "frozen_before_model_evaluation":True,"scope":"resolved-regional-native-cell-coverage",
            "area_measure":"native-polar-stereographic-map-plane","native_pixel_width_m":area_scale,
            "native_pixel_height_m":area_scale,"domains":domains,
            "timing_tolerances_s":contract.get("timing_tolerances_s",{}),"time_support_id":contract["time_support_id"]}
        support_path=output/"spatial_support.json";support_path.write_text(json.dumps(support,indent=2,sort_keys=True)+"\n")
        if not uninterpreted:spatial_support(support_path,rows)
        provenance={"schema":"srcMoon-I15-spatial-provenance-v1","processor":{"path":str(Path(__file__).resolve()),"sha256":digest(__file__)},
            "contract":{"path":str(contract_path),"sha256":digest(contract_path)},"native_files":raw_records,
            "package_metadata":metadata,"projection_parameters":labels,"raw_values_edited":False,
            "model_residuals_used":False,"synthetic_data_used":expected_hashes is not None,
            "control_fixture_not_science":expected_hashes is not None,"unknown_values_filled":False,
            "time_support_provenance":contract["time_support_provenance"],
            "processed":{"site_reference":{"path":str(table),"sha256":digest(table)},
                         "spatial_support":{"path":str(support_path),"sha256":digest(support_path)}},
            "source_locators":{"coordinates":"native MAP_PROJECTION fields / committed inverse_polar_stereographic",
                "visibility":"AVGVISIB native DN * SCALING_FACTOR + OFFSET",
                "psr":"LPSR native DN * SCALING_FACTOR + OFFSET",
                "area":"MAP_SCALE squared; explicitly projected raster-cell area, not true spherical surface area"}}
        (output/"provenance.json").write_text(json.dumps(provenance,indent=2,sort_keys=True)+"\n")
        qa={"schema":"srcMoon-I15-regional-spatial-qa-v1","scope":"native-spatial-processing-only",
            "status":"SKIPPED" if uninterpreted else "PASS","row_count":len(rows),"native_value_gaps":uninterpreted,
            "logical_mapping_status":"SUPPLIED" if mapping else "UNAVAILABLE_NOT_ASSUMED",
            "formal_i15_status":"SKIPPED","dataset_D02_promoted":False,"temporal_or_epoch_truth_generated":False}
        (output/"qa_report.json").write_text(json.dumps(qa,indent=2,sort_keys=True)+"\n")
        result.update(status=qa["status"],reason="native regional extraction complete; formal temporal/operator/predecessor gates remain separate" if not uninterpreted else
            "native values/coverage require authoritative interpretation; gaps preserved",row_count=len(rows))
    except Unavailable as error:result.update(status="SKIPPED",reason=str(error))
    except (OSError,ValueError,TypeError,KeyError,OverflowError) as error:result.update(status="ERROR",reason=str(error))
    result["artifacts"]={p.name:str(p) for p in output.iterdir() if p.is_file()}
    (output/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract",type=Path,required=True)
    parser.add_argument("--data-root",type=Path,default=Path("/data/vtenishe/moon_validation_data"))
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args();result=prepare(args.contract,args.data_root,args.output_dir)
    print(json.dumps({"status":result["status"],"reason":result["reason"]},indent=2))
    return {"PASS":0,"FAIL":1,"ERROR":2,"SKIPPED":77}[result["status"]]


if __name__=="__main__":raise SystemExit(main())
