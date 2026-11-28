#!/usr/bin/env python3
"""Prepare local measurement--basis data sets from explicitly selected sources.

This is a data-export tool, not another fit. It never selects an optimization
strategy, rescales an absolute rate, adds empirical backgrounds, or resolves
the deferred E3/E5 mass-44 correction. Its selection JSON names the complete
source files and observation blocks to use. Existing fitting outputs remain
unchanged; a new immutable directory is renamed into place only on success.
"""
from __future__ import annotations
import argparse
import copy
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import tempfile

from source_contract import HYBRID, SOLAR, source_ids


def sha(data):
    return hashlib.sha256(data).hexdigest()


def identity(text):
    parts=text.split("--")
    if len(parts)!=2 or any(not re.fullmatch(r"[A-Za-z0-9_](?:[A-Za-z0-9_-]*[A-Za-z0-9_])?",p) for p in parts):
        raise ValueError("case ID must be measurement--basis: "+text)
    source_ids(parts[1])
    return parts


def finite(text, positive=False):
    value=float(text)
    if not math.isfinite(value) or value<0 or (positive and value==0):
        raise ValueError("invalid finite nonnegative source value: "+text)
    if value==0 and Decimal(text)!=0:
        raise ValueError("source value underflows to zero: "+text)
    return value


def read_blocks(text):
    """Strict small receipt reader; no includes or inferred field values.

    Native input validation remains the final authority. This reader exists
    to merge already complete species receipts without reloading reference
    geometry or changing a source parameter after fitting.
    """
    root={"kind":"root","name":"","fields":{},"children":[]};stack=[root]
    for number,line in enumerate(text.splitlines(),1):
        line=line.split("!",1)[0].strip()
        if not line:continue
        if line.startswith("#begin "):
            body=line[7:].split(":",1);kind=body[0].strip();name=body[1].strip() if len(body)==2 else ""
            if not re.fullmatch(r"[a-z][a-z0-9_]*",kind):raise ValueError("invalid block kind")
            if any(c["kind"]==kind and c["name"]==name for c in stack[-1]["children"]):raise ValueError("duplicate block "+name)
            block={"kind":kind,"name":name,"fields":{},"children":[]}
            stack[-1]["children"].append(block);stack.append(block)
        elif line.startswith("#end "):
            if len(stack)==1 or stack[-1]["kind"]!=line[5:].strip():raise ValueError("unbalanced block at line "+str(number))
            stack.pop()
        elif line.startswith("#") or "=" not in line:
            raise ValueError("unsupported directive/line in exported receipt: "+line)
        else:
            key,value=(v.strip() for v in line.split("=",1))
            if key in stack[-1]["fields"]:raise ValueError("duplicate key "+key)
            stack[-1]["fields"][key]=value.strip('"')
    if len(stack)!=1:raise ValueError("unclosed source block")
    return root


def read_source(path):
    data=path.read_bytes();tree=read_blocks(data.decode("utf-8"))
    if tree["fields"] or len(tree["children"])!=1:raise ValueError("expected exactly one source_distribution")
    block=tree["children"][0];f=block["fields"]
    required={"schema_version","source_case","frame","mode","source_law","stripe_surface_mapping","stripe_reference_radius_m"}
    if block["kind"]!="source_distribution" or set(f)!=required:raise ValueError("incomplete/unowned source header")
    if f["schema_version"]!="2" or f["frame"]!="IAU_ENCELADUS" or f["source_law"]!="tenishev2010_point" or f["stripe_surface_mapping"]!="local_vertical_to_reference_sphere":
        raise ValueError("unsupported source schema/frame/law/mapping")
    if f["mode"] not in ("shape","absolute"):raise ValueError("invalid source mode")
    radius=finite(f["stripe_reference_radius_m"],True)
    expected=source_ids(f["source_case"])
    expected_geometry={name.split(":",1)[1]:("point_jet" if name.startswith("jet:") else "polyline_stripe") for name in expected}
    sources={};species=None
    for source in block["children"]:
        sf=source["fields"];name=source["name"]
        if source["kind"]!="source" or name not in expected_geometry or name in sources:raise ValueError("noncanonical source "+name)
        strengths={key[9:] for key in sf if key.startswith("strength.")}
        if not strengths or (species is not None and species!=strengths):raise ValueError("incomplete per-species source vector")
        species=strengths
        for symbol in species:
            if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*",symbol):raise ValueError("invalid species name")
            finite(sf["strength."+symbol])
        fields={k for k in sf if not k.startswith("strength.")}
        geometry=expected_geometry[name]
        legal={"geometry","temperature_K","bulk_speed_m_s"}
        if geometry=="point_jet":legal|={"latitude_deg","west_longitude_deg","tilt_deg","azimuth_deg"}
        if fields!=legal or sf["geometry"]!=geometry:raise ValueError("source geometry fields disagree with canonical basis")
        finite(sf["temperature_K"],True);finite(sf["bulk_speed_m_s"])
        if geometry=="point_jet":
            for key,lo,hi in (("latitude_deg",-90,90),("west_longitude_deg",0,360),("tilt_deg",0,90),("azimuth_deg",0,360)):
                v=float(sf[key])
                if not math.isfinite(v) or not lo<=v<=hi or (hi==360 and v==360):raise ValueError("invalid jet angle")
            if f["source_case"]==SOLAR and (float(sf["tilt_deg"])!=0 or float(sf["azimuth_deg"])!=0):raise ValueError("solar-41 axes must be local-normal")
            if source["children"]:raise ValueError("point source has child blocks")
        else:
            if len(source["children"])<2:raise ValueError("stripe needs vertices")
            previous=None
            for vertex in source["children"]:
                if vertex["kind"]!="vertex" or set(vertex["fields"])!={"x_m","y_m","z_m"} or vertex["children"]:raise ValueError("invalid stripe vertex")
                p=[float(vertex["fields"][v]) for v in ("x_m","y_m","z_m")]
                if not all(math.isfinite(x) for x in p) or abs(math.sqrt(math.fsum(x*x for x in p))-radius)>1e-6:raise ValueError("stripe not on declared surface")
                if p==previous:raise ValueError("zero-length stripe segment")
                previous=p
        sources[name]=source
    if set(sources)!=set(expected_geometry):raise ValueError("missing canonical source")
    for symbol in species:
        total=math.fsum(float(s["fields"]["strength."+symbol]) for s in sources.values())
        if not total>0 or not math.isfinite(total):raise ValueError("invalid source total")
        if f["mode"]=="shape" and abs(total-1)>1e-12:raise ValueError("shape is not normalized")
    return block,species,data


def physical_signature(block):
    result=copy.deepcopy(block);result["name"]=""
    # Compare physical numbers, not irrelevant decimal spelling. Strengths
    # are excluded, while every common gas parameter/vertex remains bound.
    for source in result["children"]:
        source["fields"]={k:(v if k=="geometry" else float(v)) for k,v in source["fields"].items() if not k.startswith("strength.")}
        for vertex in source["children"]:
            vertex["fields"]={k:float(v) for k,v in vertex["fields"].items()}
    result["children"].sort(key=lambda c:c["name"])
    for source in result["children"]:
        # Vertex order is physical and must not be sorted.
        for vertex in source["children"]:vertex["name"]=""
    result["fields"]["stripe_reference_radius_m"]=float(result["fields"]["stripe_reference_radius_m"])
    return result


def render_block(block,depth=0):
    indent="  "*depth
    name=": "+block["name"] if block["name"] else ""
    lines=[indent+"#begin "+block["kind"]+name]
    lines += [indent+key+"="+str(value) for key,value in block["fields"].items()]
    for child in block["children"]:lines.extend(render_block(child,depth+1))
    lines.append(indent+"#end "+block["kind"])
    return lines


def prepare_cases(selection_path,output):
    selection_path=Path(selection_path).resolve();output=Path(output).resolve()
    selection=json.loads(selection_path.read_text())
    if selection.get("schema_version")!=1 or not selection.get("cases") or not selection.get("publication_id"):raise ValueError("selection needs schema_version=1, publication_id and explicit cases")
    if set(selection)!={"schema_version","publication_id","cases"}:raise ValueError("unknown selection fields")
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*",selection["publication_id"]):raise ValueError("invalid publication_id")
    if not isinstance(selection["cases"],dict):raise ValueError("cases must be a mapping of explicit selections")
    if output.exists():raise FileExistsError("refusing existing source set: "+str(output))
    def path(value):
        p=Path(value).expanduser();return p.resolve() if p.is_absolute() else (selection_path.parent/p).resolve()
    output.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".enceladus-case-staging-",dir=output.parent) as temporary:
        staging=Path(temporary)/"data";staging.mkdir()
        catalog={"kind":"case_catalog","name":"enceladus","fields":{"schema_version":"1","publication_id":selection["publication_id"],"frame":"IAU_ENCELADUS"},"children":[]}
        receipt={"schema_version":1,"publication_id":selection["publication_id"],"selection_sha256":sha(selection_path.read_bytes()),"cases":{}}
        radius=None
        for case,settings in selection["cases"].items():
            measurement,basis=identity(case)
            if set(settings)-{"source_files","observation_input","observation_ids","flyby_trajectory_file","physical_CO2_calibration"}:raise ValueError("unknown selection fields for "+case)
            if not settings.get("source_files"):raise ValueError("explicit complete source_files required for "+case)
            merged=None;seen=set();provenance=[]
            for file in settings["source_files"]:
                source_path=path(file);block,species,data=read_source(source_path)
                if block["fields"]["source_case"]!=basis:raise ValueError("case basis differs from source file")
                if seen&species:raise ValueError("multiple alternatives selected for the same species")
                if measurement in ("E3","E5") and "CO2" in species and block["fields"]["mode"]=="absolute" and not settings.get("physical_CO2_calibration"):
                    raise ValueError("E3/E5 physical CO2 is deferred; explicit upstream physical calibration provenance required")
                if merged is None:merged=copy.deepcopy(block)
                else:
                    if physical_signature(merged)!=physical_signature(block):raise ValueError("cannot merge different source geometry, kinematics, laws or modes")
                    by_id={s["name"]:s for s in block["children"]}
                    for source in merged["children"]:
                        source["fields"].update({k:v for k,v in by_id[source["name"]]["fields"].items() if k.startswith("strength.")})
                seen|=species;provenance.append({"source_file":str(source_path),"sha256":sha(data)})
            current_radius=float(merged["fields"]["stripe_reference_radius_m"])
            if radius is not None and radius!=current_radius:raise ValueError("publication cases use different reference spheres")
            radius=current_radius
            merged["name"]=case
            folder=staging/"cases"/case;folder.mkdir(parents=True)
            source_file=folder/"sources.in"
            source_file.write_text("! Complete source data selected upstream; no fitting-method runtime control.\n"+"\n".join(render_block(merged))+"\n")
            provenance_file=folder/"provenance.json"
            provenance_file.write_text(json.dumps({"case_id":case,"upstream_sources":provenance,
                "physical_CO2_calibration":settings.get("physical_CO2_calibration"),"source_selection": "explicit source files; no automatic strategy choice"},indent=2)+"\n")
            observation_input=path(settings["observation_input"])
            obs_tree=read_blocks(observation_input.read_text())
            applications=[b for b in obs_tree["children"] if b["kind"]=="application" and b["name"]=="enceladus"]
            if len(applications)!=1:raise ValueError("observation input requires one Enceladus application")
            app=applications[0];requests=[];instrument=None
            obs_root=app["fields"].get("observation_data_root")
            obs_root=(observation_input.parent/obs_root).resolve() if obs_root else observation_input.parent
            def copy_asset(original,relative):
                target=staging/relative;target.parent.mkdir(parents=True,exist_ok=True)
                data=original.read_bytes();target.write_bytes(data)
                return relative,sha(data)
            for name in settings["observation_ids"]:
                choices=[b for b in app["children"] if b["kind"]=="observation" and b["name"]==name]
                if len(choices)!=1:raise ValueError("unknown observation "+name)
                obs=copy.deepcopy(choices[0]);f=obs["fields"]
                if f.get("event_id")!=measurement or f.get("instrument") not in ("INMS","UVIS"):raise ValueError("observation identity/instrument mismatch")
                if f.get("comparison_species") not in (None,"none") and f["comparison_species"] not in seen:raise ValueError("observation species absent from sources")
                if instrument is not None and instrument!=f["instrument"]:raise ValueError("mixed measurement instruments")
                instrument=f["instrument"]
                for key,hash_key in (("data_file","data_sha256"),("observer_geometry_file","observer_geometry_sha256")):
                    if key not in f:continue
                    value=f[key]
                    origin=observation_input.parent/value[6:] if value.startswith("input:") else obs_root/value
                    relative="observations/"+case+"/"+name+"_"+key+origin.suffix
                    f[key],f[hash_key]=copy_asset(origin,relative)
                requests.append(obs)
            if not requests:raise ValueError("case requires at least one qualified observation")
            f={"measurement_id":measurement,"instrument":instrument,"source_basis":basis,
               "source_distribution_file":str(source_file.relative_to(staging)),"source_sha256":sha(source_file.read_bytes()),
               "available_species":",".join(sorted(seen)),"provenance_file":str(provenance_file.relative_to(staging)),"provenance_sha256":sha(provenance_file.read_bytes())}
            if measurement in ("E3","E5") and "CO2" not in seen:f["unavailable_source.CO2"]="Physical CO2 rate not published; calibration deferred"
            if settings.get("flyby_trajectory_file"):
                f["flyby_trajectory_file"],f["flyby_trajectory_sha256"]=copy_asset(path(settings["flyby_trajectory_file"]),"trajectories/"+case+".csv")
            catalog["children"].append({"kind":"case","name":case,"fields":f,"children":requests})
            receipt["cases"][case]={"species":sorted(seen),"mode":merged["fields"]["mode"],"source_sha256":f["source_sha256"]}
        catalog["fields"]["reference_radius_m"]=format(radius,".17g")
        (staging/"enceladus_cases.in").write_text("\n".join(render_block(catalog))+"\n")
        (staging/"publication.json").write_text(json.dumps(receipt,indent=2)+"\n")
        shutil.copy2(selection_path,staging/"selection.json")
        # Same-filesystem rename makes the complete directory appear at once.
        # No mutable 'current' symlink or replacement of an old case set occurs.
        staging.rename(output)
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--selection",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args()
    try:result=prepare_cases(a.selection,a.output_dir)
    except (ValueError,KeyError,OSError) as error:p.exit(2,"ERROR: "+str(error)+"\n")
    print(json.dumps({"status":"READY","cases":list(result["cases"]),"directory":str(a.output_dir.resolve())},indent=2))
    return 0

if __name__=="__main__":raise SystemExit(main())
