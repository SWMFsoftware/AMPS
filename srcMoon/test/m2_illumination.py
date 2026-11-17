"""I15 geometry staging and independent scoring; never infer D02 epochs.

The production BVH/finite-disc kernel owns geometry. This module owns reviewed
reference identities, exact joins, resolution checks and unchanged scientific
gates. A listed-epoch any-disc occurrence is distinct from mean visible area.
No product, epoch, observer height, uncertainty or approval is synthesized.
"""
from __future__ import annotations

import csv
from datetime import datetime
import hashlib
import io
import json
import math
from pathlib import Path
import re


class Unavailable(RuntimeError):
    """Authoritative scientific input/review is unavailable (SKIPPED)."""


GEOMETRY_FIELDS=("query_id","site_id","latitude_deg","longitude_east_deg",
                 "epoch","time_system","observer_height_m","height_direction","include_earth")


def digest(path):
    value=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):value.update(block)
    return value.hexdigest()


def setting(path,key,default=None):
    """Read an optional flat setting without lowercasing a filesystem path.

    The native shared .in parser remains authoritative. This narrowly scoped
    preflight only handles scalar values and rejects duplicates before it can
    open a scientific contract. Quoted paths preserve spaces and literal '!'.
    """
    import shlex
    found=[]
    for line in Path(path).read_text().splitlines():
        match=re.match(r"^\s*"+re.escape(key)+r"\s*=\s*(.*)$",line)
        if not match:continue
        lexer=shlex.shlex(match[1],posix=True);lexer.whitespace_split=True;lexer.commenters="!"
        values=list(lexer)
        if len(values)!=1:raise ValueError("expected scalar "+key)
        found.append(values[0])
    if len(found)>1:raise ValueError("duplicate "+key)
    return found[0] if found else default


def resolve(owner,value):
    path=Path(value)
    return (path if path.is_absolute() else Path(owner).parent/path).resolve()


def checked_file(owner,record):
    """A missing evidence file is unavailable; altered bytes are corruption."""
    if not isinstance(record,dict) or not record.get("path") or not record.get("sha256"):
        raise ValueError("evidence requires path and SHA-256")
    if not isinstance(record["sha256"],str) or not re.fullmatch("[0-9a-f]{64}",record["sha256"]):
        raise ValueError("invalid illumination evidence SHA-256")
    path=resolve(owner,record["path"])
    if not path.is_file():raise Unavailable("required illumination evidence unavailable: "+str(path))
    if digest(path)!=record["sha256"]:raise ValueError("illumination evidence hash mismatch: "+str(path))
    return path


def json_document(path):
    """Ambiguous duplicate declarations and NaN are corruption, not defaults."""
    def pairs(items):
        out={}
        for key,value in items:
            if key in out:raise ValueError("duplicate illumination JSON key: "+key)
            out[key]=value
        return out
    def nonfinite(value):raise ValueError("non-finite illumination JSON: "+value)
    value=json.loads(Path(path).read_text(),object_pairs_hook=pairs,parse_constant=nonfinite)
    if not isinstance(value,dict):raise ValueError("illumination metadata must be a JSON object")
    return value


def rows(path,required):
    with Path(path).open(newline="") as stream:
        reader=csv.DictReader(stream)
        if (not reader.fieldnames or len(set(reader.fieldnames))!=len(reader.fieldnames) or
            not set(required)<=set(reader.fieldnames)):
            raise ValueError("missing illumination CSV fields: "+str(path))
        result=list(reader)
        if any(None in row or any(value is None for value in row.values()) for row in result):
            raise ValueError("ragged illumination CSV row: "+str(path))
    if not result:raise ValueError("empty illumination evidence: "+str(path))
    return result


def number(value,*,minimum=None,maximum=None):
    if isinstance(value,bool):raise ValueError("Boolean is not a dimensional illumination quantity")
    try:x=float(value)
    except (ValueError,TypeError) as error:raise ValueError("invalid illumination number") from error
    if not math.isfinite(x) or (minimum is not None and x<minimum) or (maximum is not None and x>maximum):
        raise ValueError("invalid finite illumination quantity")
    return x


def binary(value):
    if str(value) not in ("0","1"):raise ValueError("expected exact binary illumination classification")
    return int(value)


def epoch_value(row):
    if row["time_system"]=="TDB":return ("TDB",number(row["epoch"]))
    if row["time_system"]!="UTC":
        raise ValueError("illumination epochs require explicit UTC-Z or TDB seconds")
    match=re.fullmatch(r"(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})(?:\.(\d+))?Z",row["epoch"])
    if not match:raise ValueError("UTC must have explicit calendar/time/Z spelling")
    year,month,day,hour,minute,second=map(int,match.groups()[:6])
    if second>60 or (second==60 and (hour,minute)!=(23,59)):
        raise ValueError("invalid UTC second spelling")
    # Validate the ordinary calendar, but leave a possible leap second for
    # CSPICE's furnished leap-second kernel. No POSIX arithmetic is mistaken
    # for UTC->TDB and .0/.000 spellings cannot create duplicate epochs.
    datetime(year,month,day,hour,minute,min(second,59))
    return ("UTC",year,month,day,hour,minute,second,(match[7] or "").rstrip("0"))


def geometry_bytes(query):
    """Canonical geometry-only bytes; observation columns never reach AMPS."""
    stream=io.StringIO(newline="")
    writer=csv.DictWriter(stream,fieldnames=GEOMETRY_FIELDS,lineterminator="\n",extrasaction="ignore")
    writer.writeheader();writer.writerows(query)
    return stream.getvalue().encode("utf-8")


def native_rows(paths,query,radial_points,azimuth_points):
    """One exact native-geometry join for runtime AND resolution evidence."""
    expected={row["query_id"]:row for row in query};model={}
    required=("query_id","site_id","ephemeris_time_tdb_s","input_epoch","input_time_system",
        "latitude_deg","longitude_east_deg","terrain_visible_area_fraction","earth_visible_area_fraction",
        "jointly_visible_area_fraction","projected_incidence","sampled_any_terrain_visible",
        "sampled_any_jointly_visible","ray_count")
    for path in paths:
        for row in rows(path,required):
            key=row["query_id"]
            if key in model or key not in expected:raise ValueError("duplicate/unexpected native I15 query")
            if row["site_id"]!=expected[key]["site_id"]:raise ValueError("native I15 site identity differs")
            epoch=number(row["ephemeris_time_tdb_s"])
            if epoch_value({"epoch":row["input_epoch"],"time_system":row["input_time_system"]})!=epoch_value(expected[key]):
                raise ValueError("native I15 input epoch identity differs")
            if expected[key]["time_system"]=="TDB" and epoch!=number(expected[key]["epoch"]):
                raise ValueError("native I15 TDB epoch differs")
            for coordinate in ("latitude_deg","longitude_east_deg"):
                if number(row[coordinate])!=number(expected[key][coordinate]):
                    raise ValueError("native I15 coordinate differs")
            area=number(row["jointly_visible_area_fraction"],minimum=0,maximum=1)
            terrain=number(row["terrain_visible_area_fraction"],minimum=0,maximum=1)
            earth=number(row["earth_visible_area_fraction"],minimum=0,maximum=1)
            if area>min(terrain,earth)+1e-13:raise ValueError("joint native mask exceeds one of its component masks")
            count=number(row["ray_count"],minimum=1)
            if int(count)!=count:raise ValueError("native ray count must be an integer")
            # Unit-normal irradiance is <=1. Its native sum contains count
            # positive terms, so roundoff scales with the actual operation
            # count rather than an arbitrary constant independent of grid.
            import sys
            nu=(count+64)*sys.float_info.epsilon
            if nu>=1:raise ValueError("native irradiance summation bound is unresolved")
            number(row["projected_incidence"],minimum=0,maximum=1+nu/(1-nu))
            if binary(row["sampled_any_jointly_visible"])!=int(area>0):
                raise ValueError("native occurrence Boolean disagrees with visible area")
            if binary(row["sampled_any_terrain_visible"])!=int(terrain>0):
                raise ValueError("native terrain occurrence Boolean disagrees with terrain-visible area")
            if number(row["ray_count"])!=radial_points*azimuth_points:
                raise ValueError("native quadrature ray count differs")
            model[key]=row
    if set(model)!=set(expected):raise ValueError("missing native I15 query output")
    return model


def load_contract(path):
    path=Path(path)
    if not path.is_file():raise Unavailable("reviewed I15 contract unavailable")
    contract=json_document(path)
    if contract.get("schema")!="srcMoon-I15-reference-contract-v1":
        raise ValueError("unsupported I15 reference contract schema")
    if contract.get("status")!="READY" or contract.get("frozen_before_model_evaluation") is not True:
        raise Unavailable("I15 exact-product reference contract is not READY/frozen")
    review=contract.get("review",{})
    if review.get("approved") is not True or not review.get("reviewer") or not review.get("date_utc"):
        raise Unavailable("I15 reference contract lacks explicit investigator review")
    if datetime.fromisoformat(review["date_utc"].replace("Z","+00:00")).tzinfo is None:
        raise ValueError("I15 approval must have an explicit time zone")
    if contract.get("selection_role") not in ("calibration","hold-out"):
        raise Unavailable("I15 independent selection role must be explicitly frozen")
    source_ids=contract.get("independent_source_ids")
    if (not isinstance(source_ids,list) or not source_ids or
        any(not isinstance(x,str) or not x for x in source_ids) or len(set(source_ids))!=len(source_ids)):
        raise Unavailable("I15 independent source identities are not frozen")
    operator=contract.get("operator",{})
    required=("frame","time_system","sample_count","epoch_list_endpoint_convention",
              "cadence_description","observer_height_m","height_direction",
              "aberration_correction","earth_occultation","solar_radius_m",
              "solar_radiance_model","time_support_id")
    if any(operator.get(k) is None or operator.get(k)=="" for k in required):
        raise Unavailable("exact-product temporal/observer/ephemeris operator unresolved")
    if operator["frame"]!="MOON_ME_DE421" or operator["solar_radiance_model"]!="uniform":
        raise Unavailable("production finite-disc frame/radiance operator does not match authoritative product")
    if operator["time_system"] not in ("UTC","TDB") or operator["height_direction"] not in ("radial","facet-normal"):
        raise ValueError("unsupported explicit temporal/observer convention")
    if operator["earth_occultation"] not in ("excluded","opaque-sphere"):
        raise Unavailable("D02 Earth-occultation policy is not implemented")
    if operator["aberration_correction"] not in ("NONE","LT","LT+S","CN","CN+S"):
        raise ValueError("unsupported SPICE illumination aberration correction")
    if operator["earth_occultation"]=="opaque-sphere":
        if operator.get("opaque_earth_radius_m") is None:
            raise Unavailable("authoritative opaque Earth radius is unresolved")
        number(operator["opaque_earth_radius_m"],minimum=1e-300)
    number(operator["observer_height_m"],minimum=0);number(operator["solar_radius_m"],minimum=1e-300)
    if type(operator["sample_count"]) is not int or operator["sample_count"]<1:
        raise ValueError("explicit positive per-site sample count required")
    documents=contract.get("authoritative_documents",[])
    if not documents:raise Unavailable("producer timing/geometry documentation unavailable")
    evidence={"contract":{"path":str(path.resolve()),"sha256":digest(path)}}
    for i,record in enumerate(documents):
        checked=checked_file(path,record);evidence["documentation_"+str(i)]={"path":str(checked),"sha256":digest(checked)}
    locators=contract.get("operator_source_locators",{})
    for key in required:
        locator=locators.get(key,{})
        if (type(locator.get("document_index")) is not int or
            not 0<=locator["document_index"]<len(documents) or not locator.get("locator")):
            raise Unavailable("authoritative operator field lacks exact-product source locator: "+key)
    for key in ("geometry_queries","epoch_reference","site_reference","timing_reference","resolution_qualification","spatial_support"):
        if key not in contract:raise Unavailable("I15 requires independent "+key)
        checked=checked_file(path,contract[key]);evidence[key]={"path":str(checked),"sha256":digest(checked)}
    resolution=json_document(evidence["resolution_qualification"]["path"])
    if resolution.get("status")!="PASS" or not resolution.get("angular_and_temporal_resolution_qualified"):
        raise Unavailable("finite-disc occurrence and temporal resolution not qualified")
    # A numerical preflight must identify its real runtime evidence. A lone
    # JSON boolean never substitutes for a resolution campaign.
    if not resolution.get("native_runs"):
        raise Unavailable("illumination resolution qualification lacks native evidence")
    query=rows(evidence["geometry_queries"]["path"],GEOMETRY_FIELDS)
    seen=set();schedules={};locations={}
    for row in query:
        if not re.fullmatch(r"[A-Za-z0-9_.:-]+",row["query_id"]) or row["query_id"] in seen:
            raise ValueError("invalid/duplicate illumination query ID")
        if not re.fullmatch(r"[A-Za-z0-9_.:-]+",row["site_id"]):raise ValueError("invalid site ID")
        seen.add(row["query_id"])
        number(row["latitude_deg"],minimum=-90,maximum=90)
        lon=number(row["longitude_east_deg"],minimum=0,maximum=360)
        if lon==360:raise ValueError("longitude must be in [0,360)")
        epoch=epoch_value(row)
        if (row["time_system"]!=operator["time_system"] or
            number(row["observer_height_m"],minimum=0)!=float(operator["observer_height_m"]) or
            row["height_direction"]!=operator["height_direction"] or
            row["include_earth"]!=("true" if operator["earth_occultation"]=="opaque-sphere" else "false")):
            raise ValueError("query geometry/time conventions disagree with reviewed operator")
        location=(number(row["latitude_deg"]),lon,number(row["observer_height_m"]),row["height_direction"])
        if row["site_id"] in locations and locations[row["site_id"]]!=location:
            raise ValueError("one averaging site has changing native geometry")
        locations[row["site_id"]]=location
        schedules.setdefault(row["site_id"],[]).append(epoch)
    common=None
    for schedule in schedules.values():
        if len(schedule)!=operator["sample_count"] or len(set(schedule))!=len(schedule):
            raise ValueError("incomplete or duplicated per-site averaging epochs")
        current=tuple(sorted(schedule))
        if common is None:common=current
        elif current!=common:raise ValueError("averaging sites have different epoch support")
    return contract,evidence,query


def validate_resolution(owner,record,template,executable,query,*,receipt=None):
    """Require real, hash-bound numerical runs for the selected source/grid.

    Any-disc visibility is discontinuous for a vanishing visible sliver.
    Agreement of two sparse ray grids alone is not an error bound. Therefore
    the review must also name its independent angular/temporal assessment and
    explicitly qualify that sliver/classification error, before any formal
    observational score. This function never generates that scientific proof.
    """
    from build_m2_qualification import production_source_identity
    from integration_common import parse_receipt
    proof_path=checked_file(owner,record);proof=json_document(proof_path)
    context=("source_input_sha256","executable_sha256","production_source_signature_sha256",
             "geometry_queries_sha256","disc_radial_points","disc_azimuth_points")
    if any(proof.get(key) is None for key in context):
        raise Unavailable("illumination resolution proof lacks actual source/input/executable/grid identity")
    if executable is None or not Path(executable).is_file():
        raise Unavailable("actual I15 executable required to validate numerical proof")
    signature=hashlib.sha256(json.dumps(production_source_identity(),sort_keys=True,separators=(",",":")).encode()).hexdigest()
    radial=int(number(setting(template,"illumination_disc_radial_points"),minimum=1))
    azimuth=int(number(setting(template,"illumination_disc_azimuth_points"),minimum=4))
    expected={"source_input_sha256":digest(template),"executable_sha256":digest(executable),
        "production_source_signature_sha256":signature,
        "geometry_queries_sha256":hashlib.sha256(geometry_bytes(query)).hexdigest(),
        "disc_radial_points":radial,"disc_azimuth_points":azimuth}
    if any(proof[key]!=value for key,value in expected.items()):
        raise ValueError("I15 numerical-resolution proof differs from actual source/input/executable/query/grid")
    if receipt and (receipt.get("disc_radial_points")!=radial or receipt.get("disc_azimuth_points")!=azimuth):
        raise ValueError("native I15 disc resolution differs from its reviewed proof")
    assessment=proof.get("independent_resolution_assessment")
    if not assessment:raise Unavailable("independent angular/temporal classification error assessment unavailable")
    assessment_path=checked_file(proof_path,assessment);assessed=json_document(assessment_path)
    if (assessed.get("status")!="PASS" or assessed.get("any_disc_sliver_error_qualified") is not True or
        assessed.get("temporal_event_brackets_qualified") is not True or
        assessed.get("geometry_queries_sha256")!=expected["geometry_queries_sha256"]):
        raise Unavailable("any-disc sliver and temporal error are not independently qualified for these queries")
    if assessed.get("method") not in ("closed-form-bound","independent-numerical-bound") or not assessed.get("evidence"):
        raise Unavailable("resolution assessment lacks an independent bound and its evidence")
    for file in assessed["evidence"]:checked_file(assessment_path,file)
    runs=proof.get("native_runs",[]);grids=set()
    if len(runs)<2:raise Unavailable("I15 numerical proof needs independently retained refinement runs")
    for run in runs:
        native=parse_receipt(checked_file(proof_path,run["receipt"]))
        identity=json_document(checked_file(proof_path,run["manifest"]))
        if (native.get("test_id")!="I15" or identity.get("test_id")!="I15" or
            identity.get("executable",{}).get("sha256")!=expected["executable_sha256"] or
            native.get("mpi_ranks")!=identity.get("mpi_ranks") or identity.get("openmp_threads")!=1):
            raise ValueError("resolution evidence is not the declared actual I15 native build/layout")
        if expected["geometry_queries_sha256"] not in identity.get("inputs",{}).values() and identity.get("illumination_geometry_sha256")!=expected["geometry_queries_sha256"]:
            raise ValueError("resolution native run does not bind the actual geometry queries")
        runtime=[{"path":name,"sha256":value} for name,value in identity.get("inputs",{}).items()
                 if Path(name).suffix==".in"]
        if len(runtime)!=1:raise ValueError("resolution run must identify one actual runtime .in")
        runtime_path=checked_file(proof_path,runtime[0])
        for key in ("surface_geometry","surface_mesh_resolution_m","lola_product_id",
                    "illumination_ray_origin_offset_m","illumination_terrain_solar_model",
                    "illumination_aberration_correction","thermal_solar_disc_radius_m","thermal_opaque_earth_radius_m"):
            default="NONE" if key=="illumination_aberration_correction" else None
            if setting(runtime_path,key,default)!=setting(template,key,default):
                raise ValueError("resolution run physical configuration differs: "+key)
        for key in ("lola_image_file","lola_label_file"):
            if resolve(runtime_path,setting(runtime_path,key))!=resolve(template,setting(template,key)):
                raise ValueError("resolution run uses different terrain data")
        if (int(setting(runtime_path,"illumination_disc_radial_points"))!=native["disc_radial_points"] or
            int(setting(runtime_path,"illumination_disc_azimuth_points"))!=native["disc_azimuth_points"]):
            raise ValueError("resolution native controls disagree with actual input")
        if not run.get("model_outputs"):raise Unavailable("native resolution output tables unavailable")
        model_paths=[checked_file(proof_path,file) for file in run["model_outputs"]]
        if native.get("query_count")!=len(query):raise ValueError("resolution native query count differs")
        native_rows(model_paths,query,int(native["disc_radial_points"]),int(native["disc_azimuth_points"]))
        grids.add((int(native["disc_radial_points"]),int(native["disc_azimuth_points"])))
    if len(grids)<2 or (radial,azimuth) not in grids:
        raise ValueError("resolution proof omits refinement or the actual selected grid")
    return {"path":str(proof_path),"sha256":digest(proof_path),"bound_context":expected,
            "independent_assessment":{"path":str(assessment_path),"sha256":digest(assessment_path)}}


def stage_input(template,destination,*,executable=None):
    """Stage geometry ONLY; all thresholds/observations stay outside AMPS."""
    from integration_common import replace_moon_setting
    raw=setting(template,"m2_formal_contract_file")
    if not raw or raw.lower()=="none":raise Unavailable("I15 needs an explicit reviewed m2_formal_contract_file; D02 metadata remain unavailable")
    from m2_contracts import resolve_contract
    path,selection_evidence=resolve_contract(resolve(template,raw),"I15",unavailable=Unavailable)
    contract,evidence,query=load_contract(path)
    operator=contract["operator"]
    if setting(template,"illumination_terrain_solar_model")!="uniform-disc":
        raise Unavailable("I15 input must explicitly select uniform-disc terrain illumination")
    if (setting(template,"illumination_aberration_correction","NONE")!=operator["aberration_correction"] or
        number(setting(template,"thermal_solar_disc_radius_m"))!=float(operator["solar_radius_m"])):
        raise ValueError("I15 input and reviewed Sun/aberration operator differ")
    if (operator["earth_occultation"]=="opaque-sphere" and
        number(setting(template,"thermal_opaque_earth_radius_m"))!=float(operator["opaque_earth_radius_m"])):
        raise ValueError("I15 opaque-Earth radius differs from authoritative operator")
    resolution=validate_resolution(path,contract["resolution_qualification"],template,executable,query)
    destination=Path(destination);destination.parent.mkdir(parents=True,exist_ok=True)
    geometry=destination.parent/"illumination_geometry_queries.csv"
    geometry.write_bytes(geometry_bytes(query))
    text=Path(template).read_text()
    def set_value(key,value):
        nonlocal text
        if setting(template,key) is None:
            if text.count("#end application")!=1:raise ValueError("expected one application end for optional setting")
            text=text.replace("#end application",key+" = "+value+"\n#end application")
        else:text=replace_moon_setting(text,key,value,Path(template))
    for key in ("spice_path","lola_image_file","lola_label_file"):
        set_value(key,json.dumps(str(resolve(template,setting(template,key)))))
    set_value("m2_formal_contract_file",json.dumps(str(path)))
    set_value("m2_illumination_queries_file",json.dumps(str(geometry.resolve())))
    for key,name in (("surface_cea_file","surface.cea"),("surface_tecplot_file","surface.dat"),
                     ("surface_inventory_checkpoint_file","inventory.chk")):
        set_value(key,json.dumps(str((destination.parent/name).resolve())))
    destination.write_text(text)
    proof={"schema":"srcMoon-I15-staged-reference-v1","evidence":evidence,
           "query_geometry_sha256":digest(geometry),"operator":operator,
           "source_input":{"path":str(Path(template).resolve()),"sha256":digest(template)},
           "staged_input":{"path":str(destination.resolve()),"sha256":digest(destination)},
           "geometry_file":{"path":str(geometry.resolve()),"sha256":digest(geometry)},
           "executable":{"path":str(Path(executable).resolve()),"sha256":digest(executable)},
           "resolution":resolution,"contract_selection":selection_evidence}
    (destination.parent/"illumination_reference_manifest.json").write_text(json.dumps(proof,indent=2,sort_keys=True)+"\n")
    return proof


def spatial_support(path,site_refs):
    """Verify complete predeclared regional native-cell support, not probes.

    A sparse list of sites cannot measure PSR areas/boundaries. The support
    proof identifies rectangular native-pixel domains and a scoring rectangle
    with at least one native-cell padding ring. Area is explicitly measured
    in the native polar-stereographic map plane; no physical surface-area
    correction is guessed. Eight-neighbour native boundary cells supply the
    resolution allowance. Interior disagreement cannot be hidden by assigning
    every site an arbitrary pixel-area tolerance.
    """
    support=json_document(path)
    if (support.get("status")!="READY" or support.get("frozen_before_model_evaluation") is not True or
        support.get("scope")!="resolved-regional-native-cell-coverage"):
        raise Unavailable("resolved frozen PSR native-cell coverage unavailable; sparse sites do not qualify areas")
    if support.get("area_measure")!="native-polar-stereographic-map-plane":
        raise Unavailable("authoritative PSR area measure is unresolved")
    width=number(support["native_pixel_width_m"],minimum=1e-300)
    height=number(support["native_pixel_height_m"],minimum=1e-300)
    area=width*height
    if not math.isfinite(area):raise ValueError("non-finite native PSR map cell area")
    native={};sites={}
    for row in site_refs:
        def integer(key):
            value=number(row[key],minimum=1)
            if int(value)!=value:raise ValueError("native PSR indices must be positive integers")
            return int(value)
        key=(row["native_product_id"],integer("native_line"),integer("native_sample"))
        if key in native or row["site_id"] in sites:raise ValueError("duplicate/overlapping native PSR support")
        native[key]=row;sites[row["site_id"]]=key
        if number(row["pixel_area_m2"],minimum=1e-300)!=area:
            raise ValueError("PSR row area differs from authoritative native map-cell measure")
        binary(row["reference_psr"]);binary(row["psr_boundary_ambiguous"])
    domains=support.get("domains",[])
    if not domains:raise Unavailable("native PSR domain bounds unavailable")
    expected=set();core=set()
    for domain in domains:
        def integer(key):
            value=number(domain[key],minimum=1)
            if int(value)!=value:raise ValueError("native domain bounds must be integer indices")
            return int(value)
        lo,hi,first,last=(integer(k) for k in ("line_min","line_max","sample_min","sample_max"))
        cl,ch,cf,ct=(integer(k) for k in ("scoring_line_min","scoring_line_max","scoring_sample_min","scoring_sample_max"))
        if not (lo<cl<=ch<hi and first<cf<=ct<last):
            raise ValueError("PSR scoring rectangle needs a complete native-neighbour padding ring")
        if (hi-lo+1)*(last-first+1)>len(site_refs):
            raise ValueError("declared native PSR domain exceeds supplied complete support")
        for line in range(lo,hi+1):
            for sample in range(first,last+1):
                key=(domain["native_product_id"],line,sample)
                if key in expected:raise ValueError("overlapping native PSR domains")
                expected.add(key)
                if cl<=line<=ch and cf<=sample<=ct:core.add(key)
    if set(native)!=expected:raise ValueError("sparse/incomplete native PSR cells cannot establish area coverage")
    ambiguous=set()
    for product,line,sample in core:
        row=native[(product,line,sample)];reference=binary(row["reference_psr"])
        edge=any(binary(native[(product,line+di,sample+dj)]["reference_psr"])!=reference
                 for di in (-1,0,1) for dj in (-1,0,1) if di or dj)
        if binary(row["psr_boundary_ambiguous"])!=int(edge):
            raise ValueError("native eight-neighbour edge exclusion differs from fixed reference support")
        derived=area if edge else 0.0
        if number(row["psr_resolution_allowance_m2"],minimum=0)!=derived:
            raise ValueError("PSR resolution allowance is not its native-boundary cell area")
        if edge:ambiguous.add(row["site_id"])
    core_sites={native[key]["site_id"] for key in core}
    if not core_sites-ambiguous:raise Unavailable("PSR region has no unambiguous scoring support")
    return {"scope":support["scope"],"area_measure":support["area_measure"],
            "core_sites":core_sites,"ambiguous_sites":ambiguous,
            "boundary_area_allowance_m2":len(ambiguous)*area,"pixel_area_m2":area,
            "timing_tolerances_s":support.get("timing_tolerances_s",{})}


def score(directory,receipt,*,executable=None,render=True):
    """Score native rows with exact identity joins and original I15 gates."""
    directory=Path(directory)
    staged=json_document(directory/"illumination_reference_manifest.json")
    contract,evidence,query=load_contract(staged["evidence"]["contract"]["path"])
    if evidence!=staged["evidence"]:raise ValueError("I15 reference changed after staging")
    template=checked_file(directory/"illumination_reference_manifest.json",staged["source_input"])
    checked_file(directory/"illumination_reference_manifest.json",staged["staged_input"])
    checked_file(directory/"illumination_reference_manifest.json",staged["geometry_file"])
    bound_executable=checked_file(directory/"illumination_reference_manifest.json",staged["executable"])
    if executable is None or Path(executable).resolve()!=bound_executable:
        raise ValueError("I15 scorer actual executable differs from staged native proof")
    validate_resolution(Path(evidence["contract"]["path"]),contract["resolution_qualification"],template,executable,query,receipt=receipt)
    model=native_rows(sorted(directory.glob("receipt.txt.I15.rank*.csv")),query,
                      int(receipt["disc_radial_points"]),int(receipt["disc_azimuth_points"]))
    if len(model)!=int(receipt["query_count"]):
        raise ValueError("missing native I15 query output")
    if receipt.get("body_fixed_frame")!=contract["operator"]["frame"] or receipt.get("aberration_correction")!=contract["operator"]["aberration_correction"]:
        raise ValueError("native I15 frame/aberration receipt differs")
    references=rows(evidence["epoch_reference"]["path"],("query_id","reference_visible","source_id"))
    seen=set();correct=0;comparisons=[]
    for index,row in enumerate(references):
        key=row["query_id"]
        if key not in model or key in seen or row["source_id"] not in contract["independent_source_ids"]:raise ValueError("invalid independent epoch identity")
        seen.add(key);reference=binary(row["reference_visible"]);predicted=binary(model[key]["sampled_any_jointly_visible"])
        correct+=reference==predicted
        comparisons.append({"comparison_id":"I15_epoch_visibility","x_value":index,"x_unit":"frozen-reference-row",
            "reference_value":reference,"model_value":predicted,"y_unit":"dimensionless","uncertainty_value":"","selection":contract["selection_role"]})
    site_refs=rows(evidence["site_reference"]["path"],("site_id","reference_average_occurrence","reference_psr",
        "psr_boundary_ambiguous","pixel_area_m2","psr_resolution_allowance_m2","logical_surface_element","time_support_id",
        "native_product_id","native_line","native_sample","source_id","latitude_deg","longitude_east_deg"))
    spatial=spatial_support(evidence["spatial_support"]["path"],site_refs)
    aggregates=[];disagreement_area=0;site_ids=set();interior_disagreement=0;occurrence_differences=[]
    for row in site_refs:
        site=row["site_id"]
        if site in site_ids:raise ValueError("duplicate PSR site reference")
        site_ids.add(site)
        samples=[value for value in model.values() if value["site_id"]==site]
        if not samples:raise ValueError("missing model PSR site")
        if row["time_support_id"]!=contract["operator"]["time_support_id"]:raise ValueError("PSR temporal support differs")
        occurrence=sum(binary(value["sampled_any_jointly_visible"]) for value in samples)/len(samples)
        # Geometric terrain shadow is independent of an intervening Earth.
        # Keep the benchmark occurrence above under its reviewed Earth policy,
        # while a separate terrain-only occurrence supplies geometric_psr.
        # A finite declared sample history is not an infinite-time darkness
        # proof, nor does its geometric flag imply kinetic cold trapping.
        terrain_occurrence=sum(binary(value["sampled_any_terrain_visible"])
                               for value in samples)/len(samples)
        if row["source_id"] not in contract["independent_source_ids"]:raise ValueError("unknown independent native PSR source")
        for coordinate in ("latitude_deg","longitude_east_deg"):
            if any(number(value[coordinate])!=number(row[coordinate]) for value in samples):
                raise ValueError("native pixel/reference/query geometry differs")
        if len({number(value["ephemeris_time_tdb_s"]) for value in samples})!=len(samples):
            raise ValueError("duplicate actual native TDB epochs in a site's averaging support")
        psr=int(terrain_occurrence==0);ref=binary(row["reference_psr"])
        area=number(row["pixel_area_m2"],minimum=1e-300)
        if site in spatial["core_sites"] and psr!=ref:
            disagreement_area+=area
            if site not in spatial["ambiguous_sites"]:interior_disagreement+=1
        aggregates.append({"site_id":site,"time_support_id":row["time_support_id"],
            "logical_surface_element":row["logical_surface_element"],"frame":contract["operator"]["frame"],
            "any_disc_occurrence_fraction":occurrence,
            "terrain_any_disc_occurrence_fraction":terrain_occurrence,
            "joint_any_disc_occurrence_fraction":occurrence,
            "earth_occultation_policy":contract["operator"]["earth_occultation"],
            "mean_terrain_visible_area_fraction":sum(float(x["terrain_visible_area_fraction"]) for x in samples)/len(samples),
            "mean_visible_area_fraction":sum(float(x["jointly_visible_area_fraction"]) for x in samples)/len(samples),
            "mean_projected_incidence":sum(float(x["projected_incidence"]) for x in samples)/len(samples),
            "geometric_psr":psr,"geometric_psr_scope":"terrain-dark-over-declared-sampled-epochs",
            "epoch_count":len(samples),
            "native_product_id":row["native_product_id"],"native_line":row["native_line"],"native_sample":row["native_sample"],
            "spatial_selection":("boundary-excluded" if site in spatial["ambiguous_sites"] else
                "scoring-interior" if site in spatial["core_sites"] else "padding-only")})
        average_reference=number(row["reference_average_occurrence"],minimum=0,maximum=1)
        occurrence_differences.append(occurrence-average_reference)
        comparisons.append({"comparison_id":"I15_average_occurrence_diagnostic","x_value":len(aggregates)-1,
            "x_unit":"frozen-native-site-row","reference_value":average_reference,"model_value":occurrence,
            "y_unit":"dimensionless","uncertainty_value":"","selection":contract["selection_role"]})
    if site_ids!={r["site_id"] for r in query} or not aggregates:
        raise ValueError("incomplete/unusable native PSR support")
    timing=rows(evidence["timing_reference"]["path"],("event_id","before_query_id","after_query_id",
        "reference_epoch_tdb_s","event_kind","maximum_timing_error_s","source_id"))
    timing_ok=True;events=set();maximum_error=0
    for event in timing:
        if event["event_id"] in events:raise ValueError("duplicate illumination event")
        events.add(event["event_id"])
        before=model[event["before_query_id"]];after=model[event["after_query_id"]]
        if before["site_id"]!=after["site_id"]:raise ValueError("event brackets different sites")
        lo=number(before["ephemeris_time_tdb_s"]);hi=number(after["ephemeris_time_tdb_s"])
        if hi<=lo:raise ValueError("illumination event bracket not ordered")
        ref=number(event["reference_epoch_tdb_s"]);limit=number(event["maximum_timing_error_s"],minimum=0)
        if event["source_id"] not in contract["independent_source_ids"]:
            raise ValueError("unknown independent illumination timing source")
        if event["event_id"] not in spatial["timing_tolerances_s"]:
            raise Unavailable("independently documented frozen event-timing tolerance unavailable")
        if limit!=number(spatial["timing_tolerances_s"][event["event_id"]],minimum=0):
            raise ValueError("timing row tolerance differs from frozen native-support evidence")
        if event["event_kind"] not in ("sunrise","sunset"):raise ValueError("unsupported illumination event")
        flags=(binary(before["sampled_any_jointly_visible"]),binary(after["sampled_any_jointly_visible"]))
        worst=max(abs(lo-ref),abs(hi-ref));maximum_error=max(maximum_error,worst)
        timing_ok &= flags==((0,1) if event["event_kind"]=="sunrise" else (1,0)) and worst<=limit
    agreement=correct/len(references)
    allowance=spatial["boundary_area_allowance_m2"]
    passed=agreement>=.99 and disagreement_area<=allowance and interior_disagreement==0 and timing_ok
    def write(name,data):
        with (directory/name).open("w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(data[0]),lineterminator="\n");writer.writeheader();writer.writerows(data)
    write("comparison.csv",comparisons);write("illumination_model.csv",aggregates)
    if render:
        from integration_common import render_comparison
        render_comparison(directory,"I15 independent epoch visibility")
    return ("PASS" if passed else "FAIL",{"epoch_visibility_agreement":agreement,
        "reference_epoch_count":len(references),"psr_disagreement_area_m2":disagreement_area,
        "psr_resolution_allowance_m2":allowance,"timing_pass":bool(timing_ok),
        "maximum_event_bracket_error_s":maximum_error,"reference_evidence":evidence,
        "psr_interior_mismatch_count":interior_disagreement,"spatial_scope":spatial["scope"],
        "psr_area_measure":spatial["area_measure"],
        "geometric_psr_scope":"terrain-dark-over-declared-sampled-epochs",
        "maximum_average_occurrence_diagnostic_absolute_difference":max(map(abs,occurrence_differences)),
        "model_aggregate":str(directory/"illumination_model.csv")},
        "linked finite-disc illumination, resolved event brackets and native PSR resolution")
