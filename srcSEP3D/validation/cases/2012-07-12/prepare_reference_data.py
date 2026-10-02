#!/usr/bin/env python3
"""Download/normalize public July 2012 CME references; never synthesize a shock track.

Run from the AMPS root:
  python3 srcSEP3D/validation/cases/2012-07-12/prepare_reference_data.py

HAPI records retain native cadence and instrument fill values are masked before
SI conversion. HELCATS elongations remain measured angles: their feature and
geometric mapping have not been established as a radial shock-front observable.
The resulting reference-only bundle is discovered by the validation runners.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
from datetime import datetime, timedelta, timezone
import hashlib
import io
import json
import math
from pathlib import Path
import statistics
import tarfile
import urllib.parse
import urllib.request

VALIDATION_ROOT = Path(__file__).resolve().parents[2]
AU_M = 149597870700.0
START, END = "2012-07-12T00:00:00Z", "2012-07-17T00:00:00Z"
EVENT_IDS = ("HCME_A__20120712_02", "HCME_B__20120712_01")


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+"\n", encoding="utf-8")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_figure(fig, bundle, name):
    """Render fully before replacing a published preview and hashing its bytes."""
    for suffix in ("png", "eps"):
        stream=io.BytesIO()
        fig.savefig(stream,format=suffix,dpi=600,facecolor="white")
        path=bundle/(name+"."+suffix)
        temporary=path.with_name(path.name+".partial")
        temporary.write_bytes(stream.getvalue())
        temporary.replace(path)


def stamp(value):
    return datetime.fromisoformat(value.replace("Z", "+00:00"))


def isotime(value):
    return value.astimezone(timezone.utc).isoformat(timespec="milliseconds").replace("+00:00", "Z")


def hapi_url(endpoint, **parameters):
    return "https://cdaweb.gsfc.nasa.gov/hapi/"+endpoint+"?"+urllib.parse.urlencode(parameters)


def sources():
    """Parameter order follows HAPI metadata; different Time axes use @ IDs."""
    urls = {
        "wind_swe_hapi_info.json": hapi_url("info", id="WI_H1_SWE"),
        "wind_mfi_0_info.json": hapi_url("info", id="WI_H0_MFI@0"),
        "wind_mfi_1_info.json": hapi_url("info", id="WI_H0_MFI@1"),
        "wind_orbit_hapi_info.json": hapi_url("info", id="WI_OR_PRE"),
        "wind_plsp_info.json": hapi_url("info", id="WI_PLSP_3DP"),
        "stereo_a_info.json": hapi_url("info", id="STA_COHO1HR_MERGED_MAG_PLASMA"),
        "stereo_b_info.json": hapi_url("info", id="STB_COHO1HR_MERGED_MAG_PLASMA"),
        "ipshocks_database.csv": "https://zenodo.org/records/19730292/files/ipshocks_database.csv?download=1",
        "cfa_wind_event.html": "https://lweb.cfa.harvard.edu/shocks/wi_data/00525/wi_00525.html",
        "helcats_higeocat.json": "https://www.helcats-fp7.eu/catalogues/data/HCME_WP3_V06.json",
        "helcats_higeocat.html": "https://www.helcats-fp7.eu/catalogues/wp3_cat.html",
        "helcats_te_profiles.tar.gz": "https://www.helcats-fp7.eu/catalogues/data/tracks/HCME_WP3_V06_TE_PROFILES.tar.gz",
    }
    for name, dataset, variables, lo, hi in (
        ("wind_swe_data.json", "WI_H1_SWE", "fit_flag,Proton_V_nonlin,Proton_VX_nonlin,Proton_VY_nonlin,Proton_VZ_nonlin,Proton_W_nonlin,Proton_Np_nonlin,BX,BY,BZ,xgse,ygse,zgse", START, END),
        ("wind_mfi_1min_data.json", "WI_H0_MFI@0", "BF1,BGSE,PGSE", START, END),
        ("wind_mfi_3sec_shock_data.json", "WI_H0_MFI@1", "B3F1,B3GSE", "2012-07-14T15:00:00Z", "2012-07-14T20:00:00Z"),
        ("wind_orbit_data.json", "WI_OR_PRE", "GSE_POS,SUN_VECTOR,HEC_POS", START, END),
        ("wind_plsp_data.json", "WI_PLSP_3DP", "MOM.P.DENSITY,MOM.P.AVGTEMP,MOM.P.VELOCITY,MOM.P.VALID", START, END),
        ("stereo_a_positions.json", "STA_COHO1HR_MERGED_MAG_PLASMA", "radialDistance,heliographicLatitude,heliographicLongitude", START, END),
        ("stereo_b_positions.json", "STB_COHO1HR_MERGED_MAG_PLASMA", "radialDistance,heliographicLatitude,heliographicLongitude", START, END),
    ):
        urls[name] = hapi_url("data", id=dataset, parameters=variables,
                              **{"time.min": lo, "time.max": hi, "format": "json"})
    return urls


def download(raw, refresh=False):
    receipt_path = raw / "download_receipt.json"
    previous = json.loads(receipt_path.read_text()) if receipt_path.exists() else {}
    def one(item):
        name, url = item
        path = raw / name
        if path.exists() and not refresh and (name not in previous or previous[name].get("url")==url):
            receipt = previous.get(name, {})
            if receipt and receipt.get("sha256") != digest(path):
                raise ValueError("cached raw file changed: "+name)
            return name, dict(receipt, url=url, sha256=digest(path), bytes=path.stat().st_size,
                retrieved_utc=receipt.get("retrieved_utc", isotime(datetime.fromtimestamp(path.stat().st_mtime, timezone.utc))),
                retrieval_time_basis=receipt.get("retrieval_time_basis", "local cache file modification time"))
        request = urllib.request.Request(url, headers={"User-Agent": "AMPS-SEP3D-reference-preparer/1"})
        with urllib.request.urlopen(request, timeout=60) as response:
            data = response.read()
            headers = {k: response.headers.get(k) for k in ("Last-Modified", "ETag", "Content-Type")}
            final_url = response.url
        if "cdaweb.gsfc.nasa.gov/hapi/" in url:
            payload=json.loads(data)
            if payload.get("status",{}).get("code")!=1200:
                raise ValueError("HAPI service returned unusable data: "+str(payload.get("status")))
        temporary = path.with_name(path.name+".partial")
        temporary.write_bytes(data); temporary.replace(path)
        return name, dict(url=url, final_url=final_url, sha256=digest(path), bytes=len(data),
                          retrieved_utc=isotime(datetime.now(timezone.utc)),
                          retrieval_time_basis="HTTP retrieval", headers=headers)
    raw.mkdir(parents=True, exist_ok=True)
    result = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for name, receipt in pool.map(one, sources().items()):
            result[name] = receipt
            print("[download] "+name+" bytes="+str(receipt["bytes"]), flush=True)
    write_json(receipt_path, result)
    return result


def hapi_records(path):
    payload = json.loads(path.read_text())
    if payload.get("status", {}).get("code") != 1200 or not payload.get("data"):
        raise ValueError("HAPI download is not a successful data response: "+str(path))
    metadata = {p["name"]: p for p in payload["parameters"]}
    columns = list(metadata)
    if any(len(row) != len(columns) for row in payload["data"]):
        raise ValueError("ragged HAPI data: "+str(path))
    rows = [dict(zip(columns, row)) for row in payload["data"]]
    if any(stamp(b["Time"]) <= stamp(a["Time"]) for a, b in zip(rows, rows[1:])):
        raise ValueError("nonmonotonic HAPI times: "+str(path))
    return rows, metadata


def value(row, metadata, key, scale=1):
    candidate = row[key]
    fill = metadata[key].get("fill")
    fill = float(fill) if fill is not None else None
    def scalar(x):
        x = float(x)
        return None if not math.isfinite(x) or x == fill else x*scale
    return [scalar(x) for x in candidate] if isinstance(candidate, list) else scalar(candidate)


def write_csv(path, rows, columns):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader(); writer.writerows(rows)


def normalize(bundle):
    raw = bundle / "raw"
    # This Zenodo version is immutable: reject an unexpected catalog snapshot
    # before selecting an event or attributing the selection to its DOI.
    if hashlib.md5((raw/"ipshocks_database.csv").read_bytes()).hexdigest() != "7432ac3cd68094b32fad2a0b01e33c28":
        raise ValueError("IPShocks v1 publisher checksum does not match")
    swe, metadata = hapi_records(raw / "wind_swe_data.json")
    proton = []
    for row in swe:
        speed = value(row, metadata, "Proton_V_nonlin", 1000)
        density = value(row, metadata, "Proton_Np_nonlin", 1e6)
        width = value(row, metadata, "Proton_W_nonlin", 1000)
        flag = value(row, metadata, "fit_flag")
        good = flag is not None and flag > 1 and all(x is not None and x > 0 for x in (speed, density, width))
        # Wind's scalar isotropic thermal width is the 1/e Maxwellian speed.
        # T = m_p*w^2/(2*k_B); this is a proton temperature, not total pressure.
        temperature = 1.67262192369e-27*width*width/(2*1.380649e-23) if width is not None and width > 0 else None
        proton.append(dict(time_utc=row["Time"], proton_bulk_speed_m_s=speed if good else None,
                           proton_density_m3=density if good else None, proton_temperature_k=temperature if good else None,
                           fit_flag=flag, quality="good" if good else "bad"))
    write_csv(bundle / "wind_protons.csv", proton, list(proton[0]))
    magnetic = {}
    for name, source, bkey, vector in (("wind_magnetic_1min.csv", "wind_mfi_1min_data.json", "BF1", "BGSE"),
                                      ("wind_magnetic_3sec_shock.csv", "wind_mfi_3sec_shock_data.json", "B3F1", "B3GSE")):
        rows, metadata = hapi_records(raw / source)
        normalized = []
        for row in rows:
            magnitude, components = value(row, metadata, bkey, 1e-9), value(row, metadata, vector, 1e-9)
            good = magnitude is not None and magnitude > 0 and all(v is not None for v in components)
            normalized.append(dict(time_utc=row["Time"], magnitude_t=magnitude if good else None,
                bx_gse_t=components[0] if good else None, by_gse_t=components[1] if good else None,
                bz_gse_t=components[2] if good else None, quality="good" if good else "bad"))
        write_csv(bundle / name, normalized, list(normalized[0]))
        magnetic[name] = normalized
    rows, metadata = hapi_records(raw / "wind_orbit_data.json")
    orbit = []
    for row in rows:
        vector = value(row, metadata, "HEC_POS", 1000)
        if any(x is None for x in vector):
            raise ValueError("missing Wind heliocentric orbit data")
        orbit.append(dict(time_utc=row["Time"], x_hec_m=vector[0], y_hec_m=vector[1], z_hec_m=vector[2],
                          heliocentric_radius_m=math.sqrt(sum(x*x for x in vector))))
    write_csv(bundle / "wind_orbit.csv", orbit, list(orbit[0]))
    with (raw / "ipshocks_database.csv").open() as stream:
        matches = [r for r in csv.DictReader(stream) if r["Year"] == "2012" and r["Month (1-12)"] == "7" and
                   r["Day (1-31)"] == "14" and r["Spacecraft"] == "Wind"]
    if len(matches) != 1:
        raise ValueError("expected exactly one Wind shock on July 14")
    selected = matches[0]
    time = datetime(2012, 7, 14, int(selected["Hour (0-23)"]), int(selected["Minute (0-59)"]),
                    int(selected["Second (0-59)"]), tzinfo=timezone.utc)
    left = next(i-1 for i, r in enumerate(orbit) if stamp(r["time_utc"]) >= time)
    a, b = orbit[left], orbit[left+1]
    f = (time-stamp(a["time_utc"])).total_seconds()/(stamp(b["time_utc"])-stamp(a["time_utc"])).total_seconds()
    position = [a[key] + f*(b[key]-a[key]) for key in ("x_hec_m", "y_hec_m", "z_hec_m")]
    arrival = dict(spacecraft="Wind", feature="shock", role="validation", time_utc=isotime(time),
                   heliocentric_radius_m=math.sqrt(sum(x*x for x in position)), uncertainty_s=20.0,
                   position_hec_m=position,
                   provenance="IPShocks Zenodo v1 Wind event at 17:39:09 UTC; timing uncertainty 20 s from independent CfA Wind event 00525 (17:39:07.5 +/-20 s). Radius from linear Cartesian interpolation of NASA WI_OR_PRE (Predicted Orbit) HEC_POS, not an assumed 1 AU. Definitive Orbit did not provide a usable response for this interval.",
                   timing_uncertainty_provenance="CfA event 00525; adopted analysis uncertainty, not an uncertainty column provided by IPShocks")
    write_json(bundle / "arrival.json", arrival)
    write_json(bundle / "ipshocks_wind_event.json", selected)
    catalog = json.loads((raw / "helcats_higeocat.json").read_text())
    selected_catalog = [dict(zip(catalog["columns"], r)) for r in catalog["data"] if r[0] in EVENT_IDS]
    if len(selected_catalog) != 2:
        raise ValueError("HELCATS event identities not found")
    write_json(bundle / "helcats_event_catalog.json", dict(events=selected_catalog,
        excluded_event="HCME_A__20120712_01 precedes the July 12 afternoon eruption", role="geometric-fit context, not independent shock radius"))
    traces = []
    with tarfile.open(raw / "helcats_te_profiles.tar.gz") as archive:
        for event, pa in zip(EVENT_IDS, (75, 260)):
            member = next(m for m in archive if Path(m.name).name == event+"_PA"+str(pa).zfill(3)+".dat")
            content = archive.extractfile(member).read()
            (raw / Path(member.name).name).write_bytes(content)
            for line in content.decode("ascii").splitlines():
                trace, when, elongation, angle, spacecraft = line.split()
                traces.append(dict(event_id=event, trace_id=int(trace), time_utc=when+"Z", spacecraft=spacecraft,
                    elongation_deg=float(elongation), position_angle_deg=float(angle),
                    feature="CME brightness track; shock identification unresolved", role="diagnostic"))
    traces.sort(key=lambda r:(r["time_utc"], r["event_id"], r["trace_id"]))
    write_csv(bundle / "helcats_time_elongation.csv", traces, list(traces[0]))
    for letter in ("a", "b"):
        rows, metadata = hapi_records(raw / ("stereo_"+letter+"_positions.json"))
        positions = [dict(time_utc=r["Time"], radius_m=value(r, metadata, "radialDistance", AU_M),
                          heliographic_latitude_deg=value(r, metadata, "heliographicLatitude"),
                          heliographic_longitude_deg=value(r, metadata, "heliographicLongitude")) for r in rows]
        write_csv(bundle / ("stereo_"+letter+"_positions.csv"), positions, list(positions[0]))
    # Upstream context describes local plasma at Wind. It is not a uniquely
    # established ambient wind everywhere between the Sun and the spacecraft.
    context = {}
    for name, lo, hi in (("upstream", time-timedelta(hours=2), time-timedelta(minutes=5)),
                         ("downstream", time+timedelta(minutes=5), time+timedelta(hours=2))):
        good = [r for r in proton if r["quality"] == "good" and lo <= stamp(r["time_utc"]) <= hi]
        context[name] = dict(start_utc=isotime(lo), end_utc=isotime(hi), samples=len(good),
            **{key:statistics.median(r[key] for r in good) if good else None for key in
               ("proton_bulk_speed_m_s", "proton_density_m3", "proton_temperature_k")})
    write_json(bundle / "wind_context.json", context)
    # Ground-computed PESA Low ion moments are kept as an independent product.
    # They are not spliced into SWE or promoted into acceptance evidence. VALID
    # and fill checks belong to this product; provenance stays visible per row.
    rows, metadata=hapi_records(raw/"wind_plsp_data.json")
    ions=[]
    for row in rows:
        density=value(row,metadata,"MOM.P.DENSITY",1e6)
        temperature=value(row,metadata,"MOM.P.AVGTEMP",1.602176634e-19/1.380649e-23)
        velocity=value(row,metadata,"MOM.P.VELOCITY",1000)
        valid=value(row,metadata,"MOM.P.VALID")
        good=valid==1 and density is not None and density>0 and temperature is not None and temperature>0 and all(v is not None for v in velocity)
        ions.append(dict(time_utc=row["Time"],ion_bulk_speed_m_s=math.sqrt(sum(v*v for v in velocity)) if good else None,
                         ion_density_m3=density if good else None,ion_temperature_k=temperature if good else None,
                         instrument_valid=valid,quality="good" if good else "bad",product="WI_PLSP_3DP"))
    write_csv(bundle/"wind_pesa_ion_moments.csv",ions,list(ions[0]))
    return proton, magnetic["wind_magnetic_1min.csv"], traces, arrival


def figures(bundle, proton, magnetic, traces, arrival):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as dates
    with plt.rc_context({"font.size":9, "ps.fonttype":42}):
        fig, axes = plt.subplots(4, 1, figsize=(7.2, 7.5), sharex=True, layout="constrained")
        # Retain invalid records as NaN so the line does not bridge fit failures.
        # Insert a break across actual telemetry gaps even if both endpoints are good.
        shown=[]
        for previous,row in zip([None]+proton,proton):
            if previous is not None and (stamp(row["time_utc"])-stamp(previous["time_utc"])).total_seconds()>300:
                shown.append(dict(time_utc=isotime(stamp(previous["time_utc"])+timedelta(seconds=1)),
                                  proton_bulk_speed_m_s=None,proton_density_m3=None,proton_temperature_k=None))
            shown.append(row)
        times = [stamp(r["time_utc"]) for r in shown]
        for ax, key, scale, label in zip(axes[:3], ("proton_bulk_speed_m_s", "proton_density_m3", "proton_temperature_k"),
                                        (1e-3, 1e-6, 1e-5), ("Bulk speed (km s$^{-1}$)", "Proton density (cm$^{-3}$)", "Proton T ($10^5$ K)")):
            ax.plot(times, [r[key]*scale if r[key] is not None else float('nan') for r in shown], color="black", lw=.7)
            ax.set_ylabel(label)
        axes[3].plot([stamp(r["time_utc"]) for r in magnetic],
                     [r["magnitude_t"]*1e9 if r["quality"]=="good" else float('nan') for r in magnetic], color="black", lw=.7)
        axes[3].set_ylabel("|B| (nT)")
        for i, ax in enumerate(axes):
            ax.axvline(stamp(arrival["time_utc"]), color="#D55E00", ls="--", lw=1)
            ax.text(.98,.90,"("+chr(97+i)+")", transform=ax.transAxes, ha="right")
            ax.tick_params(direction="in", top=True, right=True)
        axes[0].set_title("Wind observations - July 2012 CME (no model curves)")
        axes[3].set_xlabel("UTC")
        locator=dates.AutoDateLocator(minticks=3,maxticks=8)
        axes[3].xaxis.set_major_locator(locator);axes[3].xaxis.set_major_formatter(dates.ConciseDateFormatter(locator,tz=timezone.utc))
        save_figure(fig,bundle,"wind-reference")
        plt.close(fig)
        fig, ax=plt.subplots(figsize=(7.2,3.6),layout="constrained")
        for event,color in zip(EVENT_IDS,("#0072B2","#D55E00")):
            for trace in sorted({r["trace_id"] for r in traces if r["event_id"]==event}):
                selected=[r for r in traces if r["event_id"]==event and r["trace_id"]==trace]
                ax.plot([stamp(r["time_utc"]) for r in selected],[r["elongation_deg"] for r in selected],"o-",ms=2,lw=.7,color=color,
                    label="STEREO-"+event[5]+" repeat "+str(trace))
        ax.set_ylabel("Elongation (degrees)");ax.set_xlabel("UTC")
        ax.set_title("HELCATS brightness tracks - not classified as shock radius")
        ax.legend(fontsize=7,framealpha=1)
        locator=dates.AutoDateLocator(minticks=3,maxticks=7)
        ax.xaxis.set_major_locator(locator);ax.xaxis.set_major_formatter(dates.ConciseDateFormatter(locator,tz=timezone.utc))
        save_figure(fig,bundle,"helcats-elongation-reference")
        plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir",type=Path,default=VALIDATION_ROOT/"reference_data/CME3D02")
    parser.add_argument("--raw-dir",type=Path,help="use/download a separate raw cache")
    parser.add_argument("--refresh",action="store_true",help="retrieve raw files again instead of using verified cache")
    parser.add_argument("--no-figures",action="store_true",help="prepare data without optional Matplotlib previews")
    args=parser.parse_args();bundle=args.output_dir.resolve();bundle.mkdir(parents=True,exist_ok=True)
    description=Path(__file__).with_name("README.md")
    if description.is_file():
        (bundle/"README.md").write_text(description.read_text(encoding="utf-8"),encoding="utf-8")
    if args.raw_dir:
        import shutil
        (bundle/"raw").mkdir(exist_ok=True)
        for name in sources():
            source=args.raw_dir/name
            if source.is_file() and source.resolve() != (bundle/"raw"/name).resolve():shutil.copyfile(source,bundle/"raw"/name)
    receipts=download(bundle/"raw",args.refresh)
    proton, magnetic, traces, arrival=normalize(bundle)
    if not args.no_figures:figures(bundle,proton,magnetic,traces,arrival)
    # Every packaged reference is checksum-owned. Large all-event catalogs and
    # track archives stay in the raw snapshot to make event selection reproducible.
    files=[]
    selected_names=set(sources())|{"download_receipt.json",*[e+"_PA"+str(pa).zfill(3)+".dat" for e,pa in zip(EVENT_IDS,(75,260))]}
    # Keep downloaded references separate from subsequently attached native
    # run products; rerunning preparation must not re-label those as observations.
    output_names={"README.md","arrival.json","ipshocks_wind_event.json","helcats_event_catalog.json",
                  "helcats_time_elongation.csv","wind_protons.csv","wind_magnetic_1min.csv",
                  "wind_magnetic_3sec_shock.csv","wind_orbit.csv","wind_context.json",
                  "wind_pesa_ion_moments.csv","stereo_a_positions.csv","stereo_b_positions.csv",
                  "wind-reference.png","wind-reference.eps","helcats-elongation-reference.png",
                  "helcats-elongation-reference.eps"}
    for path in sorted(bundle.rglob('*')):
        if not path.is_file():continue
        if path.parent == bundle/"raw":
            if path.name not in selected_names:continue
        elif path.parent != bundle or path.name not in output_names:continue
        files.append(dict(file=path.relative_to(bundle).as_posix(),sha256=digest(path),bytes=path.stat().st_size))
    write_json(bundle/"manifest.json",dict(schema="srcsep3d-swcme-reference-v1",case_id="CME3D02",
        event_name="2012 July 12 CME: shock arrival at Wind",reference_kind="observations",
        comparison_scope="arrival-only",reference_status="ready-for-arrival-diagnostic",
        native_manifest="native-manifest.json",arrival=arrival,files=files,
        counts=dict(wind_proton_records=len(proton),wind_magnetic_1min_records=len(magnetic),
                    helcats_time_elongation_records=len(traces),wind_proton_good=sum(r['quality']=='good' for r in proton)),
        shock_radius_track_available=False,
        limitations=["Native source-off support and installed-provider history exporter are still required.",
                     "HELCATS brightness tracks are not independently classified radial shock observations; do not use their full-track fit to generate holdout shock radii.",
                     "Wind plasma and IMF are reference context; current Parker-background coupling cannot validate CME plasma/ejecta fields."],
        acknowledgements=["NASA/GSFC Space Physics Data Facility and Wind MFI/SWE teams.",
                          "This paper uses data from the Heliospheric Shock Database, generated and maintained at the University of Helsinki.",
                          "HELCATS consortium and STEREO/SECCHI teams; cite Barnes et al. (2019), doi:10.1007/s11207-019-1444-4.",
                          "CfA Interplanetary Shock Database, Michael L. Stevens and Justin C. Kasper."]))
    for entry in files:
        if digest(bundle/entry['file'])!=entry['sha256']:
            raise ValueError("reference changed during preparation: "+entry['file'])
    print("reference_files_verified="+str(len(files)),flush=True)
    print("reference_bundle="+str(bundle),flush=True)
    print("arrival_utc="+arrival['time_utc']+" spacecraft_radius_au="+str(arrival['heliocentric_radius_m']/AU_M),flush=True)


if __name__=="__main__":main()
