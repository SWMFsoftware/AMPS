"""Standard-library-only synthetic sources for portable native/parser tests.

Canonical identifiers exercise membership, but positions, directions and rates
below are artificial fixtures. They do not reproduce a fitted or observed map.
"""
import json
import math
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"standalone"))
from source_contract import HYBRID,SOLAR,source_ids


def fixture_rates(case):
    return {name:0.0 if i==0 else (i+1)*1e22 for i,name in enumerate(source_ids(case))}


def fixture_config(case,authority):
    radius=252100.0
    config={"schema_version":1,"source_case":case,"frame":"IAU_ENCELADUS",
            "reference_radius_m":radius,"source_law":"tenishev2010_point",
            "molecular_species":"H2O","discrete_jets":[],"tiger_stripes":[]}
    table={}
    for i,label in enumerate(source_ids(case)):
        kind,name=label.split(":",1)
        common={"id":name,"temperature_K":180.0,"bulk_speed_m_s":1500.0}
        if kind=="jet":
            common.update(latitude_deg=-80.0,west_longitude_deg=float(i),
                          tilt_deg=0.0 if case==SOLAR else 30.0,azimuth_deg=0.0)
            config["discrete_jets"].append(common)
        else:
            points=[]
            for phi in (0.01*i,0.01*i+0.001):
                theta=0.1+0.001*i
                points.append([radius*math.sin(theta)*math.cos(phi),radius*math.sin(theta)*math.sin(phi),-radius*math.cos(theta)])
            common.update(geometry_id=name,points_m=points)
            config["tiger_stripes"].append(common);table[name]=points
    if case==HYBRID:
        authority=Path(authority)
        authority.write_text(json.dumps({"frame":"IAU_ENCELADUS","units":"m","stripes":table}))
        config["_stripe_geometry_path"]=str(authority.resolve())
        config["stripe_geometry_file"]=str(authority.resolve())
    return config
