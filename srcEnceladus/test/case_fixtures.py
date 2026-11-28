"""Synthetic case data with independent nonuniform rates and simple geometry.

No fixture is a Cassini reduction or an inferred source map. They exercise
catalog, parser, domain/refinement and native transport interfaces only.
"""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/"standalone"))
from source_case_fixtures import fixture_config,fixture_rates
from fit_source_export import export_vectors
from prepare_forward_cases import prepare_cases
from source_contract import CASES


def make_fixture(root,all_cases=True,mode="absolute",species_pair=False):
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    raw=root/"raw";raw.mkdir()
    blocks=[];source_paths={};case_settings={}
    observations=root/"observations";observations.mkdir()
    measurements=("E3","E5","E14","E17","E18","uvis_2007_zeta_orionis","uvis_2010_e10_solar") if all_cases else ("E3",)
    for basis in CASES if all_cases else CASES[:1]:
        for species in ("H2O","CO2"):
            folder=raw/(basis+"_"+species)
            index=export_vectors(fixture_config(basis,raw/"synthetic_stripes.json"),[("synthetic","fixture",species,fixture_rates(basis),species)],folder,{"fixture":"synthetic nonuniform source rates"})
            source_paths[basis,species]=folder/index["artifacts"][0]["files"][mode]["source_file"]
    for measurement in measurements:
        uvis=measurement.startswith("uvis_")
        species="H2O" if measurement in ("E3","E5") or uvis else "CO2"
        name="obs_"+measurement
        fields={"instrument":"UVIS" if uvis else "INMS","event_id":measurement,"data_file":str(observations/(name+".csv"))}
        if uvis:
            finite=measurement=="uvis_2010_e10_solar"
            fields.update(ray_geometry="finite_observer_to_sun" if finite else "infinite_chord_tangent_track_normal",
                          **{"published_h2o_column_m-2":"1e20","published_h2o_column_sigma_m-2":"1e19"})
            rows=["event_id,time_source_s,time_from_min_altitude_s,tangent_x_m,tangent_y_m,tangent_z_m,optical_depth_raw,poisson_sigma_optical_depth"]
            for t in (-1,0,1):rows.append(f"{measurement},{t},{t},400000,{t*1000},-200000,0.1,0.01")
            if finite:
                auxiliary=observations/(name+"_observer.csv")
                auxiliary.write_text("time_source_s,observer_coordinate_from_tangent_m,ray_ux_table,ray_uy_table,ray_uz_table\n"+
                    "".join(f"{t},-1000000,0,0,1\n" for t in (-1,0,1)))
                fields["observer_geometry_file"]=str(auxiliary)
        else:
            fields.update(time_column="time_s",x_column="x_m",y_column="y_m",z_column="z_m",comparison_species=species,value_column="density_m3",uncertainty_column="sigma_m3")
            rows=["time_s,x_m,y_m,z_m,density_m3,sigma_m3"]
            for t in (-1,0,1):rows.append(f"{t},400000,{t*1000},-200000,1e10,1e9")
        (observations/(name+".csv")).write_text("\n".join(rows)+"\n")
        blocks.append("#begin observation: "+name+"\n"+"\n".join(k+"="+v for k,v in fields.items())+"\n#end observation")
        for basis in CASES if all_cases else CASES[:1]:
            case=measurement+"--"+basis
            active=["H2O","CO2"] if species_pair else [species]
            settings={"source_files":[str(source_paths[basis,s]) for s in active],"observation_input":str(root/"observation_input.in"),"observation_ids":[name]}
            if measurement=="E3":
                trajectory=root/"trajectory.csv"
                trajectory.write_text("flyby,time_from_reference_s,utc,spacecraft_x_m,spacecraft_y_m,spacecraft_z_m,subpoint_x_m,subpoint_y_m,subpoint_z_m,altitude_m,west_lon_planetocentric_deg,lat_planetocentric_deg,west_lon_planetographic_deg,lat_planetographic_deg\n"+
                    "".join(f"E3,{t},synthetic,400000,{t*1000},-200000,0,0,-252100,195113,0,-30,0,-30\n" for t in (-1,0,1)))
                settings["flyby_trajectory_file"]=str(trajectory)
            case_settings[case]=settings
    (root/"observation_input.in").write_text("#begin application: enceladus\n"+"\n".join(blocks)+"\n#end application\n")
    selection=root/"selection.json"
    selection.write_text(json.dumps({"schema_version":1,"publication_id":"synthetic_v1","cases":case_settings},indent=2)+"\n")
    data_root=root/"model_data";publication=data_root/"synthetic_v1"
    prepare_cases(selection,publication)
    inputs=["#begin application: enceladus","schema_version=3","source_set=synthetic_v1",
        "reaction_products=discard","solar_activity=active","plasma_background_mode=none",
        "sample_length_steps=2","total_iterations=4","volume_field_output=false",
        "particle_splitting_merging=off","random_seed=20261113",
        "mesh.domain_min_m=-1260500,-1260500,-1260500","mesh.domain_max_m=1260500,1260500,1260500",
        "mesh.background_cell_size_m=600000","mesh.surface_cell_size_north_m=150000","mesh.surface_cell_size_south_m=150000"]
    for case in case_settings:
        measurement,basis=case.split("--")
        active=["H2O","CO2"] if species_pair else ["H2O" if measurement in ("E3","E5") or measurement.startswith("uvis_") else "CO2"]
        inputs += ["#begin case: "+case,"simulated_species="+",".join(active),"injected_species="+",".join(active),
                   "chemistry_species=none","selected_reactions=none","flyby_density_output=false","flyby_mesh_refinement=false"]
        for s in active:
            inputs += [f"primary_macroparticles_per_step.{s}=10",f"characteristic_speed.{s}_m_s-1=4000",f"courant.{s}=0.5"]
            if mode=="shape":inputs += [f"total_source_rate.{s}_s-1=1e28"]
        inputs += ["#end case"]
    inputs += ["#end application"]
    deck=root/"cases.in";deck.write_text("\n".join(inputs)+"\n")
    return deck,data_root,selection,publication
