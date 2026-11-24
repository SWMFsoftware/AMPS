"""Application-independent, explicit-input model definitions, not a reaction set.

The registry describes mathematical software and published approximations. It
contains no selected species constants, fitted surface parameters, digitized
observations or permission-restricted source data. A caller must separately
establish source rights, states, domains and an observable-specific error budget.
The publication builder includes it only on an explicit command-line request.
"""
import copy
import csv
import hashlib
import json
from pathlib import Path


TABLE_ID = "amps-reduced-model-definitions-v1"
BINDINGS_TABLE_ID = "amps-model-source-bindings-v1"
D119990_SCOPE_REVIEW_ID = "proton-sodium-ground-capture-evidence-v1"
COLUMNS = ("model_id", "implementation", "evidence_kind", "equations_json",
           "inputs_json", "outputs_json", "assumptions_json", "limitations_json", "citations_json")


def _definition(identifier, implementation, equations, inputs, outputs,
                assumptions, limitations, citations=(), evidence="mathematical identity"):
    # Every entry has the same explicit gate. Describing a reduced closure is
    # not permission to execute it for a species or an observational campaign.
    return dict(id=identifier, implementation=implementation, equations=list(equations),
                inputs=list(inputs), outputs=list(outputs), assumptions=list(assumptions),
                limitations=list(limitations), citations=list(citations), evidence_kind=evidence,
                physical_activation=False, status="DEFINITION_NOT_PHYSICAL_QUALIFICATION")


_DEFINITIONS = (
    _definition("ion-mass-balance-v1", "AMPS::ChemicalDatabase::ReducedModels::IonPhysicalMass",
        ["m_ion = m_neutral - m_electron + I/c^2"],
        ["neutral physical mass [kg]", "electron mass [kg]", "ionization energy I [J]", "speed of light c [m/s]"],
        ["singly ionized physical mass [kg]"],
        ["Same isotope and explicit parent/product internal states; I uses the same energy zero."],
        ["Does not choose NIST/CEA/other constants, establish covariances or supply a state/threshold crosswalk."],
        ["https://physics.nist.gov/cuu/Constants/"]),
    _definition("spontaneous-cascade-v1", "AMPS::ChemicalDatabase::ReducedModels::RadiativeBranches",
        ["tau = 1 / sum_j A_ij", "p_ij = A_ij / sum_j A_ij"],
        ["Complete coherent outgoing Einstein-A set for one upper state [s^-1]"],
        ["lifetime [s]", "outgoing branch probabilities [1]"],
        ["Spontaneous optically thin radiative decay; the caller supplies the complete selected channel set."],
        ["No collisional quenching, radiative transfer or invented A-value errors. Never splice inconsistent atomic sources.",
         "CHIANTI native wavelengths are vacuum angstrom; observed/theoretical level energies remain separate."],
        ["https://www.chiantidatabase.org/", "https://www.chiantidatabase.org/referencing.html"]),
    _definition("opal-energy-marginal-v1", "AMPS::ChemicalDatabase::ReducedModels::OpalSecondaryEnergy",
        ["P = (T-I)/2", "w = B*tan(U*atan(P/B))"],
        ["incident energy T [J]", "state-specific threshold I [J]", "explicit source/species width B [J]", "uniform quantile U [0,1]"],
        ["lower-energy outgoing electron energy w [J]"],
        ["Opal energy marginal; heavy-target, negligible-recoil approximation; identical outgoing electrons ordered by energy."],
        ["Not a joint daughter angle/energy distribution or a recoil-conserving three-body generator.",
         "No universal or Na width is supplied; a total Voronov rate does not determine B or a joint law."],
        ["https://doi.org/10.1063/1.1676707", "https://www.sciencedirect.com/science/article/abs/pii/S0092640X72800044"],
        "published reduced distribution"),
    _definition("two-body-conservative-v1", "AMPS::ChemicalDatabase::ReducedModels::TwoBodyProducts",
        ["V=P/(m1+m2)", "Krel=K_available-|P|^2/(2*(m1+m2))", "p=sqrt(2*mu*Krel)",
         "v1=V+p*n/m1; v2=V-p*n/m2"],
        ["two final masses [kg]", "total incoming momentum P [kg m/s]", "incoming translational/available bath energy [J]",
         "separately supplied net state-dependent Q [J]", "explicit center-of-mass product direction n [1]",
         "explicit arithmetic momentum/energy and direction-normalization bounds"],
        ["two product velocities in caller's common inertial frame [m/s]"],
        ["Nonrelativistic two-body final state; incoming photon/plasma momentum and state-dependent Q accounted for by caller."],
        ["Does not supply a measured differential cross section, scattering direction or product-state selection.",
         "Forward CEX is a declared reduced direction choice; unchanged copied velocities do not conserve energy when Q is nonzero."]),
    _definition("parent-velocity-reservoir-v1", "AMPS::ChemicalDatabase::ReducedModels::InheritParentVelocity",
        ["v_product=v_parent", "DeltaP_bath=P_before-P_after", "DeltaE_bath=E_before-E_after", "DeltaQ_bath=Q_before-Q_after"],
        ["before/after heavy-particle mass [kg] and charge [C]", "actual parent velocity in one inertial frame [m/s]",
         "incoming bath momentum [kg m/s], available energy [J] and charge [C]",
         "separately supplied net Q with a common state-energy zero [J]"],
        ["heavy product velocity [m/s]", "signed reservoir momentum [kg m/s], energy [J], charge [C]"],
        ["Explicit heavy-ion transport closure with a reservoir absorbing unresolved degrees of freedom."],
        ["No electron spectra, angular distributions or electron-heating partition is predicted.",
         "Caller must preserve actual parent time/weight, source/state mapping and quantified observable error bounds.",
         "CEX additionally requires the consumed plasma ion and fast-neutral reservoir; an electron-only ledger is insufficient."],
        ["https://doi.org/10.1029/2018JA025542", "https://doi.org/10.1029/2022JE007422"],
        "explicit reduced transport closure"),
    _definition("polanyi-wigner-fixed-sites-v1", "AMPS::ChemicalDatabase::ReducedModels::ThermalRelease",
        ["k=nu*exp(-E/(k_B*T))", "DeltaN=N*(-expm1(-k*dt))"],
        ["attempt frequency nu [s^-1]", "binding energy E [J]", "Boltzmann constant k_B [J/K]", "kinetic surface temperature T [K]",
         "finite initial inventory N [particles or statistical weight]", "interval dt [s]"],
        ["release rate constant [s^-1]", "integrated finite-inventory loss in same units as N"],
        ["First-order independent fixed sites, constant T and parameters over the integration interval; no readsorption/influx in that interval."],
        ["Caller supplies a source-consistent prefactor/energy pair and any coverage/site distribution.",
         "Not an adsorption law, equilibrium proof, multilayer model or universal regolith residence-time prescription."],
        ["https://doi.org/10.1029/2019JE006147"]),
    _definition("jones-2018-water-scenario-v1", "descriptor-only",
        ["Published water-cycle scenario; consult source for species and site assumptions."],
        ["Explicit source-consistent sticking, desorption prefactor, binding energy, illumination and photochemical branch model"],
        ["Scenario definition only; no selected runtime parameters"],
        ["The published molecular-water component uses unity sticking, a 1e13 s^-1 prefactor and 55 kJ/mol activation energy; these are source model choices, not measured universal constants.",
         "Published isotropic ground-state photofragments and unity photodissociation yield are explicit source assumptions."],
        ["No universal H2O/OH surface law or activated daughter law follows from this descriptor.",
         "The full source also contains coupled first/second-order inventories and coverage-dependent OH recombination; the generic fixed-site helper does not implement that full model.",
         "Keep the source's prefactor/energy pair intact; do not combine energies from another inversion with that prefactor."],
        ["https://doi.org/10.1029/2018GL080008"], "published scenario"),
    _definition("schorghofer-2023-diffusive-sites-v1", "descriptor-only",
        ["Site-energy distribution proportional to exp(-(E-Ep)/W) for E>=Ep; consult the published normalization and diffusion model."],
        ["Source-consistent Ep and W [energy]", "source-consistent prefactor [s^-1]", "explicit diffusion/endmember choice"],
        ["Scenario definition only"],
        ["Molecular H2O, efficient-diffusion/site-distribution endmember; not an Ar or gas-OH prescription.",
         "Illustrative source choices Ep=0.65 eV, W=0.22 eV and nu=1e16 s^-1 are kept together as published assumptions, never runtime defaults."],
        ["No default parameter selection, cross-source prefactor substitution or universal sticking law.",
         "The unbounded exponential tail yields a low-coverage divergence; it does not measure infinitely strong real binding.",
         "Actual finite diffusion, material dependence and site capacity require separate qualification."],
        ["https://doi.org/10.3847/PSJ/acf19b"], "published scenario"),
    _definition("cena-hydrogen-reflection-evidence-v1", "descriptor-only", [],
        ["Incident H plasma conditions; source-specific energy/angle/geometry cuts"],
        ["Measured H reflection constraints in the cited domain"],
        ["CENA observations concern hydrogen; empirical spectral shapes can be nonthermal."],
        ["No extrapolation to Na/He, low-energy unmeasured spectrum, or a universal joint ion-regolith reflection law.",
         "Distribution quartiles and population dispersion are not a measurement standard error."],
        ["https://doi.org/10.1029/2011JE004019", "https://doi.org/10.1029/2011GL049362", "https://doi.org/10.1002/jgra.50337"],
        "observational source descriptor"),
    _definition("state-resolved-water-products-v1", "descriptor-only", [],
        ["Explicit wavelength, precursor state, geometry and measured/model source-column role"],
        ["Native laboratory spectra/state populations, with source conditions preserved"],
        ["Narrowband laboratory measurements are not an integrated solar-spectrum daughter generator."],
        ["New O(1S)+H2 channels do not automatically extend a frozen PHIDRATES reaction family.",
         "Unknown ordinate units, fit normalization, state populations or energy-zero mappings remain unresolved.",
         "JWST spectra and ARTEMIS pickup products are external context; they are not independent lunar chemistry validation by themselves."],
        ["https://doi.org/10.1038/s41467-021-22824-7", "https://doi.org/10.1038/s41467-021-26599-9",
         "https://doi.org/10.1038/s41550-024-02203-0", "https://doi.org/10.1029/2022JE007422"],
        "laboratory/observation source descriptor"),
    _definition("effective-threshold-crosswalk-v1", "descriptor-only",
        ["E_photon_effective=I_source+E_mean_source", "K_available=E_photon_effective-I_physical (with full momentum/recoil balance)"],
        ["Explicit original rate-family threshold and units", "source-consistent mean excess energy", "physical state threshold",
         "consistent isotope masses, source spectrum and energy zero"],
        ["Declared mapping between a rate source and physical state energies"],
        ["A legacy/source-effective threshold may differ from a modern state threshold; neither is silently replaced.",
         "The frozen D04 threshold comes from its native THRESHOLD field; its cross-section range endpoint is a different quantity."],
        ["No held ASD numerical threshold or derived ion mass is included in this software definition.",
         "No guessed rounding interval or relaxed source-consistency guard is authorized.",
         "A small threshold shift alone does not bound an absolute-rate change; that requires the actual spectrum/cross-section convolution."],
        ["https://phidrates.space.swri.edu/", "https://www.nist.gov/pml/atomic-spectra-database"],
        "explicit source/state crosswalk specification"),
    _definition("source-mean-isotropic-photon-products-v1", "descriptor-only",
        ["Isotropic product direction in an explicitly declared source/center-of-mass convention; one source-mean excess energy."],
        ["Supported source/state crosswalk", "incoming photon four-momentum in a consistent inertial frame",
         "actual parent velocity/time/weight", "explicit electron and radiation ledgers", "observable-specific approximation error bounds"],
        ["Bounded experimental heavy-particle event specification"],
        ["Fixed-mean/isotropic closure is a reduced model, not a measured differential photoelectron distribution."],
        ["Moving-parent photon transformation and complete signed ledgers remain required.",
         "An unpolarized dipole angular alternative requires an energy/state-dependent anisotropy parameter; no universal He/Na/Ne/Ar value is selected.",
         "Neither a source-mean energy nor isotropy supplies electron heating, spectral widths or state populations."],
        ["https://phidrates.space.swri.edu/"], "explicit experimental closure specification"),
    # Treat the angular law as an explicitly selected mathematical model.
    # The shared two-body kinematics can use its direction, but neither
    # definition chooses a rate, species, internal state or source domain.
    _definition("isotropic-com-two-body-scattering-v1", "descriptor-only",
        ["mu = cos(theta) = 2*U - 1", "phi = 2*pi*V",
         "n = (sqrt(1-mu^2)*cos(phi), sqrt(1-mu^2)*sin(phi), mu)",
         "dP/dOmega = 1/(4*pi) [sr^-1]"],
        ["Two independent unit-uniform variates U and V in [0,1).",
         "Pair center-of-mass axes with theta measured from the incident relative direction.",
         "Separately declared reaction rate, product masses, state-energy Q [J], total incoming momentum [kg m/s] and energy [J]."],
        ["One center-of-mass direction for opposing product momenta.",
         "Compose with two-body-conservative-v1 for the momentum magnitude and event-frame velocities."],
        ["Uniform solid angle in the incoming pair center-of-mass frame.",
         "One shared sampled direction for both products, with opposite momenta and the reaction Q applied."],
        ["Experimental sensitivity model. It supplies no measured angular law, integral rate, state mapping or physical applicability bound.",
         "For forward-peaked collisions the model can substantially increase the predicted recoil energy."],
        evidence="explicit experimental sensitivity model"),
    # These are provenance findings about existing sources. A figure-scale
    # discrepancy is not an extracted numerical angular or rate table.
    _definition(D119990_SCOPE_REVIEW_ID, "descriptor-only", [],
        ["Selected reference reaction H+ + Na(3s) -> H(1s) + Na+(2p6).",
         "Existing source records iaea-d119990-20261124:live-html and chemistry-archive-20261124:knoop-h-na.",
         "Existing integral table collisiondb_119990_embedded_arrays.csv, SHA-256 2a887ad43497a7f9e2b3d7b6173ae1b2597cc893d897926762c181946c92db93."],
        ["Unresolved source/state crosswalk and angular-data gap for the declared reaction."],
        ["The native IAEA label leaves the incoming Na state and the original publication unspecified.",
         "Knoop Chapter 5 section 5.2.3 reports laboratory ground-product angular results from 7 keV/u upward. Figure 5.10 shows a peak-normalized ground-product trace at 25 keV/u. The 4 keV/u panel omits that trace.",
         "At 4 keV/u the exact D119990 value is 6.34e-15 cm^2. Knoop Figure 5.6 puts the prompt outer-shell H(n=1) point between the 1e-17 and 1e-16 cm^2 decades. This is a figure-scale consistency check."],
        ["The published angular trace supplies no numerical absolute bins, full-tail coverage or 4 keV/u ground-product distribution.",
         "Resolve the original citation, initial state and prompt or eventual product convention before combining the rate and angle evidence.",
         "The source review does not reassign the D119990 channel or qualify its rate for physical transport."],
        ["https://db-amdis.org/collisiondb/datasets/119990/",
         "https://pure.rug.nl/ws/portalfiles/portal/2855778/c5.pdf",
         "https://pure.rug.nl/ws/portalfiles/portal/2855777/c4.pdf"],
        "published source-scope review"),
)


def definitions():
    """Return a defensive copy; no registry field is a physical default."""
    return copy.deepcopy(_DEFINITIONS)


def sodium_scenarios(datasets):
    """Attach each scenario only when its reviewed source objects are in this selection.

    Use payload identities, so replaying an acquisition at a different time
    does not invalidate a scenario. A different source edition may still be
    published in its own right. Its model binding needs a new scope review.
    """
    selected = []
    # Resolve within the named source packet. A matching identifier in a
    # differently scoped dataset cannot supply the selector exported below.
    packets = {dataset["id"]: {row["id"]: row for row in dataset["files"]}
               for dataset in datasets if "source_packet" in dataset}
    pickup_files = packets.get("poppe2022-lunar-pickup-rates", {})
    poppe_table = "poppe2022:atomic-1au-rates"
    poppe_hashes = {
        poppe_table + ":csv": "76fb5cdb5181a5a4a6ceb0e9fff40da28fb3069affe107e440644863ed083b61",
        "poppe2022:ionization-rates": "a208a94f318bbe2c63c59e6e28cfd15bf78901b7f7560480c8298bb235b7826e",
        "poppe2022:article": "2c8facf7d53366bb65c9f4285ed0e18e335536527926f9c898a1718747ea1ec0",
    }
    if all(pickup_files.get(key, {}).get("sha256") == digest for key, digest in poppe_hashes.items()):
        projected = pickup_files[poppe_table + ":csv"]
        if projected.get("table", {}).get("id") != poppe_table or projected["evaluation_permission"] != "allowed":
            raise ValueError("reviewed Poppe median reference requires its accessible named CSV")
        pickup = _definition("poppe2022-na-pickup-median-v1", "descriptor-only",
            ["S_j = n_Na * nu_j", "v_birth = (0, 0, 0) in the lunar frame"],
            ["Local neutral Na number density [m^-3]", "Deposited median frequencies per neutral at 1 au [s^-1]"],
            ["Volumetric Na+ source for each selected process [m^-3 s^-1]", "Lunar-frame ion birth velocity [m/s]"],
            ["Select the Na row and the median of each process. Section 2.2 neglects parent-neutral motion at pickup.",
             "Compute each process only where its driver acts. The authors trace pickup in steady uniform crossed fields."],
            ["The proton charge-exchange frequency includes all contributing final states under the published solar wind.",
             "Adding a ground-capture source to that frequency for the same population would double count events.",
             "The minimum and maximum columns describe different conditions; they are not errors on the median."],
            # Poppe supplies the rates and birth condition. Hartle and Shen
            # independently document this birth approximation for pickup-ion
            # transport; neither supplies a Na binary differential cross section.
            ["https://doi.org/10.1029/2022JE007422", "https://doi.org/10.5281/zenodo.6994648",
             "https://doi.org/10.1029/2011JA016859", "https://doi.org/10.3847/1538-4357/ad40a4"],
            "published lunar pickup source model")
        # Select columns by name. Native decimal tokens and the CC BY 4.0
        # attribution stay with the authoritative CSV and its parent deposit.
        pickup["source_bindings"] = dict(
            rate_kind="frequency_per_neutral", unit="s^-1", heliocentric_distance_au=1,
            table_id=poppe_table, row_match={"native_species": "Na"},
            columns={process: process + "_median_per_s" for process in
                     ("photoionization", "solar_wind_proton_charge_exchange", "solar_wind_electron_impact")},
            birth_velocity_m_s=[0, 0, 0], velocity_frame="lunar", source_sha256=poppe_hashes)
        selected.append(pickup)

    theory_id = "shingal1986-na:ground-partial-benchmark"
    theory_sha = "58bc33f7f8e7dd960ff2972edf0bd3999ada697e8ec666f04da97bc08e7e8078"
    theory_packet = "sodium-ground-capture-benchmark-20261124"
    record = packets.get(theory_packet, {}).get(theory_id, {})
    if record.get("sha256") == theory_sha:
        point = record.get("quantity_metadata", {}).get("benchmark_point", {})
        if any(point.get(key) != expected for key, expected in
               (("projectile_energy", "4.0"), ("energy_unit", "keV/u"), ("cross_section_cm2", "2.2E-17"),
                ("energy_frame", "target-rest laboratory"), ("method", "A34 theory"))):
            raise ValueError("reviewed Shingal reference requires its Table 4 benchmark metadata")
        ground = _definition("hplus-na3s-ground-capture-4kev-u-isotropic-v1", "descriptor-only",
            ["nu = n_projectile * sigma_1s * v_relative for a monoenergetic beam",
             "dP/dOmega = 1/(4*pi) in the pair center of mass"],
            ["Published A34 ground-product partial at the selected target-rest energy",
             "Projectile number density [m^-3], relative speed [m/s], state-consistent masses [kg] and reaction Q [J]"],
            ["Single-channel frequency per Na target [s^-1]", "Opposing COM momenta and both event-frame product velocities [m/s]"],
            ["H+ + Na(3s) -> H(1s) + Na+(2p6). Both outgoing atoms/ions are in their ground states.",
             "The AMPS sensitivity scenario supplies isotropy. Shingal Table 4 supplies the integral cross section."],
            ["Use only at the source energy. No off-energy interpolation or solar-wind averaging is defined.",
             "The physical angular distribution and numerical cross-section uncertainty are unspecified.",
             "Hold production fixed when comparing modeled Na+ spectra or return fractions across birth-direction laws."],
            ["https://doi.org/10.1088/0022-3700/19/3/014"], "atomic theory with an approved modeled scattering law")
        # The benchmark is already inspectable publication metadata. Point to
        # that fact without copying restricted PDF bytes or creating a new
        # numerical table from them. Both products share one relative direction.
        ground["source_bindings"] = dict(
            rate_kind="single_energy_partial_cross_section", dataset_id=theory_packet,
            source_id=theory_id, metadata_path=["quantity_metadata", "benchmark_point"],
            area_field="cross_section_cm2", area_unit="cm^2", projectile_energy_kev_u="4.0",
            energy_frame="target-rest laboratory", off_energy_policy="undefined",
            incident_states=["H+", "Na(3s)"], product_states=["H(1s)", "Na+(2p6)"],
            angular_model_id="isotropic-com-two-body-scattering-v1", recoil_model_id="two-body-conservative-v1",
            source_sha256={theory_id: theory_sha})
        selected.append(ground)
    return selected


def materialize(root, release_id, amps_root, datasets=()):
    """Create authored software metadata as ordinary permission-cleared sources.

    This is NOT a download and carries no fabricated HTTP/acquisition record.
    Exact generated JSON/CSV and the repository licence become release members;
    all third-party source-data admission still uses source_packets.py.
    """
    directory = Path(root) / "releases" / release_id / "provenance" / "authored-models"
    directory.mkdir(parents=True, exist_ok=False)
    scenarios = sodium_scenarios(datasets)
    values = definitions() + tuple(scenarios)
    document = directory / "model-definitions.json"
    document.write_text(json.dumps(dict(schema="amps-chemical-model-definitions-v1", models=values,
        physical_activation=False), indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    table = directory / "model-definitions.csv"
    with table.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(COLUMNS)
        for value in values:
            writer.writerow([value["id"], value["implementation"], value["evidence_kind"]] +
                [json.dumps(value[key], ensure_ascii=False, sort_keys=True) for key in
                 ("equations", "inputs", "outputs", "assumptions", "limitations", "citations")])
    licence = Path(amps_root) / "LICENSE.txt"
    if not licence.is_file():
        raise ValueError("AMPS software licence is required for authored model definitions")
    policy = dict(storage_permission="allowed", evaluation_permission="allowed",
        redistribution="Apache-2.0 software; third-party cited data retain their own rights",
        evidence=[dict(url="https://github.com/SWMFsoftware/AMPS/blob/master/LICENSE.txt")],
        conditions="Mathematical software/author-authored source descriptions only; no third-party numerical data grant or physical activation.")
    files = []
    members = [("definitions-json", document), ("definitions-csv", table), ("software-licence", licence)]
    if scenarios:
        # Preserve the original definitions CSV schema. Native clients can join
        # this two-column table by model_id to obtain the same selectors Python
        # exposes on db.model(id). This table contains pointers and model choices.
        bindings = directory / "model-source-bindings.csv"
        with bindings.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(("model_id", "bindings_json"))
            for scenario in scenarios:
                writer.writerow((scenario["id"], json.dumps(scenario["source_bindings"], sort_keys=True)))
        members.append(("source-bindings-csv", bindings))
    for name, path in members:
        entry = dict(id="amps-models:" + name, logical_name=path.name, path=str(path),
            sha256=hashlib.sha256(path.read_bytes()).hexdigest(), role="authored-software-definition",
            source_uri="https://github.com/SWMFsoftware/AMPS", storage_permission="allowed",
            evaluation_permission="allowed", rights=policy, capabilities=["source-access"], physical_activation=False)
        if path == table:
            entry["table"] = dict(id=TABLE_ID, format="csv", schema_id="amps-reduced-model-definitions-v1",
                columns=list(COLUMNS), row_count=len(values),
                description="Explicit-input software definitions and source descriptors; not a physical reaction profile.")
        elif name == "source-bindings-csv":
            entry["table"] = dict(id=BINDINGS_TABLE_ID, format="csv", schema_id=BINDINGS_TABLE_ID,
                columns=["model_id", "bindings_json"], row_count=len(scenarios),
                description="Source selectors, physical scope and birth-model choices for explicitly named reference cases.")
        files.append(entry)
    dataset = dict(id="amps-reusable-models-v1", title="AMPS application-independent explicit-input models",
        provider="AMPS software", version="v1", description="Authored mathematical/reduced-model definitions and limitations; no selected species parameters.",
        availability="present", expected_root=str(directory), rights=policy, issues=[],
        capabilities=["source-access"], files=files)
    return dataset, values
