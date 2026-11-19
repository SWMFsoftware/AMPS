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
)


def definitions():
    """Return a defensive copy; no registry field is a physical default."""
    return copy.deepcopy(_DEFINITIONS)


def materialize(root, release_id, amps_root):
    """Create authored software metadata as ordinary permission-cleared sources.

    This is NOT a download and carries no fabricated HTTP/acquisition record.
    Exact generated JSON/CSV and the repository licence become release members;
    all third-party source-data admission still uses source_packets.py.
    """
    directory = Path(root) / "releases" / release_id / "provenance" / "authored-models"
    directory.mkdir(parents=True, exist_ok=False)
    values = definitions()
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
    for name, path in (("definitions-json", document), ("definitions-csv", table), ("software-licence", licence)):
        entry = dict(id="amps-models:" + name, logical_name=path.name, path=str(path),
            sha256=hashlib.sha256(path.read_bytes()).hexdigest(), role="authored-software-definition",
            source_uri="https://github.com/SWMFsoftware/AMPS", storage_permission="allowed",
            evaluation_permission="allowed", rights=policy, capabilities=["source-access"], physical_activation=False)
        if path == table:
            entry["table"] = dict(id=TABLE_ID, format="csv", schema_id="amps-reduced-model-definitions-v1",
                columns=list(COLUMNS), row_count=len(values),
                description="Explicit-input software definitions and source descriptors; not a physical reaction profile.")
        files.append(entry)
    dataset = dict(id="amps-reusable-models-v1", title="AMPS application-independent explicit-input models",
        provider="AMPS software", version="v1", description="Authored mathematical/reduced-model definitions and limitations; no selected species parameters.",
        availability="present", expected_root=str(directory), rights=policy, issues=[],
        capabilities=["source-access"], files=files)
    return dataset, values
