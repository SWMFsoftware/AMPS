# How EXOSPHID creates He⁺, Ne⁺ and Ar⁺ in a lunar exosphere

The optional dataset `exosphid-noble-gas-photoionization` preserves the
parameters of the EXOSPHID lunar photoionization model at commit
`62e670f9d0344ced1c1414a86b40fbba43b373ef`. Each module describes a
single photon process giving one singly charged atom and one electron. The
only reaction has probability 1 within its wavelength bin. The species
identifiers carry charge and omit the ion's electronic level. This is a
reproducible source for a bulk ionization comparison. [Author’s thesis](https://mediatum.ub.tum.de/doc/1784196/1784196.pdf)

| Parent | Quiet Sun (s⁻¹) | Active Sun (s⁻¹) | Model threshold (eV) | Upper wavelength cutoff (Å) | Pinned author source |
| --- | ---: | ---: | ---: | ---: | --- |
| He | 5.57 × 10⁻⁸ | 16.80 × 10⁻⁸ | 24.59 | 504.27 | [He.jl](https://github.com/charligolea/EXOSPHID/blob/62e670f9d0344ced1c1414a86b40fbba43b373ef/src/database/species/He.jl#L10-L58) |
| Ne | 1.72 × 10⁻⁷ | 5.74 × 10⁻⁷ | 21.57 | 574.94 | [Ne.jl](https://github.com/charligolea/EXOSPHID/blob/62e670f9d0344ced1c1414a86b40fbba43b373ef/src/database/species/Ne.jl#L10-L59) |
| Ar | 3.11 × 10⁻⁷ | 7.15 × 10⁻⁷ | 17.56 | 706.00 | [Ar.jl](https://github.com/charligolea/EXOSPHID/blob/62e670f9d0344ced1c1414a86b40fbba43b373ef/src/database/species/Ar.jl#L10-L59) |

The original modules call the rates `qr` and `ar` and write the numbers
as Julia Float32 literals. The recorded strings preserve their significant
digits. The source specifies wavelength in ångströms and energy in eV. It
gives no outcome probabilities among electronic states and no uncertainty
for any tabulated rate or threshold. The single-channel weight says nothing
about how often an excited ion forms. The shared database does not assign
an electronic level, a spectrum of outgoing electrons or activation on the
basis of these modules.

The Ar threshold deserves special care. NIST ASD gives **15.7596119 eV**
for the first ionization of ground-state neutral Ar. The EXOSPHID constant
is **17.56 eV**, higher by about **1.80 eV**. The module identifies its
product only as Ar(+). There is no source assignment for interpreting this
difference as residual excitation. Both numbers keep their own source and
meaning. An application requiring an Ar ground-state energy budget can
use NIST and evaluate the EXOSPHID threshold as a separate model case.
[NIST ASD first ionization energy](https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=argon&units=1)

## Acquisition, replay and access

Run from the AMPS root. Use fresh output directories. Online acquisition
requires curl. Preparation reads the original bytes as text without running
Julia. Each file is bound to its immutable upstream revision and SHA-256.
Only bibliographic and model metadata enter the shared catalogue. Replay
from the original receipts gives a byte-identical source packet.

```bash
work=$(mktemp -d)
recipe=utility/chemical_database/native_m3_exosphid_atoms.py
python3 -B "$recipe" --fetch --acquisition "$work/originals" --output-dir "$work/packet"
python3 -B "$recipe" --acquisition "$work/originals" --output-dir "$work/replay"
cmp "$work/packet/source_packet.json" "$work/replay/source_packet.json"
```

Pass `$work/packet/source_packet.json` to the database builder as an explicit
source packet. Query a release containing it using the standard API.

```python
for entry in database.describe("exosphid-noble-gas-photoionization")["files"]:
    facts = entry["quantity_metadata"]
    print(facts["atom"], facts["quiet_rate_s_inv"], facts["threshold_eV"])
```

Set `AMPS_CHECK_EXOSPHID_ATOMS_RELEASE=1` when running the shared database
suite against the selected lock. The source checks compare the three sets
of native constants, identify the Ar threshold discrepancy and verify
metadata-only access and the absence of assigned state probabilities.

| Local original | SHA-256 |
| --- | --- |
| He.jl | `33266e99326a251edc1cd4f54d8be2972e5f3e5e01f123f557ae7b871cbe35f3` |
| Ne.jl | `f9780f69abf03784b2d122fff1d48817220bc7f78202cf269a78279c900a4b3b` |
| Ar.jl | `f4cd3611cc0eb67ea2f4e25fd3bcb9adcc76007bbc86853d2a5e4fd4c4aca1de` |
