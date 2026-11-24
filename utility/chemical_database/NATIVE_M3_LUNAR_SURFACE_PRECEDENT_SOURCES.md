# How published lunar models treat returning OH

The shared dataset `lunar-oh-surface-model-precedents` records three published modeling choices. Applications can query them to design sensitivity studies. The available sources do not establish a joint measurement of sticking and translational energy transfer for incident OH on lunar grains. Each record declares its event, species and limitations. None of the three sources activates a chemistry channel or changes an AMPS default.

| Source | Event and choices | What a comparison can test |
| --- | --- | --- |
| Jones et al. (2018), model terms k5–k9 | The authors assume half the water photofragments escape. For those that return, k7 assumes that every incident OH radical forms bound M–OH, with temperature-independent sticking probability 1 and no barrier. Radicals land at random sites on the sunlit side in this model. | A reactive surface sink for returning radicals. The bound inventory can produce molecular water by recombination. An intact-OH re-emission kernel requires additional physics. [Primary paper](https://doi.org/10.1029/2018GL080008) |
| Smolka (2023), LPSC abstract 2089, numerical model | Full thermal accommodation and Maxwell–Boltzmann flux emission for neutral H, H₂, OH and H₂O on a spherical Moon with a stationary temperature pattern relative to the subsolar point. A three-dimensional thermal flux has speed density proportional to `v³ exp(−mv²/(2kT))`. | A fully accommodated thermal launch for comparison after surface contact. The author gives no measured OH sticking probability in this paragraph. Gas-phase photofragments need an independent recoil model before first contact. [Author’s conference report](https://www.hou.usra.edu/meetings/lpsc2023/pdf/2089.pdf) [Related journal article](https://doi.org/10.1016/j.icarus.2023.115508) |
| Laferriere et al. (2025), Figshare release v1 | In the smooth-surface calculation, OH and H₂O use the same residence function. The archived source sets `E_d = 0.6 eV`, `A = 10¹³ s⁻¹`, `f = 0.01` and computes `τ = 1/[f A exp(−E_d/(k_B T))]`. A launched particle draws a new speed from discrete weights proportional to `v² exp(−mv²/(2k_B T))`. | A mobile OH hopping scenario. Preserve all three residence parameters and the source speed law when reproducing the calculation. Landing and relaunch supply no independent energy-resolved sticking measurement. [Journal article](https://doi.org/10.1029/2024JE008628) [Versioned author archive](https://doi.org/10.6084/m9.figshare.28548857.v1) |

Energy accommodation describes the extent to which the outgoing translational energy approaches equilibrium with the surface. A probability of entering a chemically bound state, a residence time and a launch-speed law describe different events. The incident OH beam energy and regolith composition, temperature, site availability and coverage would be needed to turn a sensitivity comparison into a surface-scattering qualification. Experiments on solar-wind proton implantation and the resulting subsurface OH constrain a separate process. The evaluated `iupac2010-oh-ice:uptake` table covers gas-phase radical loss on ice under its reported laboratory conditions. [Tucker et al.](https://doi.org/10.1029/2018JE005805) [IUPAC evaluation](https://doi.org/10.5194/acp-10-9059-2010)

## Reproduce and query

Run from the AMPS root with Python 3, Poppler’s `pdftotext` and curl. Choose a fresh directory. Online acquisition checks each SHA-256 before the builder assigns a source identity. Crossref supplies a bibliographic record for Jones. The k5–k9 interpretation is curated from the primary publisher text, and its URL and paragraph locator are retained independently of the Crossref hash. The Smolka extraction verifies the author's assumptions on page 1. The Laferriere extraction inspects the named routines in the versioned Python file as text. Acquisition omits ephemeral query credentials from the final Figshare redirect. The public packet distributes scope and parameter metadata under the existing metadata-only policy.

~~~bash
workdir=$(mktemp -d)
python3 -B utility/chemical_database/native_m3_lunar_surface_precedent.py \
  --fetch --acquisition "$workdir/acquired" --output-dir "$workdir/packet"
python3 -B utility/chemical_database/native_m3_lunar_surface_precedent.py \
  --acquisition "$workdir/acquired" --output-dir "$workdir/replay"
cmp "$workdir/packet/source_packet.json" "$workdir/replay/source_packet.json"
~~~

The acquisition receipt records the HTTP status, retrieval time, source URI, final destination, size and SHA-256. The original private inputs and packet are preserved under `test_output/chemical_database/lunar-oh-precedent-twDfB4/`. The initial packet digest is `efa1e0e4df82ae09fab80542cdd5e72a9c70bf610c6d427691d3b85605098202`. A new network retrieval has a new timestamp. A changed source requires review before a new hash is accepted.

~~~bash
db_root=/data/vtenishe/CHEMICAL_DATABASE
PYTHONPATH="$db_root/python" python3 -B - "$db_root/database-lock.json" <<'PY'
import sys
from amps_chemical_database import Database

db = Database.from_lock(sys.argv[1])
for entry in db.describe("lunar-oh-surface-model-precedents")["files"]:
    print(entry["id"], dict(entry["quantity_metadata"]))
PY
~~~

## Regression contract

When the installed release includes this packet, the shared database suite selects `AMPS_CHECK_LUNAR_SURFACE_PRECEDENT_RELEASE=1`. The three checks independently verify the Jones escape and impact probabilities, both thermal launch laws and the complete archived residence parameter set, and the public metadata/access classification. The oracle comprises Jones k6/k7, the 2089 numerical-model paragraph and Figshare v1's `surftime`, `desorption_R21`, `ballistic` and `maxwell_boltz_dist`. A missing selected record raises an error, and an unexpected published value fails the test. With no selected release, these checks skip. Success verifies documentary provenance and catalogue behavior. M3 physical qualification, including the existing linked-test predecessors, remains separate.
