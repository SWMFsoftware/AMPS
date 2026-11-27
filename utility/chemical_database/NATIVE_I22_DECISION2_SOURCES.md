# I22 decision-2 sources (2026-11-26, r15)

The I22 decision-2 review (`srcMoon/test/M3_I22_DECISION2_PROPOSALS_20261126.md`)
used 19 papers on water and OH energetics, water-ion fragmentation and water
adsorption. This packet records all of them in the shared catalogue under
the [source policy](RIGHTS.md). The user approved this route on 2026-11-26.

`native_i22_decision2.py` builds the packet from an `acquire_sources.py`
directory and the local research directory. It downloads nothing and
selects no value for a model.

| Records | Policy | Payload |
|---|---|---|
| Crossref record for each of the 19 papers (`i22-decision2:<key>:crossref`) | metadata-only | held |
| Repository copies at OSTI (Harper 2019), HAL (Sage 2010), GUP (Hult Roos 2018), and the DORA record page | metadata-only: no licence statement found | held |
| Bodi et al. 2014, Chem. Sci. 5, 3057 (`bodi2014:article`) | `rsc-bodi-2014-cc-by-nc-3.0`: licence printed on every page | stored |
| `bodi2014:oh-plus-appearance-energy`: the four printed OH⁺ onsets (18.1177, 18.1183, 18.119 and the weighted 18.1182 eV, with printed uncertainties) | same as its parent | stored |

Each Crossref record carries:
- the DOI, the citation, and what the review used the paper for;
- any licence URLs that Crossref lists;
- the local copy the review used, if there was one, by path, SHA-256 and
  origin.

Thirteen publisher copies were supplied by the user from institutional
access. They have no retrieval record, so they stay outside the catalogue as
held private copies. None is fabricated into an acquisition.

The two Schörghofer PSJ articles print CC BY 4.0. Their iopscience PDFs
returned HTTP 403 to `acquire_sources.py`, which is retained in the
acquisition manifest. They are therefore metadata-only until a genuine
retrieval succeeds; the publisher's bot check was not bypassed.

The text layer of the Bodi PDF encodes the printed ± as U+0001. The
projection keeps the printed numeric tokens, and the sign was checked on the
rendered page. The table carries only the H₂O → OH⁺ + H + e threshold. It
gives no kinetic-energy release, angular distribution or H⁺/O⁺ onset.

The PHIDRATES solar flux is not included: the service was unreachable on
2026-11-26.

Evidence: `test_output/chemical_database/i22-decision2-sources-GzoR4a/`.
