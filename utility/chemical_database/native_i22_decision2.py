#!/usr/bin/env python3
"""Register the I22 decision-2 source review (2026-11-26) in the shared catalogue.

Inputs are an `acquire_sources.py` directory and the local research directory
that holds copies supplied by the user or saved during the review. Output is
one offline `amps-chemical-source-packet-v1` packet:

* a Crossref record for every paper used in the review, metadata-only;
* repository copies at OSTI, HAL and GUP, held metadata-only because no
  licence statement was found;
* Bodi et al. 2014 from the PSI DORA repository under the licence printed on
  every page (CC BY-NC 3.0 Unported). Its printed OH+ appearance energies are
  projected losslessly into one small table.

Publisher PDFs supplied by the user, and agent copies, are never admitted as
acquired sources: there is no retrieval record for them, and none is created.
Each one is recorded only by local path, SHA-256 and origin in the metadata of
its Crossref record. Failed retrievals stay in the acquisition manifest and are
not packet sources: the Schörghofer PDFs returned 403 from iopscience.

The module downloads nothing, converts no unit, and selects no value for a
model. See SOURCE_PACKET_CONTRACT.md and RIGHTS.md.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

BODI_POLICY = "rsc-bodi-2014-cc-by-nc-3.0"
BODI_URL = "https://www.dora.lib4ri.ch/psi/islandora/object/psi:9276/datastream/PDF/download"
BODI_SHA = "8882498b8d343b4603db90162a56eee4f0e424c94b2bf720e8d9461d8c1eaad2"
BODI_LICENCE = "This article is licensed under a Creative Commons Attribution-NonCommercial 3.0 Unported Licence"
CROSSREF = "https://api.crossref.org/works/"

# key: (DOI, citation, what the review used it for). Facts only; no selection.
PAPERS = {
    "tan1978": ("10.1016/0301-0104(78)85080-0",
                "K. H. Tan, C. E. Brion, Ph. E. van der Leeuw, M. J. van der Wiel (1978), Chem. Phys. 29, 299-309",
                "Tables 2-4: H2O photoabsorption and fragment-ion oscillator strengths, partial oscillator strengths of the four valence orbitals 16.5-60.5 eV, fragment appearance potentials and ion kinetic energies."),
    "kronebusch1976": ("10.1016/0020-7381(76)80088-5",
                       "P. L. Kronebusch, J. Berkowitz (1976), Int. J. Mass Spectrom. Ion Phys. 22, 283-306",
                       "D04 branch source for H2O04-06. Table 7: H2O fragmentation ratios at nine line-source energies 16.76-40.80 eV with energetics."),
    "eland1975": ("10.1016/0301-0104(75)80036-x", "J. H. D. Eland (1975), Chem. Phys. 11, 41-47",
                  "D2O+ breakdown and kinetic-energy release; D2O only."),
    "powis1991": ("10.1039/FT9918700921", "I. Powis, D. J. Reynolds (1991), J. Chem. Soc. Faraday Trans. 87, 921-926",
                  "D2O+ B 2B2 dissociation dynamics by He I PEPICO; D2O only."),
    "banna1986": ("10.1063/1.450008", "M. S. Banna, B. H. McQuaide, R. Malutzki, V. Schmidt (1986), J. Chem. Phys. 84, 4739-4744",
                  "Tables I-IV: asymmetry parameters and partial photoionization cross sections of the water orbitals, 30-140 eV."),
    "truesdale1982": ("10.1063/1.443058",
                      "C. M. Truesdale, S. Southworth, P. H. Kobrin, D. W. Lindle, G. Thornton, D. A. Shirley (1982), J. Chem. Phys. 76, 860-865",
                      "Tables I-II: photoelectron asymmetry parameters of the H2O ionic states."),
    "liu2009": ("10.1063/1.3120443",
                "J. Liu, E. J. Salumbides, U. Hollenstein, J. C. J. Koelemeij, K. S. E. Eikema, W. Ubachs, F. Merkt (2009), J. Chem. Phys. 130, 174306",
                "Ionization and dissociation energies of para-H2: D0 = 36118.06962(37) cm-1 as printed in the abstract."),
    "norwood1991": ("10.1063/1.461334", "K. Norwood, A. Ali, C. Y. Ng (1991), J. Chem. Phys. 95, 8029-8037",
                    "PEPICO breakdown diagrams of energy-selected H2O+ and D2O+ (B 2B2) at 625-700 A; Table I relative ion intensities."),
    "boyarkin2013": ("10.1016/j.cplett.2013.03.007",
                     "O. V. Boyarkin, M. A. Koshelev, O. Aseev, P. Maksyutenko, T. R. Rizzo, N. F. Zobov, L. Lodi, J. Tennyson, O. L. Polyansky (2013), Chem. Phys. Lett. 568-569, 14-20",
                     "D0(H-OH) = 41145.92 +/- 0.12 cm-1 (95 %), continuum onset in total rovibrational energy to OH[X2Pi3/2, e, J=3/2, v=0] + H."),
    "schorghofer2023": ("10.3847/PSJ/acf19b", "N. Schorghofer (2023), Planet. Sci. J. 4, 164",
                        "Water adsorption kinetics on lunar grains; constant prefactor 1e16 s-1. CC BY 4.0 printed, but the PDF returned HTTP 403 to acquisition."),
    "schorghofer2025": ("10.3847/PSJ/ade5b2", "N. Schorghofer (2025), Planet. Sci. J. 6, 164",
                        "Eq. 4 efficient-diffusion desorption rate and Eq. 9 prefactor 2.2e17/sqrt(T) s-1. CC BY 4.0 printed, but the PDF returned HTTP 403 to acquisition."),
    "joens2001": ("10.1021/jp011833u", "J. A. Joens (2001), J. Phys. Chem. A 105, 11041-11044",
                  "D0[OH(X2Pi3/2)] = 35584 +/- 10 cm-1 from thermochemical cycles."),
    "sun2020": ("10.1063/1674-0068/cjcp2003026",
                "G. Sun, X.-F. Zheng, Y. Qin, Y. Song, J. Zhang, J. M. Amero, G. J. Vazquez (2020), Chin. J. Chem. Phys. 33, 129-134",
                "D0(O-H) = 35580 +/- 15 cm-1 measured; cites ATcT 35580 +/- 3 cm-1 (its ref. 46)."),
    "ruscic2015": ("10.1021/acs.jpca.5b01346", "B. Ruscic (2015), J. Phys. Chem. A 119, 7810-7837",
                   "ATcT reference cited by Sun et al. 2020 for D0(OH); no local copy."),
    "bodi2014": ("10.1039/C4SC00791C", "A. Bodi et al. (2014), Chem. Sci. 5, 3057-3063",
                 "OH+ appearance energy from H2O by iPEPICO; the article is admitted separately under its printed licence."),
    "hultroos2018": ("10.1063/1.5055593", "A. Hult Roos, J. H. D. Eland, J. Andersson, R. J. Squibb, R. Feifel (2018), J. Chem. Phys. 149, 204307",
                     "D2O ion dissociation after valence and inner-valence ionization; C-state O+ and D+ from full atomization."),
    "harper2019": ("10.1063/1.5091966",
                   "O. J. Harper, M. Hassenfratz, J.-C. Loison, G. A. Garcia, N. de Oliveira, H. R. Hrodmarsson, S. T. Pratt, S. Boye-Peronne (2019), J. Chem. Phys. 150, 141103",
                   "Absolute OH photoionization cross section, 12.6-15 eV."),
    "sage2010": ("10.1080/00268971003596177", "A. G. Sage, T. A. A. Oliver, R. N. Dixon, M. N. R. Ashfold (2010), Mol. Phys. 108, 945-955",
                 "H2O+ photodissociation by velocity-map imaging; D0(H-OH+) and D0(H+-OH)."),
    "suarez2018": ("10.1039/C8CP03725F", "J. Suarez, L. Mendez, I. Rabadan (2018), Phys. Chem. Chem. Phys. 20, 28511-28522",
                   "Theory of nonadiabatic H2O+ fragmentation; no local copy (retrieval of a related paper returned 403)."),
}

# Local copies, relative to the research directory: (path, SHA-256, origin).
LOCAL = {
    "tan1978": ("user-supplied-20261126/1-s2.0-0301010478850800-main.pdf",
                "3e6d839dc3c8ae9fc3e2d6eb259216d28d70b1168a145538f4163a39db13efbd", "user"),
    "kronebusch1976": ("user-supplied-20261126/1-s2.0-0020738176800885-main.pdf",
                       "31f79c84c79aae8f92acb4d3298930305ef790f7b853e65d21128ca0a078a3e0", "user"),
    "eland1975": ("user-supplied-20261126/1-s2.0-030101047580036X-main.pdf",
                  "6d37fde3e2a37d214211acb4cefa7a4eecd72e5c37b64d8098e7f2dddf62e2a2", "user"),
    "powis1991": ("user-supplied-20261126/ft9918700921.pdf",
                  "35b93456e39c18e82e0f3e130f31829e18728fe012cc3577e9a6b6f243d21e8d", "user"),
    "banna1986": ("user-supplied-20261126/4739_1_online.pdf",
                  "1f4c823d1bb4be1a809bbab07b1eeefbc3bed8aebeb85f8de7371603df1343aa", "user"),
    "truesdale1982": ("user-supplied-20261126/860_1_online.pdf",
                      "c5a88dde782d4c9db51f86adb08f3f7371c3ff8bca721c527d060f561cff71f1", "user"),
    "liu2009": ("user-supplied-20261126/174306_1_online.pdf",
                "d79e96e2584203a0f92e5845f0acab2d863969f1b9631dee16bec2e0d83baf1e", "user"),
    "norwood1991": ("user-supplied-20261126/8029_1_online.pdf",
                    "05d0d978357c2f30992465f2d5bcfc3812b703650844d9a89b16f60f7795b7d8", "user"),
    "boyarkin2013": ("user-supplied-20261126/1-s2.0-S0009261413003205-main.pdf",
                     "fba516bc1ebd32d4ecfbb8dfd6b259cbda84fa22bfe816a6c2cf3943fbc4620d", "user"),
    "schorghofer2023": ("user-supplied-20261126/Schörghofer_2023_Planet._Sci._J._4_164.pdf",
                        "6b79c114d83049e3db52ee2e01fd1a6380c1cdef58dab2a7a269c963f2ca6f72", "user"),
    "schorghofer2025": ("user-supplied-20261126/Schörghofer_2025_Planet._Sci._J._6_164.pdf",
                        "307f737200733a7bc9126eac834fb0313093370e6c06791ac5a98e457e5f22a5", "user"),
    "joens2001": ("user-supplied-20261126/jp011833u.pdf",
                  "8a7516f4e558c603f5e265da48edab0a70b83ae42c23c93d7070778e7367347f", "user"),
    "sun2020": ("user-supplied-20261126/129_1_1674-0068_cjcp2003026.pdf",
                "ca4ba76ee1dca821dfb248c04a6b75819a13cf07a8a14a0e045311e081466f37", "user"),
    "bodi2014": ("downloads/dl_bodi/bodi2014.pdf", BODI_SHA, "review agent"),
    "hultroos2018": ("downloads/dl2/gu207611.pdf",
                     "11e711e10dae6fe12b26b2414640f88ec3e0fc572198263aa2aed37dca4f1626", "review agent"),
    "harper2019": ("downloads/dl2/harper_osti.pdf",
                   "66a066a60ed07ff320afa72f11ee28c45641299ccbcaf854c699d7b18fd5554c", "review agent"),
    "sage2010": ("downloads/dl3/sage2010.pdf",
                 "85afbb10a7234b21978340fb2a0c445548ac60ad12a6ee591f735587354495e9", "review agent"),
}

# Open repository copies acquired in this packet: (acquisition path, key).
# No licence statement was found in them, so they are held.
HELD_COPIES = {
    "pdf/harper2019-osti.pdf": "harper2019",
    "pdf/sage2010-hal.pdf": "sage2010",
    "pdf/hultroos2018-gup.pdf": "hultroos2018",
    "html/bodi2014-dora-record.html": "bodi2014",
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def bodi_onsets(path):
    """The printed OH+ appearance energies, from the PDF text layer (no OCR).

    The text layer encodes the printed plus-minus sign as U+0001; the sign was
    checked against the rendered page 4. Numeric tokens are kept as printed.
    """
    text = subprocess.run(["pdftotext", "-raw", "-f", "4", "-l", "4", str(path), "-"],
                          capture_output=True, encoding="utf-8", check=True, timeout=30).stdout
    pm = "[\u0001±]"
    sets = re.findall(r"determined to be (\d+\.\d+) " + pm + r" (\d+\.\d+) eV, (\d+\.\d+) " + pm +
                      r" (\d+\.\d+) eV, and\s+(\d+\.\d+) " + pm + r" (\d+\.\d+) eV, for the three data sets in Fig\. 2a", text)
    mean = re.findall(r"E0 value of (\d+\.\d+) " + pm + r" (\d+\.\d+) eV", text)
    if len(sets) != 1 or len(mean) != 1:
        raise ValueError("Bodi 2014 page 4 onset text not found exactly once")
    a = sets[0]
    rows = [("Fig. 2a data set", a[0], a[1]), ("Fig. 2b data set", a[2], a[3]),
            ("Fig. 2c data set", a[4], a[5]),
            ("reciprocal error-square weighted average", mean[0][0], mean[0][1])]
    if [r[1] for r in rows] != ["18.1177", "18.1183", "18.119", "18.1182"]:
        raise ValueError("unexpected printed Bodi onsets")
    return rows


def prepare(acquisition, research, output):
    acquisition, research, output = Path(acquisition), Path(research), Path(output)
    manifest = json.loads((acquisition / "acquisition_manifest.json").read_text())
    receipts = {r["path"]: r for r in manifest["records"]}
    output.mkdir(parents=True, exist_ok=False)
    recipe = output / Path(__file__).name
    shutil.copyfile(__file__, recipe)
    sources = []

    def acquired(member):
        r = receipts[member]
        path = acquisition / member
        if (r["status"] != "PASS" or r["http_status"] != 200 or not r.get("complete", True)
                or r["sha256"] != sha256(path) or r["byte_count"] != path.stat().st_size):
            raise ValueError("missing, failed or changed acquisition: " + member)
        target = output / member.replace("/", "-")
        shutil.copyfile(path, target)
        return r, target

    def local_copy(key):
        if key not in LOCAL:
            return None
        relative, digest, origin = LOCAL[key]
        path = research / relative
        if sha256(path) != digest:
            raise ValueError("local copy changed: " + relative)
        return dict(path=str(path), sha256=digest, supplied_by=origin, supplied_utc_date="2026-11-26",
                    storage="held: not admitted, no retrieval record; private research copy outside the catalogue")

    for key, (doi, citation, use) in PAPERS.items():
        receipt, target = acquired("crossref/" + key + ".json")
        message = json.loads(target.read_text())["message"]
        if message["DOI"].lower() != doi.lower() or receipt["url"] != CROSSREF + doi:
            raise ValueError("Crossref record does not match DOI: " + key)
        licences = sorted({item["URL"] for item in message.get("license", [])})
        sources.append(dict(
            id="i22-decision2:" + key + ":crossref", logical_name=target.name, path=target.name,
            sha256=receipt["sha256"], source_uri=receipt["url"], final_url=receipt["final_url"],
            http_status=200, retrieved_utc=receipt["retrieved_utc"], policy="metadata-only",
            role="source-metadata", rights_evidence=[], physical_activation=False,
            quantity_metadata=dict(doi=doi, citation=citation, review_use=use,
                                   crossref_licence_urls=licences, local_copy=local_copy(key),
                                   limitation="Bibliographic record only. No value in the paper is selected for a model by this record.")))

    for member, key in HELD_COPIES.items():
        receipt, target = acquired(member)
        sources.append(dict(
            id="i22-decision2:" + key + ":" + member.split("/")[1].rsplit(".", 1)[0], logical_name=target.name,
            path=target.name, sha256=receipt["sha256"], source_uri=receipt["url"], final_url=receipt["final_url"],
            http_status=200, retrieved_utc=receipt["retrieved_utc"], policy="metadata-only", role="source-metadata",
            rights_evidence=[], physical_activation=False,
            quantity_metadata=dict(doi=PAPERS[key][0], citation=PAPERS[key][1],
                                   limitation="Repository copy without a located licence statement; held metadata-only.")))

    receipt, article = acquired("pdf/bodi2014-dora.pdf")
    first_page = subprocess.run(["pdftotext", "-l", "1", str(article), "-"], capture_output=True,
                                encoding="utf-8", check=True, timeout=30).stdout
    if (receipt["url"] != BODI_URL or receipt["sha256"] != BODI_SHA
            or BODI_LICENCE not in " ".join(first_page.split())):
        raise ValueError("Bodi 2014 article or its printed licence changed")
    evidence = [dict(path=article.name, sha256=BODI_SHA, source_uri=BODI_URL)]
    metadata = dict(doi=PAPERS["bodi2014"][0], citation=PAPERS["bodi2014"][1],
                    licence="CC BY-NC 3.0 Unported", licence_url="https://creativecommons.org/licenses/by-nc/3.0/",
                    licence_evidence="Printed on the article pages: '" + BODI_LICENCE + "'.",
                    method="imaging photoelectron photoion coincidence (iPEPICO), SLS X04DB",
                    locator="Printed page 3060 (PDF page 4), Results, OH+ appearance energy",
                    local_copy=local_copy("bodi2014"),
                    limitation="Threshold of H2O -> OH+ + H + e only. No kinetic-energy release, angular distribution or H+/O+ onset.")
    parent = dict(id="bodi2014:article", logical_name=article.name, path=article.name, sha256=BODI_SHA,
                  source_uri=BODI_URL, final_url=receipt["final_url"], http_status=200,
                  retrieved_utc=receipt["retrieved_utc"], policy=BODI_POLICY, role="source-native",
                  rights_evidence=evidence, physical_activation=False, quantity_metadata=metadata)
    sources.append(parent)
    rows = bodi_onsets(article)
    table = output / "bodi2014-oh-plus-appearance-energy.csv"
    with table.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("measurement", "E0_eV_native", "uncertainty_eV_native"))
        writer.writerows(rows)
    sources.append(dict(
        id="bodi2014:oh-plus-appearance-energy:csv", logical_name=table.name, path=table.name,
        sha256=sha256(table), source_uri=BODI_URL, role="derived-table", policy=BODI_POLICY,
        derived_from=["bodi2014:article"], rights_evidence=evidence,
        recipe=dict(path=recipe.name, sha256=sha256(recipe)),
        quantity_metadata=dict(metadata, reaction="H2O + photon -> OH+ + H + e",
                               uncertainty_meaning="as printed beside each value; the paper describes no confidence level",
                               text_layer_note="The text layer encodes the printed plus-minus sign as U+0001; checked on the rendered page."),
        table=dict(id="bodi2014:oh-plus-appearance-energy", name="bodi2014-oh-plus-appearance-energy",
                   format="csv", columns=["measurement", "E0_eV_native", "uncertainty_eV_native"],
                   row_count=len(rows), schema_id="bodi2014-oh-plus-appearance-energy-v1"),
        physical_activation=False))

    failed = sorted(r["path"] + " (HTTP " + str(r.get("http_status")) + ")" for r in manifest["records"] if r["status"] != "PASS")
    packet = dict(
        schema="amps-chemical-source-packet-v1", id="i22-decision2-sources-20261126",
        title="I22 decision-2 sources: water and OH energetics, charged breakup and water adsorption",
        provider="Crossref; PSI DORA; OSTI; HAL; University of Gothenburg GUP",
        description="Bibliographic records for every paper used in the I22 decision-2 review, held repository copies, and the CC BY-NC Bodi 2014 article with its printed OH+ appearance energies.",
        version="20261126",
        limitations=["Source access only. No reaction, threshold or surface model is selected or activated.",
                     "Publisher copies supplied by the user are recorded by SHA-256 only and are not in the catalogue.",
                     "Failed retrievals kept in the acquisition manifest, not admitted: " + "; ".join(failed),
                     "The PHIDRATES solar flux is not included; the service was unreachable on 2026-11-26."],
        sources=sources, physical_activation=False)
    (output / "source_packet.json").write_text(json.dumps(packet, indent=2, ensure_ascii=False) + "\n")
    return dict(sources=len(sources), tables=1, failed_retrievals=failed)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--acquisition", type=Path, required=True)
    parser.add_argument("--research-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.acquisition, args.research_dir, args.output_dir), indent=2))


if __name__ == "__main__":
    main()
