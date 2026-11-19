"""Independent primary-publication CEX reference, never AMPS output.

Numbers are a versioned transcription of the PRINTED numerical tables/fit,
checked against the original PDF pages (ORNL A-32/PDF48; Morgan p984/PDF15).
No plotted curve, OCR-generated substitute, fit to AMPS or invented point is
used. The writer reproduces the reference artifact from this committed script
after checking both original publication hashes. Native C++ uses Clenshaw;
this oracle independently constructs the Chebyshev polynomials in Decimal.
The printed He table is retained separately from its approximate fitted law.
"""
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path


PUBLICATIONS = {
    "he": {
        "filename": "ornl_redbook_v1_osti.pdf",
        "sha256": "3bb3fdb4be8f0f806d378083be8798cf3ce347714c3d06d58f90ee483b10c030",
        "url": "https://www.osti.gov/servlets/purl/6570226",
        "locator": "ORNL-6086/V1 A-32, PDF48; Appendix I formula",
    },
    "na": {
        "filename": "nist_morgan_1985_retry2.pdf",
        "sha256": "7c412cb93f827b3412b6da000e674c63237cacb0e201733a6f5a4f431359952e",
        "url": "https://srd.nist.gov/jpcrdreprint/1.555752.pdf?download=1",
        "locator": "Morgan et al. 1985, printed984/PDF15, Table2 Na column",
    },
}

HE_COEFFICIENTS = (
    "-92.7819", "-4.69073", "-8.99666", "-0.490484", "0.935276",
    "-0.0504853", "-0.0410750", "-0.0730136", "0.0189320",
)
HE_ENERGY_LIMITS_EV_U = ("99", "11000000")
HE_PRINTED_TABLE_CM2 = (
    ("99", "4.28e-22"), ("200", "2.46e-21"), ("400", "1.77e-20"),
    ("700", "1.12e-19"), ("1000", "4.73e-19"), ("2000", "4.25e-18"),
    ("4000", "1.98e-17"), ("7000", "5.18e-17"), ("10000", "8.52e-17"),
    ("20000", "1.75e-16"), ("24000", "1.83e-16"), ("40000", "1.41e-16"),
    ("70000", "6.98e-17"), ("100000", "3.38e-17"), ("200000", "3.95e-18"),
    ("400000", "3.22e-19"), ("700000", "3.07e-20"), ("1000000", "5.42e-21"),
    ("2000000", "1.42e-22"), ("4000000", "3.29e-24"), ("7000000", "1.62e-25"),
    ("10000000", "2.09e-26"), ("11000000", "1.02e-26"),
)
NA_PRINTED_TABLE_KEV_U = (
    ("1", "43"), ("2", "68"), ("3", "75"), ("5", "66"),
    ("7", "53"), ("10", "24"), ("20", "3.5"), ("50", "0.71"),
    ("70", "0.47"), ("100", "0.32"), ("200", "0.19"),
)  # Cross-section unit: 10^-16 cm2; no populated Na knots below1keV/u.


def helium_fit_m2(energy_ev_u):
    """Original A0/2 Chebyshev series via direct T_n, 70-digit arithmetic."""
    with localcontext() as context:
        context.prec = 70
        energy = Decimal(str(energy_ev_u))
        low, high = map(Decimal, HE_ENERGY_LIMITS_EV_U)
        if not energy.is_finite() or not low <= energy <= high:
            raise ValueError("outside original He source-energy domain")
        x = (2 * energy.ln() - low.ln() - high.ln()) / (high.ln() - low.ln())
        coefficients = tuple(map(Decimal, HE_COEFFICIENTS))
        previous, current = Decimal(1), x
        logarithm = coefficients[0] / 2 + coefficients[1] * x
        for coefficient in coefficients[2:]:
            previous, current = current, 2 * x * current - previous
            logarithm += coefficient * current
        return logarithm.exp() * Decimal("1e-4")


def sodium_loglog_m2(energy_ev_u):
    """Independent log-linear interpolation of the exact positive Na knots.

    This is the selected named interpolation MODEL, not an experimental point.
    Source energy changes keV/u->eV/u once; source sigma changes10^-16cm2->m2.
    """
    with localcontext() as context:
        context.prec = 70
        energy = Decimal(str(energy_ev_u))
        knots = tuple((Decimal(e) * 1000, Decimal(s) * Decimal("1e-20"))
                      for e, s in NA_PRINTED_TABLE_KEV_U)
        if not energy.is_finite() or not knots[0][0] <= energy <= knots[-1][0]:
            raise ValueError("outside populated Na source knots")
        for index, (knot, value) in enumerate(knots):
            if energy == knot:
                return value
            if energy < knot:
                low_energy, low_value = knots[index - 1]
                fraction = (energy.ln() - low_energy.ln()) / (knot.ln() - low_energy.ln())
                return ((1 - fraction) * low_value.ln() + fraction * value.ln()).exp()
        raise AssertionError("source interval not found")


def write_reference(output_directory, research_directory):
    """Retain source evidence and raw numeric transcription without adjustment."""
    identities = {}
    for name, record in PUBLICATIONS.items():
        path = Path(research_directory) / record["filename"]
        raw = path.read_bytes()
        actual = hashlib.sha256(raw).hexdigest()
        if actual != record["sha256"] or not raw.startswith(b"%PDF-"):
            raise ValueError("original primary publication identity failed: " + str(path))
        identities[name] = dict(record, path=str(path.resolve()), verified_sha256=actual)
    result = {
        "schema": "srcMoon-primary-total-capture-transcription-v1",
        "scope": "original source numbers; independent fit/interpolation reference; not observational validation",
        "source_publications": identities,
        "helium_coefficients": HE_COEFFICIENTS,
        "helium_energy_limits_eV_u": HE_ENERGY_LIMITS_EV_U,
        "helium_printed_table_energy_eV_u_cross_section_cm2": HE_PRINTED_TABLE_CM2,
        "sodium_table_energy_keV_u_cross_section_1e-16_cm2": NA_PRINTED_TABLE_KEV_U,
        "helium_source_accuracy_note": "40% below5keV/u;20% above; fit-table RMS10%, maximum18.7%; not Gaussian1sigma",
        "writer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    path = Path(output_directory) / "primary_reference.json"
    with path.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    return result


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--research-directory", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=False)
    write_reference(args.output_directory, args.research_directory)
    print(args.output_directory / "primary_reference.json")
