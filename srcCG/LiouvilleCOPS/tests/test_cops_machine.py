"""Machine-specific "@ <machine>" input lines (tools/cops_machine.py; the same rule as src/Input.cpp)."""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import cops_machine as m  # noqa: E402
from rebuild_cops_boundary import parse_input  # noqa: E402

TEXT = "A = g1\nA = s1 @ x*\nB = g2\nK = k1\nK = k2\nK = kx @ x1, y\nC = mail@host\nD = only @ zz  # c\n"


def values(names):
    return parse_input(m.resolve_text(TEXT, names)[0])


def test_specific_lines_replace_generic_ones_lists_included():
    v = values(["x1"])
    assert v["A"] == ["s1"] and v["K"] == ["kx"] and v["B"] == ["g2"]
    assert "D" not in v                                   # a line for another machine is ignored


def test_no_matching_machine_keeps_the_generic_lines():
    v = values(["other"])
    assert v["A"] == ["g1"] and v["K"] == ["k1", "k2"]


def test_at_without_a_blank_is_part_of_the_value():
    assert values(["x1"])["C"] == ["mail@host"]


def test_patterns_are_case_insensitive_globs():
    assert m.matches(["Voyager*"], ["voyager1"]) and not m.matches(["passion"], ["voyager1"])


def test_resolved_text_has_no_qualifiers():
    text, info = m.resolve_text(TEXT, ["y"])
    assert all("@ " not in line.split("#", 1)[0] for line in text.splitlines())
    assert info["D"] == {"machines": ["zz"], "applied": False}
    assert "only for the machines zz" in m.missing_message("D", info, ["y"])


def test_at_without_a_machine_name_is_an_error():
    with pytest.raises(ValueError):
        m.resolve_text("A = v @ \n", ["x"])


def test_cops_machine_replaces_the_detected_names(monkeypatch):
    monkeypatch.setenv("COPS_MACHINE", "Passion")
    assert m.machine_names() == ["passion"]


@pytest.mark.parametrize("names, kernels", [(["voyager1"], "/home/vtenishev/"), (["gpu001", "voyager"], "/home/vtenishev/"),
                                            (["passion"], "/home/vtenishe/")])
def test_fit_base_in_paths_per_machine(names, kernels):
    v = parse_input(m.resolve_text((ROOT / "fit_base.in").read_text(), names)[0])
    assert len(v["KernelDirectory"]) == 1 and v["KernelDirectory"][0].startswith(kernels)
    assert len(v["ObservationFile"]) == 1
