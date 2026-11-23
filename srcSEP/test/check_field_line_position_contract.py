#!/usr/bin/env python3
"""SCOPE02: permit a Cartesian cache derived solely from the 1-D coordinate.

This is a source-contract audit, not a transport solver or a numerical oracle.
The native SINT08 gate separately compares actual stored Cartesian positions
with the geometry service after physical movement on real MPI decompositions.
The source gate prevents an unrelated Cartesian advance from acquiring the
same permission merely because the application must keep PIC's X view current.
"""

from pathlib import Path
import re
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
ADVANCE_FILES = (
    "sep.h", "mover_state.cpp", "focused_transport_dmumu.cpp",
    "focused_transport_mfp.cpp", "parker_mover.cpp", "transport_common.cpp",
    "production_mover_runtime.cpp",
)
# Mask comments and quoted literals before inspecting executable spellings.
# Braces inside diagnostics must not change the selected function's extent,
# and an explanatory comment must not satisfy an ownership/data-flow check.
NON_CODE = re.compile(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'', re.S)
POSITION_WRITE = re.compile(r"(?:PIC::)?(?:ParticleBuffer|PB)::SetX\b")
SIGNATURE = "StatusCommitAndAttach(constParticleContext&context){"
PROJECTION = "line.GetCartesian(position,context.state.coordinate);"
COMMIT = ("PIC::ParticleBuffer::SetFieldLineCoord(context.state.coordinate,context.data);"
          "PIC::ParticleBuffer::SetX(position,context.data);")


def code(source):
    """Ignore formatting/non-code while keeping C++ identifiers distinct.

    Whitespace-free text is sufficient for exact call/owner statements, but
    identifier-use counts below use the masked, spaced source. Otherwise
    ``double component:position`` would merge two distinct identifiers.
    """
    return NON_CODE.sub(lambda match: " " * len(match.group()), source)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def audit(sources):
    """Check every audited advance unit and the one permitted projection.

    Exactly one SetX reference is permitted across the entire advance surface
    (including address-taking, so an extra function-pointer writer rejects). It
    must live in CommitAndAttach, project the immutable owning field line at
    the very coordinate written to the same particle, and publish that view
    beside the coordinate commit. The local vector has exactly its four
    reviewed uses: declaration, geometry output, finite-value check and commit.
    Any extra access (including an alias passed to a Cartesian pusher or a
    solar-wind displacement added after projection) rejects rather than
    guessing that it leaves the field-line authority intact.
    """
    require(set(sources) == set(ADVANCE_FILES), "advance-file catalog differs")
    masked = {name: code(text) for name, text in sources.items()}
    compact = {name: re.sub(r"\s+", "", text) for name, text in masked.items()}
    for name, text in compact.items():
        for forbidden in ("Relativistic::Boris", "calculatePerpendicularDiffusion",
                          "_PIC_PARTICLE_LIST_ATTACHING_NODE_"):
            require(forbidden not in text, name + ": forbidden 3-D advance " + forbidden)
    writes = [(name, match.start()) for name, text in compact.items()
              for match in POSITION_WRITE.finditer(text)]
    require(len(writes) == 1 and writes[0][0] == "transport_common.cpp",
            "exactly one derived SetX is allowed across all advance files")
    whole = compact["transport_common.cpp"]
    require(whole.count(SIGNATURE) == 1, "missing/ambiguous CommitAndAttach owner")
    opening = whole.index(SIGNATURE) + len(SIGNATURE) - 1
    depth, closing = 0, None
    for index in range(opening, len(whole)):
        depth += (whole[index] == "{") - (whole[index] == "}")
        if depth == 0:
            closing = index
            break
    require(closing is not None, "unterminated CommitAndAttach")
    require(opening < writes[0][1] < closing, "SetX is outside the coordinate commit")
    body = whole[opening:closing + 1]
    owner = "auto&line=PIC::FieldLine::FieldLinesAll[context.state.fieldLineId];"
    require(owner in body,
            "projection does not use the particle's owning field line")
    require("line.GetSegment(context.state.coordinate)!=context.segment" in body,
            "projection no longer validates coordinate/segment ownership")
    require(body.count(PROJECTION) == 1, "projection must use the committed coordinate")
    require(COMMIT in body, "coordinate and Cartesian view must commit to the same particle")
    require(body.index(PROJECTION) < body.index(COMMIT), "projection must precede publication")
    require("doublebegin[3],end[3],position[3];" in body,
            "projection must populate its own local three-component SI position")
    require("for(doublecomponent:position)if(!std::isfinite(component))returnStatus::Error(" in body,
            "nonfinite projected Cartesian positions must reject before publication")
    # Count true identifiers in the selected body, not occurrences of a word
    # inside an error message. Four uses allow no independent Cartesian update,
    # extra pointer escape or replacement position between GetCartesian/SetX.
    spaced = masked["transport_common.cpp"]
    start = re.search(r"\bStatus\s+CommitAndAttach\s*\(\s*const\s+ParticleContext\s*&\s*context\s*\)\s*\{", spaced)
    require(start is not None, "missing typed const particle context")
    depth, end = 1, start.end()
    while end < len(spaced) and depth:
        depth += (spaced[end] == "{") - (spaced[end] == "}")
        end += 1
    require(len(re.findall(r"\bposition\b", spaced[start.end():end])) == 4,
            "projected position has an extra/missing use or independent update")
    # The line reference has only declaration, segment-count/bounds check,
    # segment-owner check and projection uses. This also rejects a hidden
    # second binding or geometry mutation while committing particle state.
    require(len(re.findall(r"\bline\b", spaced[start.end():end])) == 4,
            "field-line owner has an extra/missing binding or mutation")
    require(body.index(owner) < body.index(PROJECTION), "owning line must bind before projection")
    for statement in (owner, PROJECTION, COMMIT):
        location = body.index(statement)
        require(sum((c == "{") - (c == "}") for c in body[:location]) == 1,
                "projection/commit cannot be hidden in a conditional block")
        require(body[location - 1] in ";{}", "projection/commit is not an unconditional statement")


def negative_controls(sources):
    """Mutate only in-memory software fixtures; production stays immutable.

    Each mutation represents a separate forbidden mechanism, not a numerical
    tolerance perturbation. Failure of the source audit to reject any one is
    itself a failing SCOPE02 gate. Real physics accuracy remains owned by the
    unchanged transport/native tests rather than by these string fixtures.
    """
    original = sources["transport_common.cpp"]
    mutations = {
        "wrong-owner": original.replace("FieldLinesAll[context.state.fieldLineId]", "FieldLinesAll[0]"),
        "wrong-projection-coordinate": original.replace("GetCartesian(position,context.state.coordinate)", "GetCartesian(position,0.0)"),
        "extra-wind-displacement": original.replace("line.GetCartesian(position,context.state.coordinate);", "line.GetCartesian(position,context.state.coordinate);position[0]+=400000.0;"),
        "escaped-position-pointer": original.replace("line.GetCartesian(position,context.state.coordinate);", "line.GetCartesian(position,context.state.coordinate);AdvanceCartesian(position);"),
        "wrong-coordinate-commit": original.replace("SetFieldLineCoord(context.state.coordinate, context.data)", "SetFieldLineCoord(0.0, context.data)"),
        "wrong-position-commit": original.replace("SetX(position,context.data)", "SetX(begin,context.data)"),
        "wrong-particle-commit": original.replace("SetX(position,context.data)", "SetX(position,otherParticle)"),
        "conditional-projection": original.replace("line.GetCartesian(position,context.state.coordinate);", "if(false){line.GetCartesian(position,context.state.coordinate);}"),
        "extra-cartesian-write": original + "\nvoid invalid(){PB::SetX(extra,data);}\n",
        "aliased-cartesian-write": original + "\nauto invalid=&PIC::ParticleBuffer::SetX;\n",
        "extra-geometry-mutation": original.replace("line.GetCartesian(position,context.state.coordinate);", "line.Advance();line.GetCartesian(position,context.state.coordinate);"),
        "boris": original + "\nvoid invalid(){Relativistic::Boris();}\n",
        "perpendicular-diffusion": original + "\nvoid invalid(){calculatePerpendicularDiffusion();}\n",
        "amr-particle-list": original + "\nint invalid=_PIC_PARTICLE_LIST_ATTACHING_NODE_;\n",
    }
    for name, changed in mutations.items():
        require(changed != original, "negative control did not mutate its fixture: " + name)
        candidate = dict(sources, **{"transport_common.cpp": changed})
        try:
            audit(candidate)
        except AssertionError:
            continue
        raise AssertionError("source audit accepted forbidden mutation: " + name)
    # A second-file writer catches a whole-file exception accidentally added
    # around transport_common.cpp and proves all audited movers remain covered.
    candidate = dict(sources)
    candidate["parker_mover.cpp"] += "\nvoid invalid(){PIC::ParticleBuffer::SetX(extra,data);}\n"
    try:
        audit(candidate)
    except AssertionError:
        return len(mutations) + 1
    raise AssertionError("source audit accepted another mover's Cartesian write")


def main():
    start = time.monotonic()
    # The shell's legacy mechanism scan and this data-flow audit must cover
    # precisely the same units. A future source added to only one catalog is
    # a gate error, not permission to leave its Cartesian writes unaudited.
    if len(sys.argv) > 1:
        supplied = [Path(path).resolve() for path in sys.argv[1:]]
        require(len(supplied) == len(ADVANCE_FILES) and
                set(supplied) == {(ROOT / name).resolve() for name in ADVANCE_FILES},
                "shell/data-flow advance-file catalogs differ")
    sources = {name: (ROOT / name).read_text(encoding="utf-8") for name in ADVANCE_FILES}
    audit(sources)
    count = negative_controls(sources)
    # These controls belong to the existing SCOPE02 case; do not invent a new
    # public test identity for a source-audit implementation detail.
    print("SCOPE02 projection contract: PASS; one derived Cartesian view; negative_controls=%d duration_s=%.6f" %
          (count, time.monotonic() - start))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AssertionError, OSError) as error:
        print("FAIL SCOPE02: projection contract: " + str(error), file=sys.stderr)
        sys.exit(1)
