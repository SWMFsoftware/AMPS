#!/usr/bin/env python3
"""Convert pre-2026-10 SEP input files to the generic AMPS runtime input format.

The generic format (``#begin KIND[: NAME] ... #end KIND`` blocks, ``!``
comments, ``#include``) is specified in ``src/INPUT_FILE.md``.  This one-off
migration aid rewrites the retired formats while keeping every comment and
every assignment's spelling, so a converted file means exactly what the old
file meant.  Each old INI section is moved into the block that now owns it;
see the ``LAYOUTS`` tables below (they are the authoritative old-to-new
mapping and are reproduced in src/INPUT_FILE.md).

Modes
-----
event      reduced shock-front ``.event`` record (``key=value`` lines, ``#``
           comments) -> one ``#begin physics: shock_front`` block.
cme-event  sep_corona_swcme CME background ``event.conf`` record -> one
           ``#begin physics: cme_background`` block (spelling kept verbatim).
srcsep3d   srcSEP3D INI schema 1-5 deck (``--input``) -> shared SEP schema 6.
srcsep     srcSEP INI schema 1-4 deck (``--input``) -> shared SEP schema 6.
sep6-from-srcsep3d  srcSEP3D generic-format deck of schemas 1-5
           (``#begin application: sep3d``) -> shared SEP schema 6.
sep6-from-srcsep    srcSEP generic-format deck of schemas 1-4
           (``#begin application: sep``) -> shared SEP schema 6.

Usage::

    python3 tools/convert_amps_input.py MODE OLD_FILE [-o NEW_FILE]
    python3 tools/convert_amps_input.py MODE --in-place FILE...

Comments: an old ``#`` comment (whole line or after a value) becomes ``!``.
Comments directly above a section header or key stay attached to it.  The
tool refuses unknown sections instead of guessing an owner.

Requires Python 3.8+, standard library only.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

BlockPath = Tuple[Tuple[str, str], ...]   # ((kind, name), ...) from the root


# ---------------------------------------------------------------------------
# Old-format lexing
# ---------------------------------------------------------------------------

def split_comment(line: str, markers: str = "#") -> Tuple[str, str]:
    """Split ``line`` at the first comment marker; return (code, comment).

    The comment keeps its text without the marker; code is right-stripped."""
    positions = [line.find(m) for m in markers if m in line]
    if not positions:
        return line.rstrip(), ""
    cut = min(positions)
    return line[:cut].rstrip(), line[cut + 1:].rstrip()


def bang(comment: str) -> str:
    """Render an old comment body as a ``!`` comment (blank stays blank)."""
    return "!" + comment if comment else "!"


class Item:
    """One old assignment with the comments written directly above it."""

    def __init__(self, key: str, value: str, trailing: str,
                 before: List[str]):
        self.key, self.value, self.trailing, self.before = key, value, trailing, before


class Section:
    def __init__(self, name: str, before: List[str]):
        self.name, self.before, self.items = name, before, []  # type: ignore[var-annotated]
        self.after: List[str] = []   # comments after the last item


def lex_ini(text: str, markers: str = "#") -> Tuple[List[str], List[Section], List[str]]:
    """Return (file header comments, sections, trailing comments).

    Blank lines are kept as empty strings inside comment groups so the
    converted file keeps the original paragraphing."""
    header: List[str] = []
    sections: List[Section] = []
    pending: List[str] = []
    for raw in text.splitlines():
        code, comment = split_comment(raw, markers)
        code = code.strip()
        if not code:
            pending.append(bang(comment) if (comment or raw.strip()) else "")
            continue
        if code.startswith("[") and code.endswith("]"):
            name = code[1:-1].strip().lower()
            if not sections:
                # Split leading comments: a blank line separates the file
                # header from the comments that describe the first section.
                cut = len(pending)
                while cut > 0 and pending[cut - 1] != "":
                    cut -= 1
                header, pending = pending[:cut], pending[cut:]
            else:
                sections[-1].after = []
            sections.append(Section(name, strip_blank_edges(pending)))
            pending = []
            continue
        if "=" not in code:
            raise SystemExit("cannot convert line without '=': " + raw)
        if not sections:
            raise SystemExit("assignment before any [section]: " + raw)
        key, value = (part.strip() for part in code.split("=", 1))
        sections[-1].items.append(
            Item(key, value, bang(comment) if comment else "",
                 strip_blank_edges(pending)))
        pending = []
    trailing = strip_blank_edges(pending)
    return strip_blank_edges(header), sections, trailing


def strip_blank_edges(lines: List[str]) -> List[str]:
    while lines and lines[0] == "":
        lines = lines[1:]
    while lines and lines[-1] == "":
        lines = lines[:-1]
    return lines


# ---------------------------------------------------------------------------
# Generic-format emission
# ---------------------------------------------------------------------------

class Node:
    """A block of the output tree."""

    def __init__(self, kind: str, name: str = ""):
        self.kind, self.name = kind, name
        self.before: List[str] = []
        self.lines: List[str] = []       # rendered assignment/comment lines
        self.children: List["Node"] = []
        self.after: List[str] = []       # extra raw lines (e.g. #include)

    def label(self) -> str:
        return f"{self.kind}: {self.name}" if self.name else self.kind

    def child(self, kind: str, name: str = "") -> "Node":
        for node in self.children:
            if node.kind == kind and node.name == name:
                return node
        node = Node(kind, name)
        self.children.append(node)
        return node

    def at(self, path: BlockPath) -> "Node":
        node = self
        for kind, name in path:
            node = node.child(kind, name)
        return node


def render_item(item: Item, rename: Optional[str] = None) -> List[str]:
    lines = list(item.before)
    key = rename or item.key
    text = f"{key} = {item.value}"
    if item.trailing:
        text += "   " + item.trailing
    lines.append(text)
    return lines


def emit(node: Node, depth: int, out: List[str]) -> None:
    indent = "  " * depth
    for line in node.before:
        out.append(indent + line if line else "")
    out.append(indent + "#begin " + node.label())
    inner = "  " * (depth + 1)
    for line in node.lines:
        out.append(inner + line if line else "")
    for child in node.children:
        if node.lines or child is not node.children[0]:
            out.append("")
        emit(child, depth + 1, out)
    out.append(indent + "#end " + node.kind)
    for line in node.after:
        out.append(indent + line if line else "")


def emit_document(header: List[str], root: Node, trailing: List[str],
                  order: Sequence[Tuple[str, str]]) -> str:
    out: List[str] = list(header)
    rank = {key: index for index, key in enumerate(order)}
    children = sorted(root.children,
                      key=lambda n: rank.get((n.kind, n.name),
                                             rank.get((n.kind, "*"), len(rank))))
    for node in children:
        if out:
            out.append("")
        emit(node, 0, out)
    if trailing:
        out.append("")
        out.extend(trailing)
    return "\n".join(out).rstrip() + "\n"


# ---------------------------------------------------------------------------
# Layouts: old section -> new block path (the authoritative mapping)
# ---------------------------------------------------------------------------

def P(*parts: str) -> BlockPath:
    """``P("mesh", "refinement:solar")`` -> (("mesh",""),("refinement","solar"))."""
    out = []
    for part in parts:
        kind, _, name = part.partition(":")
        out.append((kind, name))
    return tuple(out)


SRCSEP3D_SECTIONS: Dict[str, BlockPath] = {
    "run": P("run"),
    "restart": P("run", "restart"),
    "domain": P("mesh", "domain"),
    "parker_spiral": P("mesh", "field_line"),
    "mesh": P("mesh"),
    "mesh.solar": P("mesh", "refinement:solar"),
    "mesh.tube": P("mesh", "refinement:tube"),
    "mesh.active_region": P("mesh", "active_region"),
    "memory": P("mesh", "memory"),
    "background": P("physics:background"),
    "background.parker": P("physics:background", "parker"),
    "turbulence": P("physics:turbulence"),
    "parallel_diffusion": P("physics:parallel_diffusion"),
    "swcme": P("physics:swcme"),
    "shock": P("physics:shock"),
    "transport": P("application:sep3d", "transport"),
    "population_control": P("application:sep3d", "population_control"),
    "storage": P("application:sep3d", "storage"),
    "source": P("source:shock"),
    "species": P("species"),
    "output": P("output"),
}
# [run] keys that describe the srcSEP3D application rather than the generic
# run control move into ``application: sep3d``.
SRCSEP3D_APPLICATION_RUN_KEYS = {
    "schema_version", "intent", "particle_mover", "transport",
    "background_cadence_steps", "injection_cadence_steps"}
SRCSEP3D_ORDER = [("run", ""), ("application", "sep3d"), ("mesh", ""),
                  ("physics", "background"), ("physics", "shock_front"),
                  ("physics", "swcme"), ("physics", "shock"),
                  ("physics", "turbulence"), ("physics", "parallel_diffusion"),
                  ("source", "shock"), ("species", ""), ("observer", "*"),
                  ("output", "")]

SRCSEP_SECTIONS: Dict[str, BlockPath] = {
    "run": P("run"),
    "injection": P("source:field_line"),
    "species": P("species"),
    "output": P("output"),
    "parker_spiral": P("mesh", "field_line"),
    "domain": P("mesh", "domain"),
    "mesh": P("mesh"),
    "mesh.solar": P("mesh", "refinement:solar"),
    "mesh.tube": P("mesh", "refinement:tube"),
    "background.parker": P("physics:background", "parker"),
    "swcme": P("physics:swcme"),
    "parallel_diffusion": P("physics:parallel_diffusion"),
    # Schema 2 had one unnamed observer section; schema 3/4 use observer.ID.
    "observer": P("observer"),
}
SRCSEP_APPLICATION_RUN_KEYS = {"schema_version", "particle_mover"}
SRCSEP_ORDER = [("run", ""), ("application", "sep"), ("mesh", ""),
                ("physics", "background"), ("physics", "swcme"),
                ("physics", "parallel_diffusion"), ("source", "field_line"),
                ("species", ""), ("observer", ""), ("observer", "*"),
                ("output", "")]


def convert_app(text: str, sections_map: Dict[str, BlockPath],
                app_keys: set, app_name: str,
                order: Sequence[Tuple[str, str]]) -> str:
    header, sections, trailing = lex_ini(text)
    root = Node("")
    includes: List[str] = []
    for section in sections:
        if section.name.startswith("observer."):
            path = P("observer:" + section.name.split(".", 1)[1])
        elif section.name in sections_map:
            path = sections_map[section.name]
        else:
            raise SystemExit(f"unknown section [{section.name}]; no owner known")
        node = root.at(path)
        node.before = section.before
        for item in section.items:
            target = node
            if section.name == "run" and item.key in app_keys:
                target = root.at(P("application:" + app_name))
            if section.name == "background" and item.key == "model_asset":
                # The reduced shock-front event is now its own physics block,
                # included from the converted event file.
                includes.extend(item.before)
                includes.append(f'#include "{item.value}"')
                continue
            target.lines.extend(render_item(item))
        node.lines.extend(section.after)
    if includes:
        holder = root.at(P("physics:background"))
        holder.after.extend([""] + includes)
    return emit_document(header, root, trailing, order)


# ---------------------------------------------------------------------------
# Modes
# ---------------------------------------------------------------------------

def convert_event(text: str, block: str = "shock_front",
                  verbatim: bool = False) -> str:
    """Wrap a ``key=value`` record into ``#begin physics: <block>``.

    ``verbatim`` keeps each assignment's original spelling (``key=value``
    without blanks) so that tests which edit the text by exact substring keep
    working; otherwise assignments are written as ``key = value``."""
    lines = text.splitlines()
    header: List[str] = []
    body: List[str] = []
    started = False
    for raw in lines:
        code, comment = split_comment(raw, "#")
        code = code.strip()
        if not code:
            rendered = bang(comment) if comment or raw.strip() else ""
            (body if started else header).append(rendered)
            continue
        started = True
        key, value = (part.strip() for part in code.split("=", 1))
        rendered = code if verbatim else f"{key} = {value}"
        if comment:
            rendered += "   " + bang(comment)
        body.append(rendered)
    header = strip_blank_edges(header)
    body = strip_blank_edges(body)
    out = list(header)
    if out:
        out.append("")
    out.append("#begin physics: " + block)
    out.extend(("  " + line) if line else "" for line in body)
    out.append("#end physics")
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------------------
# Shared SEP schema 6 (srcSEP and srcSEP3D, src/INPUT_FILE.md section 5)
# ---------------------------------------------------------------------------
#
# The modes sep6-from-srcsep3d and sep6-from-srcsep rewrite a generic-format
# deck of the retired content schemas (srcSEP3D 1-5 with
# '#begin application: sep3d', srcSEP 1-4 with '#begin application: sep') into
# the shared schema 6.  The rewrite works on a comment-preserving block tree:
# every comment and blank line stays where it was, except the comment lines
# directly attached to a key that schema 6 removes (they describe only that
# key).  The transformation implements decisions D1-D5 of
# src/models/sep_common/sep_input_layout.h:
#
#   D1  application: sep3d -> application: sep, schema_version -> 6; the
#       srcSEP3D key particle_mover is dropped (the family follows from
#       transport); srcSEP particle_mover becomes transport.
#   D2  physics: turbulence (srcSEP3D) moves into application: sep as
#       '#begin turbulence'; a srcSEP deck receives the complete turbulence
#       block with the values its command line used by default.
#   D4  with physics: swcme the seven ambient keys of physics.background.parker
#       are removed (SWCME owns them) and the SWCME source.* keys move to
#       source: shock in SI units (MeV -> J, mp -> kg) computed in exact
#       decimal arithmetic.
#   D5  srcSEP observers get kind = field-line-radius; energy_channels ->
#       energy_bins.
# Renamed keys: output.initialization_parker_line_tecplot_file and srcSEP's
# output product names -> initialization_{mesh,field_line,data}_tecplot_file;
# srcSEP source: field_line / macroparticles_per_step -> source: shock /
# samples_per_step; species.particle_weight -> macroparticle_weight.

from decimal import Decimal, getcontext

getcontext().prec = 40

ELEMENTARY_CHARGE = Decimal("1.602176634e-19")      # C (exact SI)
SWCME_PROTON_MASS_KG = Decimal("1.67262192595e-27")  # swcme_constants.hpp "mp"


class Entry:
    """One line of a block body: kind is 'blank', 'comment', 'assign',
    'include' or 'block'."""

    def __init__(self, kind: str, text: str = "", key: str = "",
                 value: str = "", comment: str = "", block=None) -> None:
        self.kind = kind
        self.text = text
        self.key = key
        self.value = value
        self.comment = comment
        self.block = block


class Block:
    def __init__(self, kind: str, name: str, header_comment: str = "",
                 end_comment: str = "") -> None:
        self.kind = kind
        self.name = name
        self.header_comment = header_comment
        self.end_comment = end_comment
        self.entries: List[Entry] = []

    # -- navigation -----------------------------------------------------------
    def child(self, kind: str, name: str = "") -> Optional["Block"]:
        for entry in self.entries:
            if entry.kind == "block" and entry.block.kind == kind and \
                    entry.block.name == name:
                return entry.block
        return None

    def get(self, key: str) -> Optional[Entry]:
        for entry in self.entries:
            if entry.kind == "assign" and entry.key == key:
                return entry
        return None

    def keys(self) -> List[str]:
        return [e.key for e in self.entries if e.kind == "assign"]

    # -- edits ------------------------------------------------------------------
    def _attached_start(self, index: int) -> int:
        """First index of the comment lines directly above entry ``index``
        (no blank line in between); they describe only that entry."""
        start = index
        while start > 0 and self.entries[start - 1].kind == "comment":
            start -= 1
        return start

    def remove(self, key: str) -> Optional[Entry]:
        for index, entry in enumerate(self.entries):
            if entry.kind == "assign" and entry.key == key:
                start = self._attached_start(index)
                del self.entries[start:index + 1]
                return entry
        return None

    def remove_child(self, kind: str, name: str = "") -> Optional[List[Entry]]:
        """Detach a child block with its attached comment lines."""
        for index, entry in enumerate(self.entries):
            if entry.kind == "block" and entry.block.kind == kind and \
                    entry.block.name == name:
                start = self._attached_start(index)
                taken = self.entries[start:index + 1]
                del self.entries[start:index + 1]
                # Collapse a now-doubled blank line.
                if start < len(self.entries) and start > 0 and \
                        self.entries[start].kind == "blank" and \
                        self.entries[start - 1].kind == "blank":
                    del self.entries[start]
                return taken
        return None

    def rename(self, old: str, new: str) -> None:
        entry = self.get(old)
        if entry is not None:
            entry.key = new

    def set(self, key: str, value: str) -> None:
        entry = self.get(key)
        if entry is None:
            raise SystemExit(f"internal: key {key} absent in {self.kind}")
        entry.value = value

    def append(self, entries: List[Entry]) -> None:
        # Insert before trailing blank lines so the block stays tidy.
        position = len(self.entries)
        while position > 0 and self.entries[position - 1].kind == "blank":
            position -= 1
        self.entries[position:position] = entries

    def prepend_assign(self, key: str, value: str, comments: Sequence[str] = ()) -> None:
        position = 0
        new = [Entry("comment", text=c) for c in comments]
        new.append(Entry("assign", key=key, value=value))
        self.entries[position:position] = new


def assign(key: str, value: str, comment: str = "") -> Entry:
    return Entry("assign", key=key, value=value, comment=comment)


def comment_lines(*texts: str) -> List[Entry]:
    return [Entry("comment", text=t) for t in texts]


def parse_blocks(text: str) -> Block:
    """Parse generic-format text into a comment-preserving block tree."""
    root = Block("", "")
    stack = [root]
    for number, raw in enumerate(text.splitlines(), 1):
        line = raw.rstrip()
        stripped = line.strip()
        if not stripped:
            stack[-1].entries.append(Entry("blank"))
            continue
        if stripped.startswith("!"):
            stack[-1].entries.append(Entry("comment", text=stripped[1:]))
            continue
        code, _, comment = line.partition("!")
        code = code.strip()
        if code.endswith("\\"):
            raise SystemExit(f"line {number}: continuation lines are not "
                             "supported by the schema-6 converter")
        lowered = code.lower()
        if lowered.startswith("#begin"):
            spec = code[len("#begin"):].strip()
            kind, _, name = spec.partition(":")
            block = Block(kind.strip().lower(), name.strip(), comment.rstrip())
            stack[-1].entries.append(Entry("block", block=block))
            stack.append(block)
            continue
        if lowered.startswith("#end"):
            if len(stack) == 1:
                raise SystemExit(f"line {number}: unbalanced #end")
            stack[-1].end_comment = comment.rstrip()
            stack.pop()
            continue
        if lowered.startswith("#include"):
            stack[-1].entries.append(Entry("include", text=line.strip()))
            continue
        if "=" not in code:
            raise SystemExit(f"line {number}: not an assignment: {raw}")
        key, value = (part.strip() for part in code.split("=", 1))
        stack[-1].entries.append(Entry("assign", key=key, value=value,
                                       comment=comment.rstrip()))
    if len(stack) != 1:
        raise SystemExit("unterminated block at end of file")
    return root


def render_blocks(root: Block) -> str:
    out: List[str] = []

    def emit_entries(block: Block, depth: int) -> None:
        pad = "  " * depth
        for entry in block.entries:
            if entry.kind == "blank":
                out.append("")
            elif entry.kind == "comment":
                out.append(pad + "!" + entry.text)
            elif entry.kind == "include":
                out.append(pad + entry.text)
            elif entry.kind == "assign":
                line = f"{pad}{entry.key} = {entry.value}"
                if entry.comment:
                    line += "   !" + entry.comment
                out.append(line)
            else:
                child = entry.block
                head = f"{pad}#begin {child.kind}" + \
                    (f": {child.name}" if child.name else "")
                if child.header_comment:
                    head += "   !" + child.header_comment
                out.append(head)
                emit_entries(child, depth + 1)
                tail = f"{pad}#end {child.kind}"
                if child.end_comment:
                    tail += "   !" + child.end_comment
                out.append(tail)

    emit_entries(root, 0)
    while out and out[-1] == "":
        out.pop()
    return "\n".join(out) + "\n"


def decimal_text(value: Decimal) -> str:
    """Shortest exponent notation of an exact decimal (no trailing zeros)."""
    normalized = value.normalize()
    sign, digits, exponent = normalized.as_tuple()
    mantissa = "".join(str(d) for d in digits)
    scientific = exponent + len(mantissa) - 1
    text = mantissa[0] + ("." + mantissa[1:] if len(mantissa) > 1 else "")
    return ("-" if sign else "") + text + (f"e{scientific}" if scientific else "")


def energy_to_joule(value: str) -> str:
    number, _, unit = value.strip().partition(" ")
    unit = unit.strip().lower()
    scale = {"mev": Decimal("1e6"), "kev": Decimal("1e3"), "j": None}
    if unit not in scale:
        raise SystemExit(f"cannot convert SWCME energy '{value}' to joules")
    if scale[unit] is None:
        return decimal_text(Decimal(number))
    return decimal_text(Decimal(number) * scale[unit] * ELEMENTARY_CHARGE)


def mass_to_kg(value: str) -> str:
    number, _, unit = value.strip().partition(" ")
    unit = unit.strip().lower()
    if unit == "kg":
        return decimal_text(Decimal(number))
    if unit == "mp":
        return decimal_text(Decimal(number) * SWCME_PROTON_MASS_KG)
    raise SystemExit(f"cannot convert SWCME mass '{value}' to kg")


def plain_number(value: str, unit: str) -> str:
    number, _, written = value.strip().partition(" ")
    if written.strip().lower() not in ("", unit):
        raise SystemExit(f"unexpected unit in '{value}' (expected {unit})")
    return number


PARKER_SWCME_OWNED = ("reference_radius_m", "radial_field_at_reference_t",
                      "solar_rotation_rate_rad_per_s", "solar_wind_speed_m_per_s",
                      "magnetic_polarity", "number_density_at_one_au_m3",
                      "temperature_k")


def take_swcme_source(swcme: Block) -> Dict[str, str]:
    """Remove the SWCME source.* keys (D4) and return their values."""
    values: Dict[str, str] = {}
    for key in ("source.particle_mass", "source.charge_number",
                "source.energy_min", "source.energy_max",
                "source.reference_energy", "source.injection_efficiency",
                "source.normalization", "source.reference_intensity_si"):
        entry = swcme.remove(key)
        if entry is None:
            raise SystemExit(f"physics: swcme lacks {key}")
        values[key] = entry.value
    return values


def spectrum_entries(source: Dict[str, str]) -> List[Entry]:
    normalization = source["source.normalization"].replace("_", "-")
    return comment_lines(
        " Reference spectrum of the canonical SWCME source (written once here;",
        " srcSEP and srcSEP3D synthesize the SWCME source.* keys from it).") + [
        assign("spectrum_reference_energy_j",
               energy_to_joule(source["source.reference_energy"])),
        assign("spectrum_particle_mass_kg",
               mass_to_kg(source["source.particle_mass"])),
        assign("spectrum_charge_number", source["source.charge_number"]),
        assign("spectrum_normalization", normalization),
        assign("spectrum_reference_intensity_per_m2_s_sr_j",
               plain_number(source["source.reference_intensity_si"], "si"))]


def top_block(root: Block, kind: str, name: str = "") -> Block:
    block = root.child(kind, name)
    if block is None:
        raise SystemExit(f"missing block '{kind}{': ' + name if name else ''}'")
    return block


def insert_top_block_after(root: Block, after: Block, block: Block,
                           comments: Sequence[str] = ()) -> None:
    for index, entry in enumerate(root.entries):
        if entry.kind == "block" and entry.block is after:
            new = [Entry("blank")] + comment_lines(*comments) + \
                [Entry("block", block=block)]
            root.entries[index + 1:index + 1] = new
            return
    raise SystemExit("internal: anchor block not found")


def convert_sep6_from_srcsep3d(text: str) -> str:
    root = parse_blocks(text)
    application = top_block(root, "application", "sep3d")
    application.name = "sep"
    application.set("schema_version", "6")
    application.remove("particle_mover")
    transport = application.get("transport")
    if transport is not None:
        transport.value = {"parker3d": "parker", "focused3d":
                           "focused-diffusion"}.get(transport.value.lower(),
                                                    transport.value)
    numerics = application.child("transport")
    if numerics is not None:
        for key in ("perpendicular_diffusion", "drifts"):
            entry = numerics.get(key)
            if entry is not None and entry.value.lower() == "false":
                entry.value = "none"
    # D2: the prescribed turbulence block moves into application: sep.
    taken = root.remove_child("physics", "turbulence")
    if taken is not None:
        turbulence = taken[-1].block
        turbulence.kind, turbulence.name = "turbulence", ""
        application.append([Entry("blank")] + taken[:-1] +
                           [Entry("block", block=turbulence)])
    background = top_block(root, "physics", "background")
    background.remove("external_script")
    swcme = root.child("physics", "swcme")
    if swcme is not None:
        parker = background.child("parker")
        if parker is not None:
            # The retired schema-1/2 density spelling is SWCME-owned too.
            for key in PARKER_SWCME_OWNED + ("number_density_at_reference_m3",):
                parker.remove(key)
            parker.prepend_assign(
                "validity_cadence_s", parker.remove("validity_cadence_s").value,
                [" The ambient solar wind (speed, rotation, density, temperature,",
                 " field and polarity) is owned by physics: swcme; only the",
                 " snapshot validity cadence of the analytic provider is set here."])
        source_values = take_swcme_source(swcme)
        source = top_block(root, "source", "shock")
        source.append(spectrum_entries(source_values))
    # The particle source model is explicit (user decision 2026-10-12): the
    # reduced shock front is the only model.  A converted deck that still
    # holds physics: swcme is rejected by srcSEP3D as "not ready".
    source = root.child("source", "shock")
    if source is not None and source.get("model") is None:
        source.prepend_assign(
            "model", "reduced-shock-front",
            [" Particle source model: the reduced shock front of physics:",
             " shock_front (the only model; SWCME is not ready)."])
    output = root.child("output")
    if output is not None:
        output.rename("initialization_parker_line_tecplot_file",
                      "initialization_field_line_tecplot_file")
    return render_blocks(root)


# srcSEP turbulence values used before schema 6 when no command-line option
# was given (SEP::Turbulence::Configuration defaults and the srcSEP CLI
# defaults; srcSEP/util/sep_turbulence_core.h, srcSEP/util/sep_cli.h).
SRCSEP_TURBULENCE_DEFAULTS = [
    (" Self-consistent Alfven-wave turbulence on the field line", None),
    ("source", "self-consistent-integrated"),
    ("representation", "integrated"),
    ("coupling", "streaming-energy-exchange"),
    ("advection", "true"),
    ("reflection", "true"),
    ("cascade", "true"),
    ("shock_injection", "true"),
    (" Boundary policies; the values are J (energy/reservoir) or W (flux)", None),
    ("inner_boundary", "fixed-reservoir"),
    ("inner_boundary_value_si", "0"),
    ("outer_boundary", "fixed-reservoir"),
    ("outer_boundary_value_si", "0"),
    (" Reflection/cascade closures and the spectral grid", None),
    ("reflection_coefficient", "0.6"),
    ("cascade_coefficient", "0.8"),
    ("perpendicular_correlation_length_m", "1e7"),
    ("electron_heating_fraction", "0.3"),
    ("spectral_k_min_per_m", "1e-10"),
    ("spectral_k_max_per_m", "1e-2"),
    ("spectral_bins", "128"),
    (" Operator accuracy controls", None),
    ("cfl_safety", "0.8"),
    ("operator_accuracy_safety", "0.25"),
    ("maximum_source_fraction", "0.25"),
    ("maximum_cascade_fraction", "0.2"),
    ("minimum_substep_s", "1e-12"),
    ("maximum_substeps", "1000000"),
    ("conservation_relative_tolerance", "1e-12"),
    (" Wave energy deposited by the shock and its E+ share", None),
    ("shock_source_efficiency", "0.02"),
    ("shock_source_plus_fraction", "0.5"),
    (" Spectrum diagnostic cadence in main-loop steps (0 = off)", None),
    ("diagnostic_cadence_steps", "100"),
]

SRCSEP_MOVER_TO_TRANSPORT = {"parker": "parker", "fte-dmumu": "focused-diffusion",
                             "fte-mfp": "focused-scattering"}


def convert_sep6_from_srcsep(text: str) -> str:
    root = parse_blocks(text)
    application = top_block(root, "application", "sep")
    application.set("schema_version", "6")
    mover = application.remove("particle_mover")
    # Without a file mover srcSEP ran its built-in default fte-dmumu.
    transport = SRCSEP_MOVER_TO_TRANSPORT[mover.value.lower()] if mover \
        else "focused-diffusion"
    application.append(comment_lines(
        " parker | focused-diffusion | focused-scattering (srcSEP movers parker,",
        " fte-dmumu, fte-mfp); --particle-mover overrides it.") +
        [assign("transport", transport)])
    turbulence = Block("turbulence", "")
    for key, value in SRCSEP_TURBULENCE_DEFAULTS:
        turbulence.entries.append(Entry("comment", text=key) if value is None
                                  else assign(key, value))
    application.append([Entry("blank")] + comment_lines(
        " The only application-specific block: srcSEP's turbulence model.",
        " Every key has a command-line option that overrides it.") +
        [Entry("block", block=turbulence)])
    background = top_block(root, "physics", "background")
    background.remove_child("parker")
    background.entries[0:0] = comment_lines(
        " srcSEP's ambient plasma is the standalone SWCME state along the",
        " Parker field line; physics: swcme owns the solar-wind values.") + \
        [assign("provider", "swcme")]
    swcme = top_block(root, "physics", "swcme")
    source_values = take_swcme_source(swcme)
    shock = Block("physics", "shock")
    shock.entries.append(assign("authority", "swcme"))
    insert_top_block_after(root, background, shock,
                           [" The shock is the canonical SWCME front."])
    old_source = top_block(root, "source", "field_line")
    old_source.name = "shock"
    old_source.rename("macroparticles_per_step", "samples_per_step")
    old_source.entries[0:0] = [assign("enabled", "true")]
    old_source.append([Entry("blank")] + comment_lines(
        " Kinetic-energy range and injection efficiency of the source.") + [
        assign("minimum_energy_j", energy_to_joule(source_values["source.energy_min"])),
        assign("maximum_energy_j", energy_to_joule(source_values["source.energy_max"])),
        assign("injection_efficiency", source_values["source.injection_efficiency"])] +
        spectrum_entries(source_values))
    species = top_block(root, "species")
    species.rename("particle_weight", "macroparticle_weight")
    for entry in root.entries:
        if entry.kind == "block" and entry.block.kind == "observer" and \
                entry.block.name:
            observer = entry.block
            observer.rename("energy_channels", "energy_bins")
            observer.prepend_assign("kind", "field-line-radius")
    output = top_block(root, "output")
    output.rename("mesh_tecplot_file", "initialization_mesh_tecplot_file")
    output.rename("field_line_tecplot_file", "initialization_field_line_tecplot_file")
    output.rename("data_tecplot_file", "initialization_data_tecplot_file")
    return render_blocks(root)


def convert(mode: str, text: str) -> str:
    if mode == "event":
        return convert_event(text)
    if mode == "cme-event":
        return convert_event(text, "cme_background", verbatim=True)
    # The INI modes chain into the schema-6 rewrite: no reader accepts the
    # intermediate block-format schemas any more.
    if mode == "srcsep3d":
        return convert_sep6_from_srcsep3d(convert_app(
            text, SRCSEP3D_SECTIONS, SRCSEP3D_APPLICATION_RUN_KEYS,
            "sep3d", SRCSEP3D_ORDER))
    if mode == "srcsep":
        return convert_sep6_from_srcsep(convert_app(
            text, SRCSEP_SECTIONS, SRCSEP_APPLICATION_RUN_KEYS,
            "sep", SRCSEP_ORDER))
    if mode == "sep6-from-srcsep3d":
        return convert_sep6_from_srcsep3d(text)
    if mode == "sep6-from-srcsep":
        return convert_sep6_from_srcsep(text)
    raise SystemExit("unknown mode " + mode)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("mode", choices=["event", "cme-event", "srcsep3d",
                                         "srcsep", "sep6-from-srcsep3d",
                                         "sep6-from-srcsep"])
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument("-o", "--output", type=Path,
                        help="output file (single input only); default stdout")
    parser.add_argument("--in-place", action="store_true",
                        help="rewrite each input file")
    args = parser.parse_args(argv)
    for path in args.files:
        result = convert(args.mode, path.read_text(encoding="utf-8"))
        if args.in_place:
            path.write_text(result, encoding="utf-8")
        elif args.output:
            args.output.write_text(result, encoding="utf-8")
        else:
            sys.stdout.write(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
