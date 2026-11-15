"""Static report pages generated from the results store.

    <www>/index.html                         test families x columns (machine + compiler)
    <www>/family/<family>/<column>.html      tests of one family in one column
    <www>/assets/report.css, report.js

Pages are plain HTML that work without a server; report.js only adds the status filter buttons.
"""

from __future__ import annotations

import datetime as dt
import html
import json
import os
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .records import parse_time
from .registry import OTHER_GROUP, Registry
from .store import Store, night_of

COMMIT_URL = "https://github.com/SWMFsoftware/AMPS/commit/"
ASSETS = Path(__file__).resolve().parent / "assets"
RANK = {"FAILED": 0, "SKIPPED": 1, "PASS": 2}
_SVG = 'aria-hidden="true" viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-linecap="round" stroke-linejoin="round"'
ICON = {
    "PASS": '<svg {} stroke-width="2.4"><path d="M3.5 8.5l3 3 6-7"/></svg>'.format(_SVG),
    "FAILED": '<svg {} stroke-width="2.4"><path d="M4.5 4.5l7 7M11.5 4.5l-7 7"/></svg>'.format(_SVG),
    "SKIPPED": '<svg {} stroke-width="2.4"><path d="M4 8h8"/></svg>'.format(_SVG),
    "WARN": '<svg {} stroke-width="1.8"><path d="M8 1.8l6.6 11.6H1.4z"/><path d="M8 6.3v3.2M8 11.6v.2"/></svg>'.format(_SVG),
    "BACK": '<svg {} stroke-width="2"><path d="M10 3.5L5.5 8l4.5 4.5"/></svg>'.format(_SVG),
}
_PILL = {"PASS": "st-pass", "FAILED": "st-fail", "SKIPPED": "st-skip"}
_FILL = {"PASS": "fill-pass", "FAILED": "fill-fail", "SKIPPED": "fill-skip"}


def esc(value: Any) -> str:
    return html.escape(str(value), quote=True)


def pill(status: str) -> str:
    return '<span class="st {}">{}{}</span>'.format(_PILL[status], ICON[status], status)


# ----------------------------------------------------------------------------- data model
class Column:
    www = None
    """Latest (and previous) run of one machine + compiler configuration."""

    def __init__(self, store: Store, column_id: str, now: dt.datetime, tz: Optional[dt.tzinfo], stale_hours: float):
        runs = store.latest_runs(column_id, 2)
        self.id = column_id
        self.run = runs[0]
        self.previous = runs[1] if len(runs) > 1 else None
        self.machine = self.run["machine"]
        self.label = self.run["compiler"] or self.run["config"].upper()
        # legacy imports only know when the results reached the server, not when the run ended
        self.finished_verb = ("Results collected" if self.run["producer"].startswith("import-legacy")
                              else "Run finished")
        self.tz = tz
        self.finished = _local(self.run["finished"], tz)
        age = now - parse_time(self.run["finished"])
        self.stale = age > dt.timedelta(hours=stale_hours)
        self.age_days = max(0, age.days)
        self.results = {(r["family"], r["test"]): r for r in store.results(self.run["run_id"])}
        self.previous_results = ({(r["family"], r["test"]): r for r in store.results(self.previous["run_id"])}
                                 if self.previous else {})
        self.last_pass = store.last_passes(column_id)

    def family_rows(self, family: str) -> List[Dict[str, Any]]:
        rows = [r for (f, _), r in self.results.items() if f == family]
        return sorted(rows, key=lambda r: (RANK[r["status"]], r["test"].lower()))

    def families(self) -> List[str]:
        return sorted({f for f, _ in self.results})


class Waiting:
    www = None
    """A column declared by an enabled machine file that has no run in the store yet."""

    run = None
    stale = False

    def __init__(self, declared: Dict[str, str]):
        self.id = declared["column"]
        self.machine = declared["machine"]
        self.label = declared["compiler"]
        self.toolchain = declared.get("toolchain", "")
        self.results: Dict[Tuple[str, str], Dict[str, Any]] = {}
        self.previous_results: Dict[Tuple[str, str], Dict[str, Any]] = {}

    def family_rows(self, family: str) -> List[Dict[str, Any]]:
        return []

    def families(self) -> List[str]:
        return []


def _local(utc: str, tz: Optional[dt.tzinfo]) -> dt.datetime:
    return parse_time(utc).astimezone(tz) if tz else parse_time(utc).astimezone()


def _when(t: dt.datetime) -> str:
    return "{} {}".format(t.day, t.strftime("%b %Y, %H:%M"))


def _day(t: dt.datetime) -> str:
    return "{} {}".format(t.day, t.strftime("%b %Y"))


def tally(rows: List[Dict[str, Any]]) -> Tuple[int, int, int, str]:
    """(pass, skip, fail, family status); suites with case counts contribute their cases."""
    p = s = x = 0
    for row in rows:
        n = row["cases_total"] or 1
        if row["status"] == "PASS":
            p += n
        elif row["status"] == "SKIPPED":
            s += n
        else:
            failed = min(row["cases_failed"] or n, n)
            x += failed
            p += n - failed
    status = "FAILED" if x else ("SKIPPED" if p == 0 else "PASS")
    return p, s, x, status


def counts_text(p: int, s: int, x: int) -> str:
    parts = ["{} pass".format(p)]
    if s:
        parts.append("{} skip".format(s))
    if x:
        parts.append("{} fail".format(x))
    return " · ".join(parts)


def commit_html(sha: str, inferred: bool) -> str:
    if not sha:
        return '<span class="sub">—</span>'
    title = ' title="Inferred from the checkout time of a legacy nightly run"' if inferred else ""
    mark = '<span class="inferred"{}>≈ </span>'.format(title) if inferred else ""
    return '{}<a class="mono" href="{}{}"{}>{}</a>'.format(mark, COMMIT_URL, esc(sha), title, esc(sha[:7]))


# ----------------------------------------------------------------------------- page pieces
def _page(title: str, body: str, depth: int, generated: str) -> str:
    up = "../" * depth
    return (
        "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        "<title>{title}</title>\n"
        "<link rel=\"preconnect\" href=\"https://fonts.googleapis.com\">\n"
        "<link rel=\"stylesheet\" href=\"https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700"
        "&amp;family=JetBrains+Mono:wght@400;500;600&amp;display=swap\">\n"
        "<link rel=\"stylesheet\" href=\"{up}assets/report.css\">\n"
        "<script src=\"{up}assets/report.js\" defer></script>\n"
        "</head>\n<body>\n<div class=\"wrap\">\n{body}\n"
        "<footer class=\"foot\">Generated {generated} by <span class=\"mono\">utility/nightly/amps-nightly report</span>.</footer>\n"
        "</div>\n</body>\n</html>\n"
    ).format(title=esc(title), up=up, body=body, generated=esc(generated))


def _column_header(col: Column, now_local: dt.datetime) -> str:
    run = col.run
    if run is None:
        tool = esc(col.label) + (" · " + esc(col.toolchain) if col.toolchain else "")
        return ('<div class="hcell" role="columnheader"><span class="mname">{}</span><span>{}</span>'
                '<span>No results yet</span></div>').format(esc(col.machine), tool)
    commit = commit_html(run["commit_sha"], bool(run["commit_inferred"]))
    warn = ""
    if col.stale:
        warn = '<span class="warn">{}No upload for {} day{}</span>'.format(
            ICON["WARN"], col.age_days, "" if col.age_days == 1 else "s")
    tool = esc(col.label) + (" · " + esc(run["toolchain"]) if run["toolchain"] else "")
    return ('<div class="hcell" role="columnheader"><span class="mname">{}</span><span>{}</span>'
            '<span>{} {} · {}</span>{}{}</div>').format(esc(col.machine), tool, col.finished_verb,
                                                        esc(_when(col.finished)), commit, warn, _night_link(col, ""))


# ----------------------------------------------------------------------------- pages
def _index(registry: Registry, columns: List[Column], shown: List[str], unseen: List[str],
           now: dt.datetime, now_local: dt.datetime, stale_hours: float, hidden: List[str]) -> str:
    """columns may include Waiting placeholders; statistics and changes use the columns that have a run."""
    waiting = [c for c in columns if c.run is None]
    all_columns, columns = columns, [c for c in columns if c.run is not None]
    new_failures, fixed = [], []
    for col in columns:
        for key, row in col.results.items():
            before = col.previous_results.get(key)
            if not before:
                continue
            if before["status"] == "PASS" and row["status"] == "FAILED":
                new_failures.append((col, row))
            elif before["status"] == "FAILED" and row["status"] == "PASS":
                fixed.append((col, row))
    failing_families = sum(1 for f in shown if any(tally(c.family_rows(f))[3] == "FAILED" for c in columns))
    total = sum(sum(tally(c.family_rows(f))[:3]) for c in columns for f in shown)
    fresh = [c for c in columns if not c.stale]
    stale_note = "; ".join(["{} · {}: last upload {}".format(c.machine, c.label, _day(c.finished))
                            for c in columns if c.stale] +
                           ["{} · {}: no results yet".format(c.machine, c.label) for c in waiting]) \
        or "every column uploaded within {} h".format(int(stale_hours))
    stats = [
        ("Families with a failure", "{} of {}".format(failing_families, len(shown)), "failing in at least one column"),
        ("Individual tests", "{:,}".format(total), "in the latest run of each column"),
        ("New failures", str(len(new_failures)), "passed in the previous run, fail now"),
        ("Columns up to date", "{} of {}".format(len(fresh), len(all_columns)), stale_note),
    ]
    commits = {(c.run["commit_sha"], c.run["commit_inferred"]) for c in columns}
    lead = "Latest results of {} machine and compiler configuration{}".format(len(columns), "" if len(columns) == 1 else "s")
    if len(commits) == 1:
        sha, inferred = next(iter(commits))
        if sha:
            lead += " · branch <span class=\"mono\">master</span> at " + commit_html(sha, bool(inferred))
    out = ['<div class="page">',
           '<header class="top"><div><p class="eyebrow">AMPS · SWMFsoftware/AMPS</p><h1>Nightly test results</h1>'
           '<p class="lead">{}</p></div></header>'.format(lead)]
    out.append('<div class="stats">' + "".join(
        '<div class="card stat"><span class="stat-l">{}</span><span class="stat-v">{}</span><span class="stat-n">{}</span></div>'
        .format(esc(a), esc(b), esc(c)) for a, b, c in stats) + "</div>")

    changes = [(c, r, True) for c, r in new_failures] + [(c, r, False) for c, r in fixed]
    if changes:
        items = []
        for col, row, is_new in changes[:12]:
            items.append('<a class="chg" href="family/{f}/{c}.html"><span class="st {cls}">{kind}</span>'
                         '<span class="chg-t"><span class="mono">{t}</span><span class="sub">{fam} · {m} · {l}</span></span></a>'.format(
                             f=esc(row["family"]), c=esc(col.id), cls="st-fail" if is_new else "st-pass",
                             kind="NEW FAILURE" if is_new else "FIXED", t=esc(row["test"]),
                             fam=esc(registry.describe(row["family"])["title"]), m=esc(col.machine), l=esc(col.label)))
        more = '<span class="more">and {} more</span>'.format(len(changes) - 12) if len(changes) > 12 else ""
        out.append('<section class="card changes"><h2>Changes since the previous run</h2><div class="chgs">{}{}</div></section>'
                   .format("".join(items), more))

    out.append('<div class="legend"><span>{}no failures (some tests may be skipped)</span>'
               '<span>{}at least one test failed</span><span>{}every test skipped</span>'
               '<span class="sub">Click a family or a cell to see its tests.</span></div>'.format(
                   pill("PASS"), pill("FAILED"), pill("SKIPPED")))

    width = 320 + 240 * len(all_columns)
    out.append('<section class="card scroll"><div role="table" aria-label="Test families by machine and compiler" '
               'style="min-width:{}px">'.format(width))
    out.append('<div class="hrow" role="row"><div class="hcell first" role="columnheader">Test family</div>{}</div>'.format(
        "".join(_column_header(c, now_local) for c in all_columns)))
    for group in registry.groups:
        families = [f for f in shown if registry.describe(f)["group"] == group["id"]]
        if not families:
            continue
        out.append('<div class="grow" role="row"><span class="glabel" role="rowheader">{}</span>'
                   '<span class="gsrc mono">{}</span></div>'.format(esc(group["title"]), esc(group["source"])))
        for fid in families:
            fam = registry.describe(fid)
            present = [c for c in columns if c.family_rows(fid)]
            tests = max(len(c.family_rows(fid)) for c in present)
            meta = "{} test{}".format(tests, "" if tests == 1 else "s") + (" · " + fam["runner"] if fam["runner"] else "")
            first = present[0].id
            cells = []
            for col in all_columns:
                rows = col.family_rows(fid)
                if col.run is None:
                    cells.append('<div class="cell" role="cell"><span class="na">— no results yet</span></div>')
                    continue
                if not rows:
                    cells.append('<div class="cell" role="cell"><span class="na">— not run here</span></div>')
                    continue
                p, s, x, status = tally(rows)
                cells.append('<div class="cell" role="cell"><a class="cellbtn" href="family/{f}/{c}.html" '
                             'aria-label="{aria}">{pill}<span class="counts">{counts}</span></a></div>'.format(
                                 f=esc(fid), c=esc(col.id), pill=pill(status), counts=esc(counts_text(p, s, x)),
                                 aria=esc("{} on {} {}: {}, {}".format(fam["title"], col.machine, col.label, status,
                                                                      counts_text(p, s, x)))))
            out.append('<div class="row" role="row"><div class="cell first" role="rowheader"><a class="fambtn" '
                       'href="family/{f}/{c}.html"><span class="famname">{t}</span><span class="fammeta">{m}</span></a></div>{cells}</div>'
                       .format(f=esc(fid), c=esc(first), t=esc(fam["title"]), m=esc(meta), cells="".join(cells)))
    out.append("</div></section>")
    if unseen:
        out.append('<p class="sub">No results yet for: {}.</p>'.format(
            esc(", ".join(registry.describe(f)["title"] for f in unseen))))
    if hidden:
        out.append('<p class="sub">Not shown, because no enabled machine file declares them: {}.</p>'.format(
            esc(", ".join(hidden))))
    out.append("</div>")
    return "\n".join(out)


def _family_page(registry: Registry, fid: str, col: Column, columns: List[Column]) -> str:
    fam = registry.describe(fid)
    group = registry.group(fam["group"])
    rows = col.family_rows(fid)
    out = ['<div class="page tight">',
           '<div class="bar"><a class="back" href="../../index.html">{}All test families</a></div>'.format(ICON["BACK"]),
           '<header class="fhead"><p class="eyebrow">{}</p><h1>{}</h1><p class="src"><span class="mono">{}</span></p></header>'
           .format(esc(group["title"]), esc(fam["title"]), esc(fam["runner"] or group["source"]))]
    tabs = []
    for other in columns:
        other_rows = other.family_rows(fid)
        current = ' aria-current="true"' if other.id == col.id else ""
        if other_rows:
            p, s, x, status = tally(other_rows)
            summary = '<span class="tab-sum">{}<span class="counts">{}</span></span>'.format(pill(status), esc(counts_text(p, s, x)))
        else:
            summary = '<span class="counts na-text">not run here</span>'
        tabs.append('<a class="tab" href="{c}.html"{cur}><span class="tab-name">{m}</span><span class="tab-cfg">{l}</span>{s}</a>'
                    .format(c=esc(other.id), cur=current, m=esc(other.machine), l=esc(other.label), s=summary))
    run = col.run
    warn = ('<span class="warn">{}No upload for {} day{}</span>'.format(ICON["WARN"], col.age_days, "" if col.age_days == 1 else "s")
            if col.stale else "")
    out.append('<div><nav class="tabs" aria-label="Machine and compiler">{}</nav>'
               '<div class="card runline"><span><strong>{}</strong> · {}</span><span>{} {}</span>'
               '<span>Commit {}</span>{}{}</div></div>'.format(
                   "".join(tabs), esc(col.machine), esc(col_label(col)), col.finished_verb, esc(_when(col.finished)),
                   commit_html(run["commit_sha"], bool(run["commit_inferred"])), warn, _night_link(col, "../../")))
    if not rows:
        out.append('<div class="card notcfg">This family has no results in the latest run of <strong>{} · {}</strong>.</div>'
                   .format(esc(col.machine), esc(col.label)))
        out.append("</div>")
        return "\n".join(out)

    counts = {s: sum(1 for r in rows if r["status"] == s) for s in RANK}
    out.append('<nav class="chips" aria-label="Filter by status" data-filters hidden>' + "".join(
        '<button type="button" class="chip" data-filter="{}" aria-pressed="{}">{} <span class="n">{}</span></button>'.format(
            key, "true" if key == "all" else "false", label, n)
        for key, label, n in (("all", "All", len(rows)), ("FAILED", "Failed", counts["FAILED"]),
                              ("SKIPPED", "Skipped", counts["SKIPPED"]), ("PASS", "Passed", counts["PASS"]))) + "</nav>")
    out.append('<section class="card scroll"><div role="table" aria-label="Tests" style="min-width:820px">'
               '<div class="thead" role="row"><div class="c-name" role="columnheader">Test name</div>'
               '<div class="c-status" role="columnheader">Current status</div>'
               '<div class="c-date" role="columnheader">Last date passed</div>'
               '<div class="c-commit" role="columnheader">Git commit when last passed</div></div>')
    latest_day = col.finished.date()
    for row in rows:
        sub = ""
        n = row["cases_total"]
        if n:
            unit = "case" if n == 1 else "cases"
            sub = ("{} of {} cases failed".format(row["cases_failed"] or n, n) if row["status"] == "FAILED"
                   else "{} {} skipped".format(n, unit) if row["status"] == "SKIPPED" else "{} {}".format(n, unit))
        detail = ""
        if row["status"] != "PASS":
            reason = row["reason"] or ("failed at the {} stage".format(row["stage"]) if row["stage"] else "")
            if reason:
                detail += '<span class="sub">{}</span>'.format(esc(reason))
            links = []
            if row["log"]:
                links.append('<a href="../../{}/{}">log</a>'.format(esc(run["path"]), esc(log_href(col.www, run["path"], row["log"]))))
            if links:
                detail += '<span class="links">{}</span>'.format("".join(links))
        last = col.last_pass.get((row["family"], row["test"]))
        if last:
            when = _local(last["finished"], col.tz)
            if last["run_id"] == run["run_id"]:
                ago = "this run"
            else:
                days = (latest_day - when.date()).days
                ago = "same day, earlier run" if days <= 0 else "{} day{} before this run".format(days, "" if days == 1 else "s")
            date_html = '<span>{}</span><span class="sub">{}</span>'.format(esc(_day(when)), esc(ago))
            commit = commit_html(last["commit_sha"], bool(last["commit_inferred"]))
        else:
            date_html = '<span>never</span><span class="sub">no passing run stored</span>'
            commit = '<span class="sub">—</span>'
        out.append('<div class="trow{fail}" role="row" data-status="{st}"><div class="c-name" role="rowheader">'
                   '<span class="mono tname">{name}</span>{sub}</div><div class="c-status" role="cell">{pill}{detail}</div>'
                   '<div class="c-date" role="cell">{date}</div><div class="c-commit" role="cell">{commit}</div></div>'.format(
                       fail=" trow-fail" if row["status"] == "FAILED" else "", st=esc(row["status"]), name=esc(row["test"]),
                       sub='<span class="sub">{}</span>'.format(esc(sub)) if sub else "", pill=pill(row["status"]),
                       detail=detail, date=date_html, commit=commit))
    out.append('<div class="empty" data-empty hidden>No tests with this status.</div></div></section></div>')
    return "\n".join(out)


def _night_link(col: Any, up: str) -> str:
    """Link to the night's summary page in the archive, when there is one."""
    path = col.run["path"] if col.run else ""
    if path and col.www is not None and (Path(col.www) / path / "index.html").is_file():
        return '<span><a href="{}{}/index.html">Night summary and logs</a></span>'.format(up, esc(path))
    return ""


def log_href(www: Optional[Path], run_path: str, log: str) -> str:
    """The log as archived: plain, or compressed by `amps-nightly prune` (log.gz)."""
    if www is not None and not (Path(www) / run_path / log).exists() and (Path(www) / run_path / (log + ".gz")).exists():
        return log + ".gz"
    return log


def write_run_index(directory: Path) -> Path:
    """index.html of one archived night of one column: runners, failed tests and links to every log."""
    directory = Path(directory)
    run = json.loads((directory / "run.json").read_text())
    rows = [json.loads(l) for l in (directory / "results.jsonl").read_text().splitlines() if l.strip()]
    runners = json.loads((directory / "runners.json").read_text()) if (directory / "runners.json").is_file() else []

    def link(path: str, text: str) -> str:
        if not path:
            return ""
        href = path if (directory / path).exists() else path + ".gz"
        return '<a href="{}">{}</a>'.format(esc(href), esc(text)) if (directory / href).exists() else ""

    def count(status: str, source: str = "") -> int:
        return sum(1 for r in rows if r["status"] == status and (not source or r.get("source") == source))

    sha = run.get("commit", "")
    head = ['<div class="page tight">',
            '<div class="bar"><a class="back" href="../../../../../index.html">{}All test families</a></div>'.format(ICON["BACK"]),
            '<header class="fhead"><p class="eyebrow">Night of {} · {} · {}</p><h1>{}-{}</h1>'
            '<p class="src">Commit {} · started {} · finished {} · {} pass, {} failed, {} skipped{}</p></header>'.format(
                esc(night_of(run)), esc(run["machine"]), esc(run.get("compiler", "")),
                esc(run["machine"]), esc(run["config"]),
                commit_html(sha, bool(run.get("commit_inferred"))) if sha else "unknown",
                esc(run.get("started", "")), esc(run.get("finished", "")), count("PASS"), count("FAILED"), count("SKIPPED"),
                " · " + link("logs/orchestrator.log", "orchestrator log") if (directory / "logs").is_dir() else "")]
    if run.get("complete") is False and run.get("notes"):
        head.append('<p class="src"><span class="warn">{}Incomplete run: {}</span></p>'.format(ICON["WARN"], esc(run["notes"])))
    elif run.get("notes"):
        # companion commits (share, util, BATL, SWMF_data) and warnings of the ~/bin script (--note)
        parts = [p.strip() for p in run["notes"].split(";") if p.strip()]
        plain = [p for p in parts if not p.startswith("WARNING")]
        if plain:
            label = "With" if run.get("producer") == "amps-nightly run" else "Note:"
            head.append('<p class="src">{} {}</p>'.format(label, esc(", ".join(plain))))
        for warning in (p for p in parts if p.startswith("WARNING")):
            head.append('<p class="src"><span class="warn">{}{}</span></p>'.format(ICON["WARN"], esc(warning)))
    if runners:
        head.append('<section class="card night"><h2>Runners</h2><table class="tbl"><tr><th>Runner</th><th>State</th>'
                    '<th>Minutes</th><th>Pass</th><th>Failed</th><th>Skipped</th><th>Log</th><th>Note</th></tr>')
        for info in runners:
            minutes = "{:.0f}".format(info["elapsed_s"] / 60) if info.get("elapsed_s") is not None else ""
            note = info.get("problem") or info.get("note") or ("time limit reached" if info.get("timed_out") else "")
            head.append('<tr><td class="mono">{}</td><td>{}</td><td class="num">{}</td><td class="num">{}</td>'
                        '<td class="num">{}</td><td class="num">{}</td><td>{}</td><td>{}</td></tr>'.format(
                            esc(info["runner"]), esc(info.get("state", "")), minutes, count("PASS", info["runner"]),
                            count("FAILED", info["runner"]), count("SKIPPED", info["runner"]),
                            link(info.get("log", ""), "runner log"), esc(note)))
        head.append("</table></section>")
    failed = [r for r in rows if r["status"] == "FAILED"]
    head.append('<section class="card night"><h2>Failed tests ({})</h2>'.format(len(failed)))
    if failed:
        head.append('<table class="tbl"><tr><th>Family</th><th>Test</th><th>Stage</th><th>Reason</th><th>Log</th></tr>')
        for r in sorted(failed, key=lambda r: (r["family"], r["test"])):
            head.append('<tr><td>{}</td><td class="mono">{}</td><td>{}</td><td>{}</td><td>{}</td></tr>'.format(
                esc(r["family"]), esc(r["test"]), esc(r.get("stage", "")), esc(r.get("reason", "")),
                link(r.get("log", ""), "log")))
        head.append("</table>")
    head.append("</section></div>")
    title = "{}-{} night of {} · AMPS nightly tests".format(run["machine"], run["config"], night_of(run))
    (directory / "index.html").write_text(_page(title, "\n".join(head), 5, _when(dt.datetime.now().astimezone())))
    return directory / "index.html"


def col_label(col: Column) -> str:
    return col.label + (" · " + col.run["toolchain"] if col.run["toolchain"] else "")


# ----------------------------------------------------------------------------- build
def build(store: Store, registry: Registry, www: Path, now: Optional[dt.datetime] = None,
          tz: Optional[dt.tzinfo] = None, window_days: int = 60, stale_hours: float = 36,
          declared: Optional[List[Dict[str, str]]] = None) -> Dict[str, Any]:
    """Regenerate index.html, family pages and assets under www; runs/ is left untouched.

    declared: the columns of the enabled machine files (machines.declared_columns), in display order;
    each is shown, with "No results yet" until its first run. Stored columns that are not declared are
    hidden (their history stays in the store). None: show every column with a run in the last
    window_days days.
    """
    www = Path(www)
    now = (now or dt.datetime.now(dt.timezone.utc)).astimezone(dt.timezone.utc)
    now_local = now.astimezone(tz) if tz else now.astimezone()
    cutoff = now - dt.timedelta(days=window_days)
    stored = store.columns()
    hidden: List[str] = []
    if declared is None:
        all_columns: List[Any] = [Column(store, c["column_id"], now, tz, stale_hours) for c in stored
                                  if parse_time(c["last_finished"]) >= cutoff]
    else:
        have = {c["column_id"] for c in stored}
        wanted = [d["column"] for d in declared]
        all_columns = [Column(store, d["column"], now, tz, stale_hours) if d["column"] in have else Waiting(d)
                       for d in declared]
        hidden = [c["column_id"] for c in stored
                  if c["column_id"] not in wanted and parse_time(c["last_finished"]) >= cutoff]
    columns = [c for c in all_columns if c.run is not None]
    for col in all_columns:
        col.www = www
    present = set(f for c in columns for f in c.families())
    order = [f["id"] for f in registry.families]
    shown = [f for f in order if f in present] + sorted(present - set(order))
    unseen = [f for f in order if f not in present and registry.describe(f)["group"] != OTHER_GROUP]
    generated = _when(now_local)

    build_dir = www / ".build-{}".format(os.getpid())
    if build_dir.exists():
        shutil.rmtree(str(build_dir))
    (build_dir / "family").mkdir(parents=True)
    shutil.copytree(str(ASSETS), str(build_dir / "assets"))
    pages = 0
    if columns:
        body = _index(registry, all_columns, shown, unseen, now, now_local, stale_hours, hidden)
    else:
        lead = ("No results yet from: {}.".format(", ".join("{} · {}".format(c.machine, c.label) for c in all_columns))
                if all_columns else "No runs finished in the last {} days.".format(window_days))
        body = ('<div class="page"><header class="top"><div><p class="eyebrow">AMPS · SWMFsoftware/AMPS</p>'
                '<h1>Nightly test results</h1><p class="lead">{}</p></div></header></div>'.format(esc(lead)))
    (build_dir / "index.html").write_text(_page("AMPS nightly test results", body, 0, generated))
    for fid in shown:
        family_dir = build_dir / "family" / fid
        family_dir.mkdir()
        for col in columns:
            title = "{} · {} {} · AMPS nightly tests".format(registry.describe(fid)["title"], col.machine, col.label)
            (family_dir / "{}.html".format(col.id)).write_text(
                _page(title, _family_page(registry, fid, col, columns), 2, generated))
            pages += 1

    www.mkdir(parents=True, exist_ok=True)
    for name in ("assets", "family"):
        old = www / ".old-{}-{}".format(name, os.getpid())
        if (www / name).exists():
            os.replace(str(www / name), str(old))
        os.replace(str(build_dir / name), str(www / name))
        if old.exists():
            shutil.rmtree(str(old))
    os.replace(str(build_dir / "index.html"), str(www / "index.html"))
    shutil.rmtree(str(build_dir))
    return {"columns": [c.id for c in columns], "waiting": [c.id for c in all_columns if c.run is None],
            "hidden": hidden, "families": shown, "family_pages": pages}
