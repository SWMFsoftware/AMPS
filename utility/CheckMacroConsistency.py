#!/usr/bin/env python3
"""Check that every preprocessor macro has one value in all translation units of an AMPS build.

AMPS selects its configuration with macros.  ampsConfig.pl and the application's configuration
scripts write the application input's settings into the build copy of the headers (picGlobal.dfn,
Exosphere.dfn, meshAMRdef.h, ...).  If such a setting lands outside a header's include guard, every
repeated #include of the header re-applies it in the middle of a translation unit (TU) and undoes
whatever was set in between, so different TUs of the same executable can see different values of the
same macro.  This script detects that.

Usage, from the AMPS directory of a configured application (after ./Config.pl -application=... and
'make build' or a complete 'make amps'):

    utility/CheckMacroConsistency.py [--jobs N] [--only REGEX] [--ignore REGEX] [--all] [--strict]
                                     [--json FILE] [--verbose]

1. 'make -n -B -w amps_after_build' (no file is changed) gives the real compile command of every
   C/C++ TU of the application.
2. Every TU is preprocessed with its own flags plus '-dM -E'; the final value of each macro is
   collected.  Every macro that ends with different values in different TUs is reported: 'effective'
   when TUs with different values refer to the macro, 'latent' when only the TUs with one of the values
   refer to it (a file included by the TU names the macro other than in its #define/#undef/#ifdef).
3. Every TU is also preprocessed with '-dD -E'.  This gives the file:line of each final definition and
   finds the values that a repeated #include of a header changes in the middle of a TU (definitions
   outside the header's include guard).  Such a TU compiles the code before and after the repeated
   #include with different values, even when all TUs end with the same value.
Macros whose definitions all come from system headers or the compiler are ignored unless --all.
Limitation: include paths that a makefile passes through the environment (srcSEP3D/makefile sets
CPLUS_INCLUDE_PATH per target) do not appear in the dry run, so such TUs fail to preprocess (exit 2).

Exit status: 0 = one value per macro, 1 = effective differences or values changed by a repeated
#include (with --strict also latent differences), 2 = error.
"""

import argparse
import collections
import concurrent.futures
import json
import os
import re
import shlex
import signal
import subprocess
import sys

SOURCE_EXTENSIONS = ('.cpp', '.cc', '.cxx', '.C', '.c')
# Macros whose value depends on the TU by design and that are not configuration settings:
# SRCSEP_SEP_COMMON_HEADER (srcSEP/util/sep_common_header_path.h) selects the include root under which the
# TU finds the same src/models/sep_common header.
DEFAULT_IGNORE = r'^SRCSEP_SEP_COMMON_HEADER$'
ENTER_RE = re.compile(r"^\S*make(?:\[\d+\])?: Entering directory [`'](.*)'\s*$")
LEAVE_RE = re.compile(r"^\S*make(?:\[\d+\])?: Leaving directory [`'](.*)'\s*$")
MARKER_RE = re.compile(r'^# (\d+) "((?:[^"\\]|\\.)*)"((?: \d+)*)\s*$')
DEFINE_RE = re.compile(r'^#\s*define\s+([A-Za-z_]\w*)(.*)$')
UNDEF_RE = re.compile(r'^#\s*undef\s+([A-Za-z_]\w*)')


def die(message):
    print(message, file=sys.stderr)
    sys.exit(2)


def norm_value(rest):
    """Value of '#define NAME<rest>': the parameter list (if any) and the body, blanks collapsed."""
    if rest.startswith('('):
        close = rest.find(')')
        params = re.sub(r'\s+', '', rest[:close + 1])
        return params + ' ' + ' '.join(rest[close + 1:].split())
    return ' '.join(rest.split())


def split_shell(line):
    """Split a shell line into simple commands (lists of words) at ';', '&&' and '||'."""
    lexer = shlex.shlex(line, posix=True, punctuation_chars=True)
    lexer.whitespace_split = True
    commands, current = [], []
    for token in lexer:
        if token in (';', '&&', '||', '&', '|'):
            if current:
                commands.append(current)
            current = []
        else:
            current.append(token)
    if current:
        commands.append(current)
    return commands


def compile_commands(amps_dir, target, make):
    """Run make dry-run and return [(cwd, source, argv)] for every C/C++ compile command."""
    # The recipes run sub-makes as 'cd build/X; $(MAKE) ...'.  When build/X is missing, the sub-make runs in
    # the AMPS directory again and recurses without end: the dry run is stopped (as a process group) as soon
    # as a sub-make of level 2 or deeper enters the AMPS directory.
    proc = subprocess.Popen([make, '-n', '-B', '-w', target], cwd=amps_dir, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, universal_newlines=True, start_new_session=True)
    lines = []
    for line in proc.stdout:
        lines.append(line)
        m = re.match(r"^\S*make\[(\d+)\]: Entering directory [`'](.*)'", line)
        if m and int(m.group(1)) >= 2 and os.path.realpath(m.group(2)) == amps_dir:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
            die('CheckMacroConsistency: "%s -n -B -w %s" recurses into %s (a build/ subdirectory is missing?)'
                     % (make, target, amps_dir))
    proc.wait()
    output = ''.join(lines)
    if proc.returncode != 0:
        die('CheckMacroConsistency: "%s -n -B -w %s" failed:\n%s' % (make, target, output[-3000:]))
    dirs, result, seen = [amps_dir], [], set()
    for line in output.splitlines():
        m = ENTER_RE.match(line)
        if m:
            dirs.append(m.group(1))
            continue
        if LEAVE_RE.match(line):
            if len(dirs) > 1:
                dirs.pop()
            continue
        if ' -c ' not in line:
            continue
        try:
            commands = split_shell(line)
        except ValueError:
            continue
        cwd = dirs[-1]
        for argv in commands:
            if argv[0] == 'cd' and len(argv) > 1:
                cwd = os.path.normpath(os.path.join(cwd, argv[1]))
                continue
            if '-c' not in argv or 'fortran' in os.path.basename(argv[0]):
                continue
            sources = [a for a in argv[1:] if a.endswith(SOURCE_EXTENSIONS) and not a.startswith('-')]
            if len(sources) != 1:
                continue
            path = os.path.normpath(os.path.join(cwd, sources[0]))
            if path not in seen:
                seen.add(path)
                result.append((cwd, path, argv))
    return result


def preprocess_argv(argv, mode):
    """The compile command turned into a preprocessor run that writes to stdout."""
    out, skip = [], False
    for a in argv:
        if skip:
            skip = False
            continue
        if a in ('-o', '-MF', '-MT', '-MQ'):
            skip = True
            continue
        if a in ('-c', '-MD', '-MMD', '-MP', '-M', '-MM') or (a.startswith('-o') and len(a) > 2):
            continue
        out.append(a)
    return out + [mode, '-E']


def run(argv, cwd):
    proc = subprocess.run(argv, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True, errors='replace')
    if proc.returncode != 0:
        raise RuntimeError('preprocessing failed in %s:\n  %s\n%s' % (cwd, ' '.join(argv), proc.stderr[-2000:]))
    return proc.stdout


def final_values(text):
    values = {}
    for line in text.splitlines():
        m = DEFINE_RE.match(line)
        if m:
            values[m.group(1)] = norm_value(m.group(2))
    return values


def trace(text, cwd):
    """Replay the #define/#undef lines of '-dD -E' output.

    Returns (final, origin, reapplied, files): final[name] = value at the end of the TU, origin[name] =
    'file:line' of that definition, reapplied = net changes of macro values made by a second or later
    inclusion of the same file (definitions outside the file's include guard), files = files included."""
    state, origin, reapplied = {}, {}, []
    entries = collections.Counter()
    stack = []            # [file, inclusion number, 'file:line' it was included from, {macro: (old, old_from)}]
    lineno = 0

    def close(frame):
        for name, (old, old_from) in frame[3].items():
            if old is not None and state.get(name) != old:
                reapplied.append({'macro': name, 'at': origin.get(name, frame[0]), 'inclusion': frame[1],
                                  'included_from': frame[2], 'old': old, 'old_from': old_from,
                                  'new': state.get(name)})

    for line in text.split('\n'):
        if not line.startswith('#'):
            lineno += 1
            continue
        m = MARKER_RE.match(line)
        if m:
            name = m.group(2)
            if not name.startswith('<'):
                name = os.path.normpath(os.path.join(cwd, name))
            flags = m.group(3).split()
            if '1' in flags:
                parent = '%s:%d' % (stack[-1][0], lineno) if stack else ''
                entries[name] += 1
                stack.append([name, entries[name], parent, {}])
            elif '2' in flags:
                while stack and stack[-1][0] != name:
                    close(stack.pop())
                if not stack:
                    stack.append([name, 1, '', {}])
            elif not stack or stack[-1][0] != name:
                if len(stack) == 1:
                    close(stack.pop())
                stack.append([name, 1, '', {}])
            lineno = int(m.group(1))
            continue
        md = DEFINE_RE.match(line)
        if md:
            name, value = md.group(1), norm_value(md.group(2))
        else:
            mu = UNDEF_RE.match(line)
            name, value = (mu.group(1), None) if mu else (None, None)
        if name is not None:
            if stack and stack[-1][1] > 1 and name not in stack[-1][3]:
                stack[-1][3][name] = (state.get(name), origin.get(name))
            if value is None:
                state.pop(name, None)
                origin.pop(name, None)
            else:
                state[name] = value
                origin[name] = '%s:%d' % (stack[-1][0], lineno) if stack else '?:%d' % lineno
        lineno += 1
    while stack:
        close(stack.pop())
    return state, origin, reapplied, set(entries)


_REFERENCES = {}
IDENT_RE = re.compile(r'[A-Za-z_]\w*')
COMMENT_RE = re.compile(r'/\*.*?\*/|//[^\n]*', re.S)
DIRECTIVE_NAME_RE = re.compile(r'^\s*#\s*(?:define|undef|ifdef|ifndef)\s+([A-Za-z_]\w*)')


def references(path):
    """Identifiers that a source file refers to, other than the name of a #define/#undef/#ifdef/#ifndef."""
    if path not in _REFERENCES:
        names = set()
        try:
            with open(path, errors='replace') as f:
                text = COMMENT_RE.sub(' ', f.read())
        except OSError:
            text = ''
        for line in text.splitlines():
            m = DIRECTIVE_NAME_RE.match(line)
            if m:
                line = line[m.end():]
            names.update(IDENT_RE.findall(line))
        _REFERENCES[path] = names
    return _REFERENCES[path]


def check_tu(item, do_trace):
    cwd, source, argv = item
    result = {'source': source}
    result['final'] = final_values(run(preprocess_argv(argv, '-dM'), cwd))
    if do_trace:
        final, origin, reapplied, files = trace(run(preprocess_argv(argv, '-dD'), cwd), cwd)
        mismatch = [n for n in set(final) | set(result['final'])
                    if final.get(n) != result['final'].get(n) and not n.startswith('__')]
        result.update(origin=origin, reapplied=reapplied, files=files | {source}, trace_mismatch=sorted(mismatch))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--amps-dir', default='.', help='AMPS directory (default: current directory)')
    parser.add_argument('--target', default='amps_after_build', help='make target whose compile commands are checked')
    parser.add_argument('--make', default='make')
    parser.add_argument('--jobs', '-j', type=int, default=8)
    parser.add_argument('--no-trace', action='store_true', help='skip the -dD pass (no origins, no re-applied definitions)')
    parser.add_argument('--all', action='store_true', help='also report macros defined only by system headers or the compiler')
    parser.add_argument('--strict', action='store_true', help='exit 1 also for latent differences')
    parser.add_argument('--only', metavar='REGEX', help='check only the TUs whose source path matches REGEX')
    parser.add_argument('--ignore', metavar='REGEX', default=DEFAULT_IGNORE,
                        help='macros not to compare (default: %(default)s, TU-dependent by design)')
    parser.add_argument('--json', help='write the full result to this file')
    parser.add_argument('--verbose', '-v', action='store_true', help='list all TUs for every value')
    args = parser.parse_args()

    amps_dir = os.path.realpath(args.amps_dir)
    build_dir = os.path.join(amps_dir, 'build')
    missing = [d for d in ('general', 'meshAMR', 'pic', 'species', 'main') if not os.path.isdir(os.path.join(build_dir, d))]
    if missing:
        die('CheckMacroConsistency: %s/{%s} missing; configure the application and run "make build" first'
                 % (build_dir, ','.join(missing)))
    rel = lambda p: os.path.relpath(p, amps_dir) if p.startswith(amps_dir + os.sep) else p

    items = compile_commands(amps_dir, args.target, args.make)
    if not items:
        die('CheckMacroConsistency: no C/C++ compile commands found for target %s' % args.target)
    compiled = set(source for _, source, _ in items)
    total = len(items)
    if args.only:
        items = [item for item in items if re.search(args.only, item[1])]
    not_compiled = sorted(rel(os.path.join(d, f)) for d, _, files in os.walk(build_dir) for f in files
                          if f.endswith(SOURCE_EXTENSIONS) and os.path.join(d, f) not in compiled)

    results, errors = [], []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = {pool.submit(check_tu, item, not args.no_trace): item for item in items}
        for future in concurrent.futures.as_completed(futures):
            try:
                results.append(future.result())
            except RuntimeError as error:
                errors.append(str(error))
    results.sort(key=lambda r: r['source'])

    def in_tree(location):
        return location is not None and (location.startswith(amps_dir + os.sep) or location.startswith('<command-line>'))

    # macro -> value -> [TU]
    values = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in results:
        for name, value in r['final'].items():
            values[name][value].append(r)
    differing = []
    for name in sorted(values):
        if len(values[name]) < 2 or (args.ignore and re.search(args.ignore, name)):
            continue
        origins = [r['origin'].get(name) for v in values[name] for r in values[name][v] if 'origin' in r]
        if not args.all:
            if args.no_trace:
                if name.startswith('__') or not re.match(r'_?[A-Z]', name):
                    continue
            elif not any(in_tree(o) for o in origins):
                continue
        differing.append(name)

    reapplied = collections.OrderedDict()
    for r in results:
        for event in r.get('reapplied', []):
            if not args.all and not in_tree(event['at']):
                continue
            key = (event['macro'], event['at'], event['old'], event['new'])
            tus = reapplied.setdefault(key, {'event': event, 'tus': []})['tus']
            if r['source'] not in tus:
                tus.append(r['source'])

    def referenced(name, tu):
        return any(name in references(f) for f in tu.get('files', ()) if f.startswith(amps_dir + os.sep))

    print('CheckMacroConsistency: %s' % amps_dir)
    app = ''
    try:
        with open(os.path.join(amps_dir, 'Makefile.local')) as f:
            settings = collections.OrderedDict(l.strip().split('=', 1) for l in f
                                               if l.startswith(('APPLICATION=', 'InputFileAMPS=', 'TESTMODE=')))
            app = ' '.join('%s=%s' % kv for kv in settings.items())
    except OSError:
        pass
    print('  configuration: %s' % app)
    print('  translation units: %d compiled by "make %s"%s; %d other C/C++ sources in build/ are not compiled'
          % (total, args.target, ', %d of them checked (--only)' % len(results) if args.only else '', len(not_compiled)))
    for error in errors:
        print('  ERROR: ' + error)
    mismatches = sorted(set(n for r in results for n in r.get('trace_mismatch', [])))
    if mismatches:
        print('  WARNING: -dM and -dD disagree for: ' + ' '.join(mismatches))

    # A difference is 'effective' when TUs with different values refer to the macro (code is compiled
    # with different settings), 'latent' when only the TUs with one of the values refer to it.
    effective, latent = [], []
    for name in differing:
        if args.no_trace:
            effective.append(name)
            continue
        groups = sum(1 for tus in values[name].values() if any(referenced(name, t) for t in tus))
        (effective if groups > 1 else latent).append(name)

    def show(tus):
        shown = [rel(t if isinstance(t, str) else t['source']) for t in tus]
        if not args.verbose and len(shown) > 4:
            shown = shown[:4] + ['... (%d more)' % (len(tus) - 4)]
        return ', '.join(shown)

    for title, names in (('Macros with different values in TUs that refer to them', effective),
                         ('Macros with different values, referred to only in TUs with one of the values (latent)', latent)):
        print('\n  %s: %d' % (title, len(names)))
        for name in names:
            print('    %s' % name)
            for value, tus in sorted(values[name].items(), key=lambda kv: -len(kv[1])):
                where = sorted(set(rel(t['origin'].get(name, '?')) for t in tus if 'origin' in t))
                print('      = %s' % value)
                if args.no_trace:
                    print('          in %d TUs: %s' % (len(tus), show(tus)))
                else:
                    used = sum(1 for t in tus if referenced(name, t))
                    print('          in %d TUs (%d of them refer to it): %s' % (len(tus), used, show(tus)))
                if where:
                    print('          last defined at: %s' % ', '.join(where))
            undefined = len(results) - sum(len(t) for t in values[name].values())
            if undefined:
                print('      (undefined in %d TUs)' % undefined)

    if not args.no_trace:
        print('\n  Macro values changed in the middle of a TU by a repeated #include of a header'
              ' (definitions outside its include guard): %d' % len(reapplied))
        for (name, at, old, new), entry in reapplied.items():
            event, tus = entry['event'], entry['tus']
            print('    %s: %s (inclusion #%d of the file, from %s)'
                  % (name, rel(at), event['inclusion'], rel(event['included_from'] or '?')))
            print('        resets %s (set at %s)' % (old, rel(event['old_from'] or '?')))
            print('        to     %s' % new)
            print('        in %d TUs: %s' % (len(tus), show(tus)))

    if args.json:
        with open(args.json, 'w') as f:
            json.dump({'amps_dir': amps_dir, 'configuration': app, 'target': args.target,
                       'translation_units': [rel(r['source']) for r in results],
                       'not_compiled': not_compiled, 'errors': errors,
                       'differing': {n: {v: [rel(t['source']) for t in tus] for v, tus in values[n].items()}
                                     for n in differing},
                       'effective': effective, 'latent': latent,
                       'reapplied': [dict(e['event'], tus=[rel(t) for t in e['tus']]) for e in reapplied.values()]},
                      f, indent=1)

    if errors:
        return 2
    if effective or reapplied or (args.strict and latent):
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
