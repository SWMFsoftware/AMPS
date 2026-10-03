#!/usr/bin/env python3
"""Create a filtered source archive using the Python 3.8 standard library."""

import argparse
import os
from pathlib import Path
import re
import stat
import sys
import tarfile
import tempfile


DEFAULT_EXTENSIONS = {".md", ".cpp", ".h", ".hpp", ".cc", ".f", ".f90", ".for"}

HELP = """Selection rules:
  Directories are searched recursively, including hidden directories.
  Include .md, .cpp, .h, .hpp, .cc, .f, .f90, and .for files by default.
  Include Makefile, makefile, GNUmakefile, Makefile.* and GNUmakefile.*,
  plus files ending in .mk or .mak. Names/extensions are case-insensitive.
  Additional extensions supplement (never replace) the default selection.
  Explicitly listed files also must satisfy these selection rules.
  Overlapping inputs are deduplicated. Symbolic links and special files
  are skipped with warnings; directory links are never followed.
  Parent directories of selected files are added once as directory entries.
  Empty directories are omitted. The archive root itself is not included.
  Statistics report files, directories, total source bytes and archive bytes.
  Dry-run reports selected file/directory counts without writing an archive.

Additional extensions (-e):
  Supply one comma-separated or quoted space-separated list per -e option.
  Repeat -e as needed. Leading dots are optional; wildcards are not needed.
  Examples: -e py       -e .py,.ini       -e 'py sh txt'       -e py -e ini
  Compound extensions such as .in.template are supported.
  For several unquoted extensions, use commas or repeat -e; -e consumes
  exactly one argument, leaving subsequent arguments as input paths.

Archive paths:
  Normally paths are relative to the current working directory.
  Example: src/models/foo.cpp remains src/models/foo.cpp in the archive.
  If any input is outside that directory, use the common ancestor of the
  current directory, input directories and input file parents as the root.
  The root is printed before writing. No absolute paths or '..' are stored.
  Only the listed files/directories are scanned, regardless of root choice.

Compression and existing files:
  .tar.gz / .tgz       gzip-compressed tar
  .tar.bz2 / .tbz2     bzip2-compressed tar
  .tar.xz / .txz       xz-compressed tar
  Any other name      uncompressed tar
  Existing output is refused unless --force is supplied. The archive is
  written to a temporary file and published only after successful writing.
  The output archive itself is always excluded from the selected inputs.
  Missing/unreadable inputs or no matching files cause a nonzero exit.

Examples:
  # Archive two source trees and their Makefile.
  tar_sources -f sources.tar srcSEP srcSEP3D Makefile

  # Include Python scripts, INI inputs and shell scripts as well.
  tar_sources -f sources.tar.gz -e py,ini,sh src/models srcSEP3D

  # Repeated -e options; uppercase Fortran .F90 is included automatically.
  tar_sources -f sources.tgz -e .py -e .txt src test

  # Preview selected archive members without creating an archive.
  tar_sources -f sources.tar -e 'py ini' --dry-run .

  # Paths containing spaces must be quoted.
  tar_sources -f sources.tar 'project one/src' 'project two/Makefile'

  # Replace an existing archive and print each included path.
  tar_sources -f sources.tar.gz --force --verbose .

  # Use -- before paths beginning with a dash.
  tar_sources -f sources.tar -- -project

  # Run without installation or execute permission.
  python3.8 ./tar_sources -f sources.tar -e py src

Installation (Python 3.8 or later; no third-party packages):
  # Direct execution uses python3.8 from PATH via the shebang.
  # With a newer interpreter, run: python3 ./tar_sources ...
  chmod +x tar_sources
  mkdir -p "$HOME/.local/bin"
  cp tar_sources "$HOME/.local/bin/"
  # Ensure $HOME/.local/bin is on your PATH, or use ./tar_sources directly.

Exit status: 0 = success; 1 = filesystem/archive error; 2 = invalid CLI.
"""


def build_parser():
    parser = argparse.ArgumentParser(
        prog="tar_sources",
        description="Archive selected source files while preserving their directory tree.",
        epilog=HELP,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=False,
    )
    parser.add_argument("-f", "--file", required=True, metavar="TAR_FILE",
                        help="output tar archive; compression inferred from its name")
    parser.add_argument("-e", "--extensions", action="append", default=[], metavar="EXTENSIONS",
                        help="additional comma/space-separated extensions; repeatable")
    parser.add_argument("-v", "--verbose", action="store_true", help="print each archived path")
    parser.add_argument("--dry-run", action="store_true", help="list selected paths without writing")
    parser.add_argument("--force", action="store_true", help="replace an existing output archive")
    parser.add_argument("paths", nargs="+", metavar="PATH", help="files and/or directories to archive")
    return parser


def parse_extensions(values, parser):
    extensions = set(DEFAULT_EXTENSIONS)
    for value in values:
        tokens = re.split(r"[,\s]+", value.strip())
        if not any(tokens):
            parser.error("-e requires at least one extension")
        for token in filter(None, tokens):
            extension = token.lower().lstrip(".")
            if not extension or any(c in extension for c in "/\\*?[]") or extension.endswith("."):
                parser.error("invalid extension {!r}; use e.g. 'py,ini'".format(token))
            extensions.add("." + extension)
    return extensions


def selected_name(name, extensions):
    name = name.lower()
    makefile = (name in {"makefile", "gnumakefile"}
                or name.startswith(("makefile.", "gnumakefile."))
                or name.endswith((".mk", ".mak")))
    return makefile or name.endswith(tuple(extensions))


def absolute_path(value):
    # Normalize '..' without resolving symlinks: links are skipped, not followed.
    return Path(os.path.abspath(os.path.expanduser(value)))


def collect_files(inputs, output, extensions):
    """Walk only requested trees; propagate access errors rather than silently omit files."""
    files = set()

    def consider(path):
        mode = path.lstat().st_mode
        if stat.S_ISLNK(mode):
            print("warning: skipping symbolic link: {}".format(path), file=sys.stderr)
        elif stat.S_ISREG(mode):
            # samefile also detects a hard link to an existing output archive.
            is_output = path == output or (output.exists() and os.path.samefile(path, output))
            if not is_output and selected_name(path.name, extensions):
                files.add(path)
        elif not stat.S_ISDIR(mode):
            print("warning: skipping special file: {}".format(path), file=sys.stderr)

    def walk_error(error):
        raise error

    for path in inputs:
        mode = path.lstat().st_mode  # Also rejects missing explicitly listed inputs.
        if not stat.S_ISDIR(mode):
            consider(path)
            continue
        for directory, directories, filenames in os.walk(path, onerror=walk_error, followlinks=False):
            base = Path(directory)
            # os.walk puts directory symlinks in directories even with followlinks=False.
            for name in sorted(directories):
                if (base / name).is_symlink():
                    consider(base / name)
                    directories.remove(name)
            directories.sort()
            for name in sorted(filenames):
                consider(base / name)
    return files


def compression_mode(output):
    name = output.name.lower()
    for suffixes, mode in [((".tar.gz", ".tgz"), "w:gz"),
                           ((".tar.bz2", ".tbz2"), "w:bz2"),
                           ((".tar.xz", ".txz"), "w:xz")]:
        if name.endswith(suffixes):
            return mode
    return "w"


def directory_members(root, members):
    """Include each ancestor once, excluding the archive root and empty trees."""
    names = set()
    for _, name in members:
        for parent in Path(name).parents:
            if parent != Path("."):
                names.add(parent.as_posix())
    # Parents precede their children, independent of filesystem iteration order.
    return [(root / name, name) for name in sorted(names, key=lambda n: (n.count("/"), n))]


def write_archive(output, members, directories, force, verbose):
    """Publish atomically; a failed read never leaves a partial final archive."""
    descriptor, temporary = tempfile.mkstemp(prefix=".tar_sources-", dir=str(output.parent))
    os.close(descriptor)
    source_bytes = 0
    try:
        with tarfile.open(temporary, compression_mode(output), dereference=False) as archive:
            for path, name in directories:
                info = archive.gettarinfo(str(path), arcname=name)
                if not info.isdir():
                    raise OSError("archive parent is no longer a directory: {}".format(path))
                archive.addfile(info)
                if verbose:
                    print(name + "/")
            for path, name in members:
                # Reject a file that changed into a link/special file after selection.
                info = archive.gettarinfo(str(path), arcname=name)
                if not info.isfile():
                    raise OSError("selected input is no longer a regular file: {}".format(path))
                with path.open("rb") as source:
                    archive.addfile(info, source)
                source_bytes += info.size
                if verbose:
                    print(name)
        if force:
            os.replace(temporary, output)
        else:
            # Atomic exclusive publication: do not overwrite even if another process
            # created the requested output while this archive was being written.
            os.link(temporary, output)
        return source_bytes
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    extensions = parse_extensions(args.extensions, parser)
    output = absolute_path(args.file)
    inputs = [absolute_path(value) for value in args.paths]
    try:
        if not args.dry_run and os.path.lexists(output) and not args.force:
            raise OSError("output already exists (use --force to replace): {}".format(output))
        if output.is_symlink() or output.is_dir():
            raise OSError("output must not be a symbolic link or directory: {}".format(output))
        files = collect_files(inputs, output, extensions)
        if not files:
            raise OSError("no matching source files found; no archive created")
        # Keep the working-directory prefix when possible. For external inputs,
        # broaden the root without creating absolute or traversal member names.
        root = Path(os.path.commonpath(
            [str(Path.cwd())] + [str(p if p.is_dir() else p.parent) for p in inputs]
        ))
        members = sorted(((path, path.relative_to(root).as_posix()) for path in files),
                         key=lambda item: item[1])
        directories = directory_members(root, members)
        print("Archive root: {}".format(root))
        if args.dry_run:
            for _, name in directories:
                print(name + "/")
            for _, name in members:
                print(name)
            print("Selected {} file(s), {} directory/directories; no archive written."
                  .format(len(members), len(directories)))
        else:
            source_bytes = write_archive(output, members, directories, args.force, args.verbose)
            print("Created {}".format(output))
            print("Files added:       {}".format(len(members)))
            print("Directories added: {}".format(len(directories)))
            print("Source bytes:      {:,}".format(source_bytes))
            print("Archive bytes:     {:,}".format(output.stat().st_size))
        return 0
    except (OSError, ValueError, tarfile.TarError) as error:
        print("tar_sources: error: {}".format(error), file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("tar_sources: interrupted", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
