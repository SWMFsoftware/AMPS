# AMPS user and developer manual

**AMPS.tex** is the main file. The title page uses the supplied, unmodified
AMPS_clean_sweep_white.jpg. This edition covers core/application structure,
installation and execution, the two input layers, kinetic numerics, AMR,
movers and fields, collisions/internal energy, photochemistry, dust, surfaces
and arbitrary shapes, planetary/comet and heliospheric models, SWMF coupling,
sampling/observations/restart, testing, end-to-end examples, development and
troubleshooting. The appendix includes every registered test-suite owner.

The PDF records implementation and qualification boundaries explicitly.
It does not claim that every application or combination of physics is qualified.

## Build

From the AMPS root:

~~~sh
make PDF
~~~

The existing root PDF target delegates to this directory. The equivalent
documentation-only commands, which also work without AMPS configuration, are:

~~~sh
make -C doc PDF
make -C doc MANUAL
~~~

Inside a standalone copy of doc, use make PDF. Output is **AMPS.pdf in
that directory**. Compiler transcripts, dependency records and
build-receipt.json are retained in build-manual/.
None of these commands builds, configures or runs the AMPS physics code.

Requirements:

- Python 3.8 or newer and make.
- Either a working pdfLaTeX distribution or Tectonic.
- Ordinary TeX packages: fontenc/inputenc, Latin Modern, geometry, graphicx,
  amsmath/amssymb, booktabs, longtable, tabularx/array, xcolor, listings,
  microtype, hyperref, bookmark and xurl. They are standard distribution
  packages; no SWMF style, custom bibliography style or shell escape is used.

The build driver prefers pdflatex, then tectonic in PATH, then a local
.build-tools/tectonic if present. An explicit executable can be selected:

~~~sh
make -C doc PDF LATEX_ENGINE=/path/to/pdflatex
make -C doc PDF LATEX_ENGINE=/path/to/tectonic
~~~

Tectonic's first use retrieves ordinary TeX resources. Its caches/configuration
stay under this directory in .tex-cache/ and .tex-config/; subsequent
fully cached builds can use:

~~~sh
cd doc
python3 -B build_manual.py --engine tectonic --only-cached
~~~

The optional local compiler used for this edition was downloaded from the
[official Tectonic release](https://github.com/tectonic-typesetting/tectonic/releases/tag/tectonic%400.17.0);
its URL and SHA-256 are in .build-tools/download-receipt.json.
That generated compiler/cache is not a required source asset. A recipient can
use their own TeX installation. The maintained
[Tectonic installation instructions](https://github.com/tectonic-typesetting/tectonic/blob/master/docs/src/installation/index.md)
describe platform-specific distributions.

## Self-contained source files

Chapter source filenames omit the AMPS_ prefix. AMPS.tex remains the main
file, and the supplied logo keeps its original filename. All actual input
files and figures used by AMPS.tex live here:

- AMPS.tex, HEADER.tex, the chapter .tex files and suite appendix.
- AMPS_clean_sweep_white.jpg.
- figures/kinetic_reference.pdf, a supplied analytical teaching figure with
  explicitly labelled assumptions. Its PNG and JSON companion retain the
  visual/numeric reference; make_reference_figures.py regenerates them when
  NumPy/matplotlib are installed. It is not simulation output.
- Embedded bibliography References.tex; no BibTeX pass or external
  agu.bst is needed. The original AMPS.bib is retained for provenance.
- Makefile, build_manual.py, check_manual.py and this README.

Repository paths printed in prose and example commands are reader references,
not build dependencies. Copying this directory's maintained manual sources
and logo to another location does not require copying src, test, SWMF,
SPICE, production runs or observation archives to compile the PDF.

The optional update_test_catalog.py editorial tool needs the surrounding
AMPS registry when refreshing the suite appendix. Normal PDF builds consume
the local appendix and do not run that generator.

## Validation and maintenance

~~~sh
make -C doc check
python3 -B doc/update_test_catalog.py --check
~~~

The first checks the actual include graph, local graphics, duplicate/missing
labels and citations using only the standard library. The second verifies the
suite snapshot against the expanded registry in a complete AMPS checkout.
After changing registrations:

~~~sh
python3 -B doc/update_test_catalog.py
~~~

The compile driver performs three pdfLaTeX passes, or Tectonic's automatic
reruns. It rejects unresolved references/citations and publishes AMPS.pdf
only after success. Inspect layout warnings and representative rendered pages;
a structurally valid TeX source can still have poor tables, clipping or
misleading scientific prose. Build validation is documentation evidence,
not a new execution of any printed test command.

Update STATUS.md with exact validation results. Maintain source pointers and
explicit units, assumptions, valid domains, independent references and
unsupported behavior. Keep command excerpts clearly distinct from complete
runnable decks. Do not copy a current run output into the maintained source
reference set.

## Generated files and cleanup

The clean target removes only build-manual/. The distclean target also removes
the generated AMPS.pdf. Both preserve all TeX sources, the logo, the historical
notes, local compiler and package caches. The driver refuses a symlink as its
build directory.

Optional compiler/cache files, build logs, rendered previews and PDFs are
generated artifacts. They do not need to accompany the maintained source
bundle when a recipient has a suitable TeX distribution.

## Historical material

The original idf.tex, ecsim.tex, focused-transport, Parker, shock and
turbulence notes remain intact as supplementary research material. They are
not automatically included in the new edition because their historical
conventions and implementation scope require separate reconciliation.
Four replaced operator chapters are preserved under legacy/; those CVS,
old platform and old syntax instructions are not current operation guidance.
