#!/bin/bash
# Check that every version-controlled copy listed in RUN_COPIES.sha256 equals its run copy once the header block is
# removed (the header lines between the two "==== ... ====" marker lines are the only addition). Marker forms:
#     "# ==== version-controlled copy: header added ..."     ... "# ==== end of header ===="       (Python, shell)
#     "// ==== version-controlled copy: header added ..."    ... "// ==== end of header ===="      (kernels/*.cu)
#     "<!-- ==== version-controlled copy: header added ..." ... "==== end of header ==== -->"     (Markdown)
# A file without a header (grain_*.py/.sbatch, committed byte-identical) is checked as it is.
#     model_ext/check_copies.sh            (from anywhere; prints OK/DIFFERS per file, exit 1 if any differs)
cd "$(dirname "$0")" || exit 1
bad=0
while read -r sha here run; do
    case $sha in \#*|"") continue;; esac
    got=$(sed -e '/^# ==== version-controlled copy: header added/,/^# ==== end of header ====/d' \
              -e '\%^// ==== version-controlled copy: header added%,\%^// ==== end of header ====%d' \
              -e '/^<!-- ==== version-controlled copy: header added/,/^==== end of header ==== -->/d' \
              "$here" | sha256sum | cut -d' ' -f1)
    if [ "$got" = "$sha" ]; then echo "OK      $here"; else echo "DIFFERS $here"; bad=1; fi
done < RUN_COPIES.sha256
exit $bad
