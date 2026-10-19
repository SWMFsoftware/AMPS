# Third-party material

## polyhedral-gravity-model 3.3 (GNU GPL v3)

| Path | Content | SHA-256 |
|---|---|---|
| `download_polyhedral_gravity_v3.3/polyhedral_gravity-3.3.tar.gz` | PyPI sdist, the canonical source for the C03 translation | `1476ab9d638538f70b4d31e7cb1efce032c56c95b04264d88d005929769e7715` |
| `download_polyhedral_gravity_v3.3/polyhedral_gravity-3.3-cp38-cp38-manylinux_2_17_x86_64.manylinux2014_x86_64.whl` | PyPI wheel; its `.so` is byte-identical to the campaign binary in `/tmp/liouville_cops_polyhedral_gravity` | `a3b7fbe47ff82d9caa46944744aaf86d9feb11b35832a86afa09e27c857be671` |
| `src_pypi/polyhedral_gravity-3.3/` | the sdist, unpacked as is (read only, never built here) | — |

Both digests match the values PyPI publishes for the release (checked 2026-10-19 with
`sha256sum -c`). The license is in `src_pypi/polyhedral_gravity-3.3/LICENSE`. See
`../GRAVITY_REFERENCE.md` §7 for the obligations of the CUDA translation. The translated files are
`../cops_cuda/gravity.py` and `../cops_cuda/kernels/gravity.cu`; both carry the GPLv3 notice.

This is downloaded material. It is read for translation and comparison and is never
executed from this directory. The CPU reference package is installed from the
hash-locked wheel (`../requirements-cuda.lock`) or is the campaign's own binary.
