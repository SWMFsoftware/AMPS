# DIM1 comet-result visualization

`plot_1d.py` reads the quoted variable names in AMPS Tecplot output instead of
assuming fixed column positions. It plots one selected physical quantity from
one or more files on common axes. The X axis is always the logarithmic
cometocentric radius. The Y scale must be selected explicitly as `linear`,
`log`, or `symlog`.

The script requires Python 3, NumPy, and matplotlib. Saving with `--output`
selects matplotlib's non-interactive backend automatically, so no display is
needed on a batch node.

List the available physical quantities:

```sh
python3 srcHartley2/visualization/plot_1d.py \
  PT/plots/pic.H2O.s=13.out=0.dat --list-quantities
```

Plot water number density with logarithmic axes and radius in kilometres:

```sh
python3 srcHartley2/visualization/plot_1d.py \
  PT/plots/pic.H2O.s=13.out=0.dat \
  -q dens --yscale log --radius-unit km --rmax 1000 \
  --output PT/plots/h2o_number_density
```

Compare the same physical quantity from several output files:

```sh
python3 srcHartley2/visualization/plot_1d.py \
  run-a/pic.H2O.s=13.out=0.dat run-b/pic.H2O.s=13.out=4.dat \
  -q temp --yscale linear --output temperature-comparison
```

`-q` is the short form of `--quantity`. The standard aliases are:

| Alias | AMPS quantity |
| --- | --- |
| `r` | Radius `X` |
| `level`, `id`, `rank` | AMR level, temporary mesh ID, and MPI rank/thread |
| `dens`, `part` | Number density and model-particle number |
| `vel` or `vx`, `vy`, `vz` | Radial/X, Y, and Z velocity components |
| `speed` | Velocity magnitude `\|V\|` |
| `temp`, `tpar`, `ttan` | Total, parallel, and tangential translational temperature |
| `dt`, `weight` | Local timestep and local particle weight |
| `coll-SPECIES` | Collision frequency with `SPECIES`, such as `coll-h2o` or `coll-co2` |

Run `--help` for the complete species-specific alias list and copyable examples
for listing quantities, plotting one file, comparing runs, and batch output.
Exact quantity names remain case-insensitive, and an unambiguous prefix is
accepted. Ambiguous selections are rejected and report the matching columns.

Limit the plotted cometocentric distance with one of these mutually exclusive
controls:

- `--max-distance MAX` or `--rmax MAX` plots from the available inner radius
  through `MAX`;
- `--distance-range MIN MAX` or `--rrange MIN MAX` plots an inclusive radial
  interval.

All distance limits use the unit selected by `--radius-unit`. For example,
`--radius-unit km --rrange 6 100` selects cometocentric radii from 6 km through
100 km. Both limits must be positive and a range must have `MIN < MAX`.

Each curve's legend label is the input file's basename. With logarithmic Y
scale, the script skips non-positive points, reports their count on standard
error, and breaks the curve at each skipped point. Use
`--yscale symlog --linthresh VALUE` when zero or signed values must remain
visible. Without a distance limit, the full radial extent remains visible when
logarithmic Y scaling skips zero-valued outer cells.

If a file has no positive values left for a logarithmic Y plot, its entire
curve is skipped with a diagnostic and the remaining input files are still
plotted. The command fails only when none of the requested files has plottable
data. This permits shell patterns such as `PT/plots/*.dat` when some species
are absent in an early sampling window.

`--output NAME` always saves `NAME.png` and `NAME.eps`. If the argument has a
different extension, such as `.pdf`, that format is saved in addition to PNG
and EPS. After saving, the script opens the PNG in the system's desktop image
viewer. Use `--no-show` for a batch job. When no graphical display is
available, the script reports that condition, saves the files, and skips the
viewer automatically. If `--output` is omitted, matplotlib's interactive
viewer is used directly.
