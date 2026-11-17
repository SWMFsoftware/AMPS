# SWMF field-line producer contract audit

This is a source audit of a specific local SWMF/BATSRUS producer, not a
qualification of the coupled application. The user-defined local wave state
continues to store magnetic energy alone. No input is rescaled or assigned a
polarization from its column name.

The inspected clean producer files belong to local SWMF revision
`9258817b91d9a299ac8b046593f210ed115efdbf`, under
`/data/vtenishe/AMPSTEST/dev-sepswmf/AMPS/SWMF-IH/SWMF`.
The separate numerical MFLAMPA comparator is pinned to
`38e9defc9c525bd3d49669bc2e7d8014987de488`; that regression does not identify
the historical producer revision of every archived MH file.

## Branch mapping established for the inspected standard-energy path

`IH/BATSRUS/src/ModWaves.f90` assigns `AlfvenPlusFirst_ = WaveFirst_` and
the minus group to the last half of the wave variables.
`IH/BATSRUS/src/ModPhysicalFlux.f90` then uses
`U_normal + B_normal/sqrt(rho)` for the plus group and the corresponding
minus sign for the minus group. These are signs relative to the actual B,
including its polarity, rather than a permanent radial direction.

`IH/BATSRUS/srcInterface/IH_wrapper.f90::IH_get_for_mh` exports
`State_V(WaveFirst_:WaveLast_)` in that order with an energy-density unit
conversion. `CON/Coupler/src/CON_bline.f90::BL_put_from_mh` copies those
two entries to `Wave1_:Wave2_` without a branch swap. The inspected AMPS
`srcInterface/amps2swmf.cpp::amps_get_bline_c_` also copies the two entries
without reordering.

Thus, **when this exact standard-energy producer contract applies**, the
direction mapping is:

| Imported entry | Local propagation branch | Physical characteristic |
|---|---|---|
| Wave1 / WaveFirst | E+ | U + vA b |
| Wave2 / WaveLast | E− | U − vA b |

No additional swap follows from reversing B or the geometric line ordering.
A transformation between common Cartesian bases rotates B and U; it does
not change which wave propagates along B. Circular polarization and spectral
shape are not supplied by these two integrated entries.

This mapping does **not** establish the semantics of arbitrary files called
Wave1/Wave2. BATSRUS also supports representative variables and boundary
thread major/minor amplitudes. `ModFieldLineThread.f90` explicitly converts
major/minor order using the sign of B dot x. Such private arrays must not be
imported as though they were the standard energy variables. The active
producer mode, its conversion path, and the exact exported schema must be
recorded before an automatic adapter can use the table above.

## Energy normalization remains an explicit interface issue

The adopted srcSEP magnetic-energy conversion is `deltaB² = 2 mu0 W_B`.
The inspected producer's `ModTurbulence.f90` around line 371 explicitly
identifies half its wave energy as the magnetic oscillation energy. Its
standard variables therefore carry total Alfvén-wave energy in that mode.
The pinned numerical MFLAMPA implementation instead computes
`deltaB² = mu0 (Wave1 + Wave2)` in `src/ModDiffusion.f90`. These conversions
are incompatible if the same numeric columns are interpreted as magnetic
energy in both programs. The numerical comparator's energy convention must
therefore be translated through an explicit producer/consumer contract; the
user's local magnetic-energy definition is not changed to accommodate it.

Units of J/m³ alone do not distinguish magnetic energy from total Alfvén-wave
energy. A future adapter must state that distinction and its conversion in
metadata. It must also reject absent L/R spectra when a polarization-resolved
law is requested; assigning equal L/R is authorized for the declared **local
initial condition**, not for an arbitrary imported observation or state.

The user selected this standard BATSRUS total-energy producer after reviewing
the discrepancy. Its raw entries therefore require an explicit factor of one
half in the local magnetic-energy view. Keep raw imported arrays immutable
and record the conversion in provenance. A customized magnetic-only export
would require a separate declared contract; the implementation must not guess
between them from a Wave1/Wave2 label.

## Density and temperature findings

The inspected `BL_put_from_mh` exports `Rho = rho_SI / m_proton`, labelled
amu/m³ by its source. Its temperature is derived from pressure and mass
density using the coupler's declared `EnergySi2Io`. The default is
`1/k_B`; the current PT wrapper calls `BL_set_grid` without overriding it,
so that inspected live PT path carries kelvins. In contrast, the archived
MFLAMPA files explicitly carry a temperature energy in keV and their reader
performs its own conversion. A file and the live PT interface cannot share
an implicit temperature conversion based only on the name T.

The current AMPS importer copies Rho directly into a number-density slot and
forms pressure as `Rho * k_B * T`. This is a proton-number-density
interpretation only under an explicit composition and pressure-species
contract. No general mixture, electron pressure, proper density or missing
velocity-gradient tensor follows from those columns. The full imported
focused-transport provider remains unsupported without its required tensor.
This audit does not enable that provider or declare the copied primitives
physically qualified.

## Reproducible source identity

| Source file | SHA256 |
|---|---|
| IH/BATSRUS/src/ModWaves.f90 | `d3106b0ef5538bb7f11b1397516cb03eeefa916a8ecd44a6af7677ed050d8cf6` |
| IH/BATSRUS/src/ModTurbulence.f90 | `587353f50376376f6fb3091b3a12011810faced5ba3020494ab1c55f72c1256f` |
| IH/BATSRUS/src/ModPhysicalFlux.f90 | `6edd514d58bc79d7904a00a9f41be9fcba6e54c850029f5b8e36c5e56b50a7f3` |
| IH/BATSRUS/srcInterface/IH_wrapper.f90 | `7ee30736e7a5808a6ba385a677d142a3ebff7ab9cc6b4885dd7a7a46c9f6962e` |
| IH/BATSRUS/src/ModFieldLineThread.f90 | `c254b29168f4aa42878e55dd4ac1c8525a993461e34b1a4163461474166784c6` |
| CON/Coupler/src/CON_bline.f90 | `0a70b260dd6660eb838cb027c2f572c980d10e58879bea4fd96c5e8398f543a0` |
| AMPS srcInterface/amps2swmf.cpp | `22404e5e9a0db280c839d08e87e91778464e9c0de6a7cfbce09fbb909f2f7549` |
| AMPS srcInterface/PT_wrapper.f90 | `ffc6ba33184b0a7330a40559998c438bf73b921d4df0a8f056186196fd94620e` |

See [STATUS.md](../../STATUS.md) for implementation gates and
[reference_sources.json](../reference_sources.json) for the numerical
MFLAMPA archive and comparator provenance.
