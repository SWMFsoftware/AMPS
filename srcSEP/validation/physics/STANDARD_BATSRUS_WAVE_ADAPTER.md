# Standard BATSRUS wave-energy adapter

The user selected **standard BATSRUS total Alfvén-wave energy** exports. This
overrides the earlier SWMF-raw-is-magnetic statement; local evolved magnetic
energy remains unchanged. The exact producer audit, revision and source hashes
are in [SWMF_FIELD_LINE_CONTRACT_AUDIT.md](SWMF_FIELD_LINE_CONTRACT_AUDIT.md).
The conversion here is application-owned; no shared PIC storage or other
application is rescaled.

The standard two integrated fields have the adopted meaning:

| Raw entry [J/m³] | Private magnetic view [J/m³] | Direction |
|---|---|---|
| Wave1 total Alfvén energy | E+=Wave1/2 | along actual B |
| Wave2 total Alfvén energy | E−=Wave2/2 | opposite actual B |

Magnetic fluctuation variance uses the existing common magnetic helper:
`deltaB±²=2 mu0 E±`. Independently, the producer totals therefore imply
`deltaB²=mu0(Wave1+Wave2)`. No swap follows from reversing magnetic polarity or
line-coordinate order. Representative amplitudes/private major-minor arrays
are different producer contracts and must not be guessed from a column name.
The two totals supply neither L/R handedness nor a spectral shape. The user's
equal **initial local** L/R choice does not supply that missing imported state.

`util/sep_standard_batsrus_waves.{h,cpp}` defines typed
`RawTotalAlfvenEndpoint` and `MagneticBranchView`. It validates both raw
endpoints, even a zero-weight one, before making their half-energy views and
convex interpolation. Finite nonnegative raw energy is required. Positive
half-energy, interpolated energy or canonical variance that underflows returns
`StepUnderflow`; invalid/nonfinite data returns a typed error. No clipping,
flooring, extrapolation, endpoint masking, L/R partition or raw mutation occurs.
Wider convex arithmetic avoids overflowing an unnecessary raw endpoint sum.
Output replacement happens only after every check succeeds.

The contract string has one authority, `ContractIdentity()`. Its version,
both complete raw endpoints, interpolation fraction and nonzero generation
bind a canonical scalar hash. `Representation()` states the total-to-magnetic
half conversion and B-directed labels. `circularPolarizationsKnown=false`
remains explicit. The SWMF background fingerprint and coefficient source
identity include this contract; a model cannot resume an earlier convention
while presenting the same metadata.

`standard_batsrus_wave_adapter.{h,cpp}` supplies
`Transport::PICAdapter::ReadStandardBatsrusWaveView`. It reads the actual
allocated `DatumAtVertexPlasmaWaves` arrays through const pointers. Access
requires an active immutable `Swmf/ImportedReadOnly` snapshot, exact context
generation, valid in-range coordinate/matching segment and two components at
both vertices. It calls the same pure converter used by tests. The production
SWMF branch of `coefficient_providers.cpp` uses this accessor, its known
directions and its representation/checksum. Other coefficient sources retain
their magnetic-energy conversion. The overall imported plasma provider still
rejects its unresolved gradient/composition/frame contracts; this bounded
wave adapter does not weaken those gates.

Run the strict-warning/ASan/UBSan dependency-light factory with:

```sh
SRCSEP_REPORT_DIR=/tmp/sep-standard-batsrus sh srcSEP/test/run_standard_batsrus_wave_tests.sh
```

SBW01–05 in group `standard-batsrus-waves` pass **5/5**: independent
`mu0*rawSum` reference/local-unchanged control, endpoint reversal, invalid and
underflowing endpoints, interpolation domain/representability and contract/
generation/raw identity with absent L/R. The same callbacks can run natively.

SBWN01–02 in group `standard-batsrus-waves-native` pass **2/2 on two MPI ranks
in build8**; SBW01–05 also pass natively in that build. The binary SHA256 is
`049a1959918dd154e5a81d75ed7ebb18a57764f29c09c3917c81b4068f9c58eb`;
STATUS tracks its captured source point and structured reports. They use
actual allocated native vertex arrays and the
production accessor under explicit synthetic importer-owned read-only
generations. They check polarity reversal, endpoint/fraction reversal,
foreign-provider/read-phase/generation rejection, invalid unweighted endpoints
and byte-exact nonmutation of the complete vertex associated buffers. They
restore the original complete buffers and exact immutable snapshot pointer;
PB counts and legacy RNG/cache remain unchanged. This is a native accessor
fixture, not a live BATSRUS coupling or observational comparison.

Still separate: live imported field-line plasma/gradient/composition/temperature
semantics, full guiding-center source/particle frames, a four-spectrum imported
authority, material geometry and an exactly-once local-evolution handoff.
The coefficient radius for offset analytic origins and guessed unused nucleon
metadata are tracked independent followups; no frame or isotope metadata is
invented in this bounded wave conversion.
