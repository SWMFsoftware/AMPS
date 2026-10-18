# Stage 7 qualification contract

`stage7_qualification_plan.v1.json` freezes the project-owner decisions used by
the Stage 7 software implementation. It is an acceptance contract, not an
external measurement.

The committed AMS-02 and GOES-16 bytes remain under
`tests/reference_spectra/`. Stage 7 does not copy or reinterpret those files.
The AMS A=1 proton component is admitted only as an incomplete GCR component
test. The GOES fixture remains ingestion/time-integration only because neither
the sensor-to-simulation transform nor its angular response has been supplied.

The representative matrix deliberately reverses the same Al/HDPE layers and
uses independently defined natural-silicon and water slabs. Quantity
availability follows the already validated Stage 4--6 material/species
boundaries. In particular, proton irradiation cannot be reported as ASTM
neutron-equivalent fluence, and water has no Stage-5 damage response.

Response statistics use explicitly separated pilot/allocation phases. The
12-node base grid, 11 even-index refinement energies, and final 22 odd-index
energies receive separate independent 256-history pilots and frozen
allocations. Each of five production replicates receives exactly 73,728 base
histories, 56,320 first-refinement histories, and 112,640 exact-refinement
histories, with a 512-history cell floor (1,213,440 production histories per
stack). The final supplement was triggered because the completed 23-node
HDPE-to-Al comparison passed the 10% relative gate but returned `z=2.078`,
outside the unchanged `1.96` gate. Earlier outcomes cannot alter their stopping
counts after either refinement decision. The
larger of the propagated within-replicate and empirical between-replicate
standard errors is used for the direct-versus-matrix gate. This prevents a
single stream that under-samples the rare TID tail from claiming excessive
precision. Production seeds are disjoint and pilot moments are
excluded from the response estimate; the pilot therefore controls allocation
without becoming a data-dependent optional-stopping contribution. The
allocator cannot enable Geant4 physics biasing.

The conservative maximum-within/between rule applies to both sides of the
comparison. At the exact 45-node grid, the reversed-stack five direct seed
means had a between-seed standard error 20.8% larger than their pooled
event-level estimate. Applying the existing rare-tail rule symmetrically does
not change the mean or the registered 10%/95% gates; it prevents the direct
side from claiming more precision than its independent seeds demonstrate.

Full generated response cells are validated before they are compacted to the
two selected scores and corresponding covariance block. A content-hashed
completion marker preserves restart closure. This storage operation does not
alter score moments or the response key and is separate from acceptance.

`SHIELDSIM_STAGE7_WORKERS` changes only how many independent serial cases are
launched concurrently. It does not change this contract, any seed/history
count, or the deterministic evidence order. MPI and Geant4 event-level
multithreading are not used by this campaign.

Each direct comparison row pools five independent 8,000-history seeds. The
40,000-history total is required because finite-slab TID retains a heavy
nuclear/secondary tail; the earlier 4,000-history development pilot was not a
stable equivalence estimator and is not an accepted campaign size.

Direction uses three Gauss-Legendre nodes for the isotropic-plane `mu`
density and three Gauss-Legendre azimuth nodes on the square-symmetry interval
`[0, pi/4]`. A two-endpoint azimuth development rule differed from the
continuous-source HDPE→Al row by 12%; it was replaced by the higher-order rule
instead of weakening the registered 10%/95% gates.

All generated cases, response cells, folds, logs, numerical evidence, and
figures belong in the configured build/test-output tree. They are not reference
data and must not be committed as replacement acceptance values.
