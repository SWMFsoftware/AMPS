# Stage 6 LET reference data

The proton and He-4 reference matrix reuses the checksum-pinned NIST
PSTAR/ASTAR tables under `tests/reference_tid/nist`. The selected heavy-ion
screen is fully stripped C-12 in natural silicon at 10, 100, and 1000 MeV/u.

The three raw HTML responses were downloaded from the official GSI ATIMA web
calculator on 2026-10-17 using ATIMA v1.41's published endpoint:

`https://web-docs.gsi.de/~weick/atima/atima14_script.cgi?zp=6&massp=12&energy=E&Zt=14&at=0&thickness=0.001`

The web result labels its computational backend “Atima 1.4.” The CSV converts
the reported MeV/(mg/cm2) values to MeV cm2/g by multiplying by 1000 and
converts energy per nucleon to total C-12 kinetic energy by multiplying by 12.
The raw responses are authenticated by `manifest.sha256`; no interpolation is
performed in the registered test.

No redistribution license was identified on the calculator page. These small
machine responses are retained only as pinned validation evidence with source
and query provenance. They are not a production dependency. The comparison is
a 15% implementation/reference screen, not an uncertainty statement and not
experimental validation.
