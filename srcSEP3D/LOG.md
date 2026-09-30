commit: e4a4d44668121d1a37655e5bf67fd74f21472918
Full domain:
PBS r1i7n7 176> mpiexec -n 10 ./amps --input srcSEP3D/examples/sep3d_analytic_parker.in --initialization-only --initialization-output-dir sep3d_mesh_preview 

Refined corridor:
mpiexec -n 10 ./amps --input srcSEP3D/examples/sep3d_analytic_parker_active_tube.in --initialization-only --initialization-output-dir sep3d_mesh_preview

Result:
initialization_parker_line=sep3d_mesh_preview/sep3d-initialization-parker-line.dat


Run tests:
srcSEP3D:
PBS r3i3n3 124> mpiexec -n 4 ./amps --test SCCM3D01 --test SCCM3D02 --test SCCM3D03 --test SCCM3D04 --test SCCM3D05 --test SCCM3D06 --test SCCM3D07 --test-input srcSEP3D/examples/sep3d_analytic_parker.in --test-steps 0 --expect-mpi-ranks 4 --test-json test_output/coupled-sep-corona/native.json --artifact-directory test_output/coupled-sep-corona/artifacts


env MAKEFLAGS="-j16" test/run_tests.py --all --amps-source .. --make-config ../Makefile.conf --output-dir test_output/all --rebuild


srcSEP:
test/run_tests.py --amps ../amps --all --output-dir test_output/all
