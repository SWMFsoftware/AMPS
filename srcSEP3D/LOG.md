commit: e4a4d44668121d1a37655e5bf67fd74f21472918
PBS r1i7n7 176> mpiexec -n 10 ./amps --input srcSEP3D/examples/sep3d_analytic_parker.in --initialization-only --initialization-output-dir sep3d_mesh_preview 

initialization_parker_line=sep3d_mesh_preview/sep3d-initialization-parker-line.dat

commit: fb53a819fc449b39d52d20ac3ca0cfca34560500 
mpiexec -n 10 ./amps --input srcSEP3D/examples/sep3d_analytic_parker_active_tube.in --initialization-only --initialization-output-dir sep3d_mesh_preview

