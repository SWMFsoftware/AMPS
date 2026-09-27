commit: e4a4d44668121d1a37655e5bf67fd74f21472918
PBS r1i7n7 176> mpiexec -n 10 ./amps --input srcSEP3D/examples/sep3d_analytic_parker.in --initialization-only --initialization-output-dir sep3d_mesh_preview 

initialization_parker_line=sep3d_mesh_preview/sep3d-initialization-parker-line.dat


commit: 0511c745a7afe138ea251753b636d1340162f6f2 
a part of the domain outside of the a corridor around a parker spiral is turned off -- still there are some issues
mpiexec -n 10 ./amps --input srcSEP3D/examples/sep3d_analytic_parker_active_tube.in  --initialization-only --initialization-output-dir sep3d_mesh_preview
 
