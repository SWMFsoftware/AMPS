  Run the qualified case with:

  mpiexec -n 4 ./amps \
    --input srcSEP/examples/reduced-shock/field_line.in \
    --reduced-shock-event srcSEP/examples/reduced-shock/positive_1au.event \
    --total-iterations 208 \
    --coupling off --cascade off --reflection off --shock-injection off

run test:
srcSEP/test/run_step1_tests.sh

