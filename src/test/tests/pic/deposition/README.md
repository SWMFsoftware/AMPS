# PIC deposition tests

`CORE-PIC-DEPOSITION-001` qualifies the dependency-light deposition contract
shared with `CornerBased::InitStencil`.  Coordinates are dimensionless cell
coordinates; corner `ix+2*iy+4*iz` receives the tensor-product CIC shape.
Outputs are represented nodal charge and `q*v`, before division by a mesh dual
volume.

The cell-crossing case uses a finite-volume segment current.  If `Q_i` is CIC
charge on node `i`, its edge current is built from the prefix of
`Q_i(new)-Q_i(old)`, giving
`dQ_i/dt + I_(i+1/2)-I_(i-1/2)=0` to roundoff even when the segment crosses
several cells.  This current is intentionally distinguished from ECSIM's
point-current/susceptibility deposition; the test does not claim that an
implicit point current is an Esirkepov current.

Run from the repository root:

```sh
python3 src/test/run_tests.py --test CORE-PIC-DEPOSITION-001
```
