# Decision boundary: regional cross-ray evolution

Date: 2026-11-10. **The user selected Path A: retain the analytical/semi-empirical
regional closure and develop conservative lateral material transport without
an MHD runtime.** This decision preserves the runtime prohibition.

This is a scope/closure decision document, not implementation
or qualification evidence. The completed bounded transport, geometry, CT,
fluid-update and corotating-contact kernels remain useful for either path.

The user explicitly retained the 3-D ambient, selected conservative cross-ray
transport, selected a field-threaded material contact, and selected a 1200-s
corotation-to-fixed-orientation transition. Those choices resolve ambient,
startup and magnetic-interface issues. They do not unambiguously select the
kind of regional evolution algorithm that may replace the radial tubes.

The authoritative addendum
[`BG3D4_CODEX_REVIEW_AND_CORRECTIONS.md`, section 6](BG3D4_CODEX_REVIEW_AND_CORRECTIONS.md#6-a-velocity-closure-menu-is-not-a-closed-physical-algorithm)
says **“Do not introduce an MHD runtime.”** The repository G1 contract also
calls for a prescribed analytical/semi-empirical model, rather than a global
MHD/eruption-instability solver. Full-CME resumption does not explicitly revoke
that algorithm boundary. The user explicitly selected the analytical/semi-empirical path; no exception
to that algorithm boundary is authorized.

## Path A: retain the analytical/semi-empirical regional-model boundary

Develop a complete lateral material/velocity closure, initialized from the
retained time-dependent ambient, and evolve its conservative mass/flux map.
The prescribed contact and compatible outer/admission states must determine
normal motion. Lateral transport, the mass distribution, thermodynamics and
magnetic mapping require a jointly specified closure, not just a menu of
velocity prescriptions. A metric-consistent surface/volume formulation and
its actual boundary conditions must close integral compatibility and avoid
invented mass sources or flank exits. Required force/heating/work remains
explicit and must meet the frozen physical residual gates independently.

The current HLL/CT/static-update kernels supply independent conservation
references and bounded transport tools. Their existence does not select an
MHD PDE time integrator as the production model. This path preserves the
existing prohibition and needs a documented analytical closure before new
production evolution is claimed. It may still fail validity or interface
qualification; a prescribed map is not accepted merely for conserving mass.

## Path B: explicitly authorize numerical 3-D MHD in the finite CME region

Evolve regional mass, all momentum components, total energy and CT magnetic
flux using a qualified MHD time integrator, driven by the one prescribed
contact. The eruption/driver is still prescribed. This is limited to the
sheath/compression/background region with explicit undisturbed buffers;
it is not a global solar-eruption MHD simulation or particle-feedback model.

This path requires an explicit exception to the addendum's MHD-runtime rule.
The next implementation would add compatible curved-face/time quadrature,
ideal upwind edge EMF, canonical multidimensional CFL control and conservative
moving geometry/remap with the actual contact. Then qualify unforced ambient
continuity/induction and its prescribed force/heating/work, controlled domain
boundaries/buffers, computed shock, original mass/time/ray refinement gates,
validity, assembly and the remaining G1/native stages. Generic PIC and legacy
providers would retain their existing algorithms and identities.

The currently implemented atomic static update is not this integrator: its
face flux and edge EMF are supplied inputs. It does not select those fluxes,
advance a production model autonomously, move the domain, or declare physical
boundary/source laws. No native runtime or new production selector is present.

## Common acceptance boundary

Neither path may bypass the original contact, divergence, conservation,
shock-branch, convergence, force/work, unsupported-area or runtime gates. The
frozen radial fixture and failure receipts remain reproducible. A distinct
expanded closure identity must precede production selection. BG3D-4 and G1
remain unqualified until actual required evidence exists; reduced-front and
new bounded kernel PASS results cannot substitute for it.
