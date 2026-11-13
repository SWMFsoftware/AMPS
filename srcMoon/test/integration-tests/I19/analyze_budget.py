#!/usr/bin/env python3
"""Independent signed Na/Ar accounting, not an I19 predecessor waiver.

Read actual production event weights and independently summed inventories.
Decimal arithmetic avoids losing small residuals by subtracting large double
values a second time. No observational values or theoretical source totals
are supplied to the ledger. Numerical deletion is never a physical escape.
"""
import argparse
import csv
from decimal import Decimal, localcontext
import json
from pathlib import Path


def evaluate(path):
    """Reconstruct every gas, surface and combined identity and count receipt."""
    with path.open(newline="") as stream:rows=list(csv.DictReader(stream))
    if not rows: raise ValueError("empty production M2 ledger")
    species={r["species"] for r in rows}
    if species!={"Na","Ar"}: raise ValueError("ledger must contain both Na and Ar")
    failures=[];maximum=Decimal(0)
    with localcontext() as context:
        context.prec=60
        for r in rows:
            d=lambda name: Decimal(r[name])
            if any(not d(k).is_finite() for k in r if k!="species"):
                raise ValueError("nonfinite M2 ledger value")
            count=d("gas_count")-d("initial_count")-d("emission_count")+d("adsorption_count")+d("photolytic_loss_count")-d("photolytic_product_count")+d("escape_count")+d("numerical_removal_count")
            gas=d("gas_atoms")-d("initial_gas_atoms")-d("emission_atoms")+d("adsorption_atoms")+d("photolytic_loss_atoms")-d("photolytic_product_atoms")+d("escape_atoms")+d("numerical_removal_atoms")
            surface=d("surface_atoms")-d("initial_surface_atoms")-d("surface_deposit_atoms")+d("surface_release_atoms")-d("reservoir_supply_atoms")
            bound=d("roundoff_bound_atoms")
            if bound<0 or d("mass_kg")<=0:raise ValueError("invalid ledger bound/mass")
            # CSV serialization contributes at most one last printed digit
            # per operand, independently bounded here, rather than making a
            # receipt mismatch look like a physical conservation failure.
            printed=sum(Decimal(10)**d(k).as_tuple().exponent for k in
                        ("gas_atoms","initial_gas_atoms","emission_atoms","adsorption_atoms",
                         "surface_atoms","initial_surface_atoms","surface_deposit_atoms","surface_release_atoms"))
            # The RECEIPT residual is itself computed in C++ long double.
            # Decimal reconstructs the serialized operands exactly, so its
            # difference from that receipt also needs the declared arithmetic
            # gamma_n bound, not merely CSV last-digit error. This tight
            # expression bound is separate from the much larger cumulative
            # double surface-update bound and the physical 0.5% criterion.
            epsilon=d("arithmetic_epsilon")
            if not 0<epsilon<Decimal("0.001"):raise ValueError("invalid arithmetic epsilon")
            operands=("gas_atoms","initial_gas_atoms","emission_atoms","adsorption_atoms",
                "photolytic_loss_atoms","photolytic_product_atoms","escape_atoms","numerical_removal_atoms",
                "surface_atoms","initial_surface_atoms","surface_deposit_atoms","surface_release_atoms","reservoir_supply_atoms")
            arithmetic=20*epsilon/(1-20*epsilon)*sum(abs(d(k)) for k in operands)
            if count!=d("count_residual"):raise ValueError("reported integer residual disagrees with events")
            if abs(gas-d("gas_residual_atoms"))>printed+arithmetic or abs(surface-d("surface_residual_atoms"))>printed+arithmetic:
                raise ValueError("reported weighted residual disagrees with events")
            if abs(gas+surface-d("combined_residual_atoms"))>printed+arithmetic:
                raise ValueError("reported combined residual disagrees with events")
            dominant=max(d("emission_atoms")-d("surface_release_atoms"),d("surface_deposit_atoms"),d("surface_release_atoms"),d("escape_atoms"),d("numerical_removal_atoms"))
            relative=abs(gas+surface)/dominant if dominant else Decimal(0)
            maximum=max(maximum,relative)
            passed=(count==0 and abs(gas)<=bound+printed and abs(surface)<=bound+printed and
                    abs(gas+surface)<=2*bound+printed and relative<Decimal("0.005") and
                    d("photolytic_loss_atoms")==0 and d("photolytic_product_atoms")==0 and
                    d("numerical_removal_count")==0)
            if not passed:failures.append({"species":r["species"],"step":r["completed_steps"],"count_residual":str(count),"gas_residual_atoms":str(gas),"surface_residual_atoms":str(surface)})
    finals={s:[r for r in rows if r["species"]==s][-1] for s in species}
    # Require actual coverage, not a vacuous zero-event conservation PASS.
    for s,r in finals.items():
        for e in ("emission","adsorption","surface_deposit","escape"):
            if int(r[e+"_count"])==0:failures.append({"species":s,"missing_event":e})
    if int(finals["Ar"]["surface_release_count"])==0:
        failures.append({"species":"Ar","missing_event":"surface_release"})
    return {"status":"FAIL" if failures else "PASS","scope":"supporting-I19-Na-Ar-accounting",
            "formal_i19_status":"SKIPPED","formal_reason":"I13 qualification is not supplied by this supporting check",
            "rows":len(rows),"maximum_relative_imbalance":str(maximum),"failures":failures,"final_rows":finals}


def render(path,output):
    """Plot unchanged production rows; label numerical evidence explicitly."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with path.open(newline="") as stream:rows=list(csv.DictReader(stream))
    fig,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for col,species in enumerate(("Na","Ar")):
        group=[r for r in rows if r["species"]==species]
        x=[int(r["completed_steps"]) for r in group]
        for field in ("gas_atoms","surface_atoms"):
            axes[0,col].plot(x,[float(r[field]) for r in group],label=field)
        for field in ("gas_residual_atoms","surface_residual_atoms","roundoff_bound_atoms"):
            axes[1,col].plot(x,[float(r[field]) for r in group],label=field)
        axes[0,col].set(title=species,ylabel="atoms")
        axes[1,col].set(xlabel="completed synchronized particle steps",ylabel="atoms",yscale="symlog")
        for ax in axes[:,col]:ax.legend()
    fig.suptitle("SUPPORTING M2 ACCOUNTING — NOT FORMAL I19 / NOT OBSERVATIONAL VALIDATION")
    for suffix in ("png","eps"):fig.savefig(output/f"budget.{suffix}",dpi=150)
    plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("ledger",type=Path);p.add_argument("--output-dir",required=True,type=Path)
    a=p.parse_args()
    if a.output_dir.exists() and any(a.output_dir.iterdir()):p.error("output directory must be new or empty")
    a.output_dir.mkdir(parents=True,exist_ok=True)
    result=evaluate(a.ledger);render(a.ledger,a.output_dir)
    (a.output_dir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2));raise SystemExit(0 if result["status"]=="PASS" else 1)
