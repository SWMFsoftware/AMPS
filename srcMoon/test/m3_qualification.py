#!/usr/bin/env python3
"""Approved M3 planning and independent cohort scoring, never synthetic data.

Pure scorers evaluate explicitly supplied original-parent records. They do not
turn an analytical fixture, an unbound CSV, or an approval boolean into formal
qualification. The linked owner must additionally bind actual production
execution, source/state providers, independent numerical/budget evidence and
same-campaign predecessors. The original preparation protocol stays untouched.
"""
from __future__ import annotations
import csv
from decimal import Decimal,InvalidOperation,localcontext
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path
import sys

TEST=Path(__file__).resolve().parent
ADDENDUM=TEST/'M3_QUALIFICATION_ADDENDUM_20261118.json'
APPROVED_SHA256='551bb63b5de5be5766137bafb456aa140ae70fba5c38b4d9d8e9f70e74c15dbd'
BASE_SHA256='29e183978e4c336891ea957726ddeb32e62898e7cb9eab8181909910b520ffaa'
ENSEMBLE_LIMIT={'I21':.02,'I22':.02,'I23':.03,'I24':.03}
NUMERICAL_LIMIT={'I21':.002,'I22':.002,'I23':.003,'I24':.003}
PARENT_COLUMNS=('original_parent_id','status','relative_event_time_s','channel_id','initial_represented_number')


class MissingPrerequisite(RuntimeError):
    """A genuine missing dependency/proof, not corrupt source or failed math."""


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest_syntax(value):
    return isinstance(value,str) and len(value)==64 and all(c in '0123456789abcdef' for c in value)


def strict_object(path):
    def pairs(items):
        result={}
        for key,value in items:
            if key in result:raise ValueError('duplicate JSON key: '+key)
            result[key]=value
        return result
    def nonfinite(value):raise ValueError('nonfinite JSON token: '+value)
    value=json.loads(Path(path).read_text(),object_pairs_hook=pairs,parse_constant=nonfinite)
    if not isinstance(value,dict):raise ValueError('JSON document must be an object')
    return value


def finite(value,label,positive=False):
    if isinstance(value,bool):raise ValueError(label+' cannot be boolean')
    try:number=float(value)
    except (TypeError,ValueError,OverflowError):raise ValueError(label+' must be a representable finite number')
    if not math.isfinite(number) or (positive and number<=0):raise ValueError(label+' must be finite'+(' and positive' if positive else ''))
    if isinstance(value,str) and number==0:
        try:
            if Decimal(value)!=0:raise ValueError(label+' underflows binary64; it cannot silently become zero')
        except InvalidOperation:raise ValueError(label+' is not a finite decimal number')
    return number


def integer(value,label,positive=True):
    if type(value) is not int or value<(1 if positive else 0):raise ValueError(label+' must be an explicit '+('positive' if positive else 'nonnegative')+' integer')
    return value


def approved_protocol(path=ADDENDUM):
    """Authenticate the reviewed addendum AND unchanged legacy contract bytes.

    The pinned digest is a code-reviewed identity, not an editable readiness
    flag. Scientific/native prerequisites remain independent even after this
    protocol approval. A missing file is SKIPPED; changed bytes are ERROR.
    """
    path=Path(path)
    if not path.is_file():raise MissingPrerequisite('approved M3 addendum is absent')
    if sha(path)!=APPROVED_SHA256:raise ValueError('approved M3 addendum identity changed')
    value=strict_object(path)
    base=TEST/'M3_QUALIFICATION_PROTOCOL.json'
    if not base.is_file():raise MissingPrerequisite('original frozen M3 protocol is absent')
    if sha(base)!=BASE_SHA256:raise ValueError('original frozen M3 protocol bytes changed')
    if value['approval']['status']!='APPROVED':raise ValueError('approved addendum status contradicts its frozen identity')
    return value


def clopper_pearson(successes,trials,alpha):
    """Exact two-sided binomial interval, preserving zero/full-count limits.

    SciPy supplies the inverse regularized beta function; no Gaussian/Poisson
    replacement silently changes rare-channel coverage. Endpoint formulas use
    expm1/log to avoid subtractive cancellation. Alpha is the already allocated
    per-claim risk, not the family risk selected after observed counts.
    """
    integer(trials,'binomial trials');integer(successes,'success count',False)
    if successes>trials:raise ValueError('success count exceeds independent trials')
    alpha=finite(alpha,'allocated alpha',True)
    if alpha>=1:raise ValueError('allocated alpha must be below1')
    if successes==0:return (0.,-math.expm1(math.log(alpha/2)/trials))
    if successes==trials:return (math.exp(math.log(alpha/2)/trials),1.)
    try:from scipy.special import betaincinv
    except ImportError as error:raise MissingPrerequisite('SciPy exact beta-quantile dependency unavailable') from error
    low=float(betaincinv(successes,trials-successes+1,alpha/2))
    high=float(betaincinv(successes+1,trials-successes,1-alpha/2))
    if not (math.isfinite(low) and math.isfinite(high) and 0<=low<=high<=1):
        raise ArithmeticError('exact binomial interval is not representable')
    return low,high


def finite_grid_tail_bound(maximum_draws):
    """Report finite uniform support; NOT a universal model-error guarantee.

    The support omits2^-53 at each end of an ideal continuous uniform. A
    union bound concerns ONLY these omitted tails, not every observable or
    chaotic trajectory error. The qualification owner must separately admit
    the actual maximum draw family/numerical effect before physical claims.
    """
    integer(maximum_draws,'declared maximum draw family',False)
    smallest=2**-53
    with localcontext() as context:
        context.prec=80
        tail=min(Decimal(1),Decimal(maximum_draws)/(Decimal(2)**52))
    return {'maximum_draws':maximum_draws,'minimum_open_uniform':smallest,
        'maximum_open_uniform':1-smallest,'minimum_optical_depth':-math.log1p(-smallest),
        'maximum_optical_depth':-math.log(smallest),'omitted_end_support_mass_per_draw':2**-52,
        'omitted_tail_union_upper':float(tail),'omitted_tail_union_upper_decimal':str(tail),
        'scope':'finite-grid tail support only; not a complete RNG/physics error or independence proof'}


def survival_plan(rate_s,horizon_s,relative_gate,comparisons,risk,estimator):
    """Conservative independent-parent Hoeffding planning at ACTUAL lambda*T.

    Rate and lifetime gates have DIFFERENT accepted survival intervals. This
    planner never substitutes assumed exact log(2) for represented rate/time.
    The family union bound requires no independence between matched dt/MPI
    comparisons, while trials WITHIN each cohort remain independent parents.
    """
    rate=finite(rate_s,'reference rate',True);time=finite(horizon_s,'physical horizon',True)
    gate=finite(relative_gate,'relative ensemble gate',True);risk=finite(risk,'planning risk',True)
    integer(comparisons,'frozen comparison count')
    if gate>=1 or risk>=1:raise ValueError('relative gate and planning risk must be below1')
    a=rate*time
    if not math.isfinite(a) or a<=0:raise ValueError('reference cumulative hazard is unrepresentable')
    probability=math.exp(-a)
    if estimator=='rate':low,high=math.exp(-a*(1+gate)),math.exp(-a*(1-gate))
    elif estimator=='lifetime':low,high=math.exp(-a/(1-gate)),math.exp(-a/(1+gate))
    else:raise ValueError('primary estimator must explicitly be rate or lifetime')
    margin=min(probability-low,high-probability)
    if not math.isfinite(margin) or margin<=0:raise ArithmeticError('survival precision cannot be resolved at this represented horizon')
    needed=math.log(2*comparisons/risk)/(2*margin*margin)
    if not math.isfinite(needed):raise ArithmeticError('planned cohort size is unrepresentable')
    count=8*math.ceil(needed/8)
    return {'independent_original_parents':count,'reference_survival_probability':probability,
        'represented_reference_cumulative_hazard':a,'accepted_survival_interval':[low,high],
        'family_planning_risk':risk,'comparison_count':comparisons,'primary_estimator':estimator,
        'method':'two-sided Hoeffding sufficient family bound; no reduced physical acceptance gate'}


def branch_precision_plan(probability,reactions,claims,confidence=.95,relative_halfwidth=.10,planning_risk=.01):
    """Check PREDECLARED reaction count against exact interval-width planning.

    This function does not guess a maximum cohort or drop a rare branch. It
    returns an exact-binomial high-probability count range and evaluates EVERY
    possible count in that range with the actual exact interval. Explicit
    resource controls belong to the caller; the vector size is reported before
    allocation. A large design must be sent to a cluster rather than shortened.
    """
    p=finite(probability,'source branch probability',True)
    integer(reactions,'planned independent removals');integer(claims,'frozen branch claims')
    if p>1:raise ValueError('source branch probability exceeds1')
    confidence=finite(confidence,'family branch confidence',True)
    precision=finite(relative_halfwidth,'relative interval halfwidth',True)
    risk=finite(planning_risk,'precision planning risk',True)
    if confidence>=1 or risk>=1:raise ValueError('confidence/risk must lie strictly below1')
    try:from scipy.stats import binom
    except ImportError as error:raise MissingPrerequisite('SciPy exact binomial planning dependency unavailable') from error
    # Exact source-distribution quantiles avoid a rare-branch Gaussian law or
    # a needlessly p^-2 Hoeffding allocation. Each omitted tail has probability
    # at most risk/(2*claims); discreteness makes this conservative.
    tail=risk/(2*claims)
    lower=int(binom.ppf(tail,reactions,p));upper=int(binom.isf(tail,reactions,p))
    if lower>upper:raise ArithmeticError('planned binomial count envelope is unrepresentable')
    alpha=(1-confidence)/claims
    return {'reference_probability':p,'independent_removals':reactions,'allocated_alpha':alpha,
        'count_envelope':[lower,upper],'counts_to_certify':upper-lower+1,
        'relative_halfwidth_limit':precision,'precision_probability_at_least':1-risk/claims,
        'status':'REQUIRES_EXACT_WIDTH_CERTIFICATION','method':'exact binomial count envelope + exact CP widths, no rare-channel omission'}


def certify_branch_precision(plan,maximum_counts_to_check):
    """Certify every width in a planned envelope; caller sets work ceiling."""
    integer(maximum_counts_to_check,'explicit planning work ceiling')
    if plan['counts_to_certify']>maximum_counts_to_check:
        raise MissingPrerequisite('exact branch precision planning exceeds declared work ceiling; use explicit cluster plan, do not reduce branch family')
    p=plan['reference_probability'];n=plan['independent_removals'];maximum=0.
    for k in range(plan['count_envelope'][0],plan['count_envelope'][1]+1):
        low,high=clopper_pearson(k,n,plan['allocated_alpha']);maximum=max(maximum,(high-low)/(2*p))
    return dict(plan,status='PASS' if maximum<=plan['relative_halfwidth_limit'] else 'FAIL',
        maximum_certified_relative_halfwidth=maximum)


def complete_cohort_plan(test_id,rate_s,horizon_s,branch_probabilities,rate_comparisons,
                        branch_claims,minimum_removals,maximum_parents,maximum_search_steps):
    """Sufficient fixed-N plan for rate AND rare-channel interval precision.

    Every resource/search limit and complete family size is caller-declared.
    Uniform Bonferroni planning allocation covers one endpoint-error condition,
    one removal-count envelope per cohort, and each branch-count envelope.
    Exact binomial quantiles and monotonic CP endpoints bound every possible
    interval within those envelopes, without assuming Gaussian rare counts or
    independent daughters. Doubling supplies a sufficient design, not a claim
    of minimum N. Exceeding an explicit cap requests a cluster plan, never a
    smaller physical criterion. The original final acceptance is unchanged.
    """
    if test_id not in ENSEMBLE_LIMIT:raise ValueError('existing I21-I24 ID required')
    integer(rate_comparisons,'all frozen endpoint comparisons');integer(branch_claims,'all frozen branch claims')
    integer(minimum_removals,'declared minimum removal events',False)
    integer(maximum_parents,'explicit parent ceiling');integer(maximum_search_steps,'explicit planning search limit')
    probabilities=[finite(p,'complete source branch probability') for p in branch_probabilities]
    probability_sum=math.fsum(probabilities)
    # Source rates/total are represented binary64 ratios. Their complete
    # mathematical family can sum to the adjacent float below1 (real D04
    # molecular families do). Permit only the derived per-ratio rounding
    # envelope, RETAIN every original probability and report the discrepancy;
    # this is neither source renormalization nor a physical gate tolerance.
    rounding_bound=len(probabilities)*sys.float_info.epsilon
    if not probabilities or any(p<0 or p>1 for p in probabilities) or abs(probability_sum-1)>rounding_bound:
        raise ValueError('complete nonnegative source branching must sum exactly to1; no renormalization')
    if branch_claims<len(probabilities):raise ValueError('frozen claim family omits source branches')
    try:from scipy.stats import binom
    except ImportError as error:raise MissingPrerequisite('exact binomial planning dependency unavailable') from error
    d=ENSEMBLE_LIMIT[test_id];b=NUMERICAL_LIMIT[test_id]
    # (1+b)*(1+s)-1=d: reserve the APPROVED numerical bias before
    # stochastic planning; this is not a tighter final acceptance criterion.
    stochastic_margin=(d-b)/(1+b)
    conditions=2*rate_comparisons+branch_claims;allocated=.01/conditions
    initial=survival_plan(rate_s,horizon_s,stochastic_margin,1,allocated,
        'rate' if test_id=='I22' else 'lifetime')
    count=initial['independent_original_parents'];removal_p=1-initial['reference_survival_probability']
    alpha=.05/branch_claims;history=[]
    for _ in range(maximum_search_steps):
        if count>maximum_parents:break
        tail=allocated/2
        low=int(binom.ppf(tail,count,removal_p));high=int(binom.isf(tail,count,removal_p))
        widths=[]
        for p in probabilities:
            if p==0:widths.append(None);continue
            if low==0:widths.append(math.inf);continue
            # CP intervals are exactly symmetric under successes/failures.
            # Using min(p,1-p) avoids destroying dominant-branch precision by
            # bounding two highly correlated near-total counts separately.
            rare=min(p,1-p)
            if rare==0:
                span=clopper_pearson(0,low,alpha)[1]
            else:
                k_low=int(binom.ppf(tail,count,removal_p*rare))
                k_high=int(binom.isf(tail,count,removal_p*rare))
                if k_high>low:span=1.
                else:
                    # At fixed successes, increasing trials shifts BOTH beta
                    # quantiles down; increasing successes shifts them up.
                    # Opposite rectangle corners therefore enclose every
                    # possible conditional interval in these count ranges.
                    lower=clopper_pearson(k_low,high,alpha)[0]
                    upper=clopper_pearson(k_high,low,alpha)[1]
                    span=upper-lower
            widths.append(span/(2*p))
        adequate=low>=minimum_removals and all(w is None or w<=.10 for w in widths)
        history.append({'original_parents':count,'removal_count_envelope':[low,high],
            'relative_halfwidth_upper_bounds':[None if w is None or not math.isfinite(w) else w for w in widths],
            'adequate':adequate})
        if adequate:
            return {'status':'PASS','scope':'pre-execution sufficient planning only; no native/physical qualification',
                'independent_original_parents':count,'allocated_planning_risk':allocated,
                'total_declared_planning_conditions':conditions,'family_planning_risk':.01,
                'stochastic_relative_planning_margin':stochastic_margin,
                'unchanged_final_relative_gate':d,'approved_numerical_bias_budget':b,
                'branch_interval_alpha':alpha,'history':history,'endpoint_plan':initial,
                'source_probabilities':probabilities,'source_probability_sum':probability_sum,
                'source_probability_roundoff_bound':rounding_bound,'source_renormalization':False}
        count*=2
    return {'status':'SKIPPED','scope':'planning only; no physical test ran',
        'reason':'required fixed cohort exceeds explicit planning resources; declare cluster package without shortening criteria',
        'next_candidate_original_parents':count,'history':history}


def validate_plan(plan):
    """Validate an explicit original-parent cohort and its permitted estimator."""
    approved_protocol()
    if not isinstance(plan,dict) or plan.get('schema')!='srcMoon-M3-cohort-plan-v1':raise ValueError('unsupported M3 cohort plan')
    if plan.get('approved_addendum_sha256')!=APPROVED_SHA256:raise ValueError('cohort plan is not bound to the approved protocol')
    identity=plan.get('test_id')
    if identity not in ENSEMBLE_LIMIT:raise ValueError('cohort plan requires existing I21-I24 identity')
    integer(plan.get('original_parents'),'original independent parent count')
    finite(plan.get('initial_represented_number'),'initial represented number',True)
    finite(plan.get('horizon_s'),'common physical horizon',True)
    # The original v1 format omitted this field for the homogeneous case.
    # Preserve that receipt contract while requiring explicit opt-in to a
    # different estimator and treating unknown declared kinds as an error.
    kind=plan.get('test_kind','homogeneous-fixed-position')
    if kind not in ('zero-hazard-invariant','homogeneous-fixed-position','trajectory-integrated-hazard'):
        raise ValueError('unsupported cohort test_kind')
    zero=kind=='zero-hazard-invariant'
    varying=kind=='trajectory-integrated-hazard'
    if varying and identity!='I21':
        raise ValueError('the explicit straight-line exposure oracle belongs to I21')
    total=finite(plan.get('reference_total_rate_s'),'independent total rate',not zero)
    if zero and total!=0:raise ValueError('zero-hazard invariant requires exactly zero reference rate')
    integer(plan.get('branch_claims'),'whole frozen branch family count')
    if plan.get('seed') not in (731993,104729,130363) or type(plan.get('seed')) is not int:raise ValueError('cohort seed not in approved seed family')
    if not isinstance(plan.get('cohort_id'),str) or not plan['cohort_id']:raise ValueError('explicit semantic cohort ID required')
    channels=plan.get('source_channel_rates_s')
    if not isinstance(channels,list) or not channels:raise ValueError('complete declared source channel list required')
    identifiers=[];rates=[]
    for channel in channels:
        if not isinstance(channel,dict) or not isinstance(channel.get('id'),str) or not channel['id']:raise ValueError('nonempty source channel identity required')
        identifiers.append(channel['id']);rate=finite(channel.get('rate_s'),'source partial rate')
        if rate<0:raise ValueError('negative source partial rate')
        rates.append(rate)
    if len(identifiers)!=len(set(identifiers)):raise ValueError('duplicate source channel IDs')
    if math.fsum(rates)!=total:raise ValueError('source partial rates do not sum exactly to declared fsum total; no renormalization allowed')
    if plan['branch_claims']<len(channels):raise ValueError('whole branch claim family smaller than declared source family')
    expected='integrated-exposure' if varying else 'zero-hazard' if zero else 'rate' if identity=='I22' else 'lifetime'
    if plan.get('primary_estimator')!=expected:raise ValueError('primary estimator differs from preapproved per-ID choice')
    if not zero and not varying:
        cumulative=total*float(plan['horizon_s'])
        if not math.isfinite(cumulative) or cumulative<=0 or not 0<math.exp(-cumulative)<1:
            raise ValueError('declared rate/horizon cannot resolve a nontrivial endpoint survival; explicitly choose a representable physical horizon')
        if expected=='lifetime' and not math.isfinite(1/total):raise ValueError('reference lifetime is not representable')
    return plan


def read_original_parents(path,plan):
    """Validate actual root histories; children/splits cannot become trials."""
    validate_plan(plan)
    with Path(path).open(newline='') as stream:
        reader=csv.DictReader(stream)
        if tuple(reader.fieldnames or ())!=PARENT_COLUMNS:raise ValueError('original-parent CSV schema mismatch')
        rows=list(reader)
    if len(rows)!=plan['original_parents']:raise ValueError('actual independent parent count differs from frozen cohort')
    seen=set();allowed={c['id'] for c in plan['source_channel_rates_s']};horizon=float(plan['horizon_s'])
    weight=float(plan['initial_represented_number'])
    for row in rows:
        key=row['original_parent_id']
        if not key or not key.isascii() or not key.isdigit() or str(int(key))!=key or key in seen:raise ValueError('original parent IDs must be unique canonical nonnegative integers')
        seen.add(key)
        if finite(row['initial_represented_number'],'actual initial represented number',True)!=weight:
            raise ValueError('unequal initial history weights require another reviewed statistical estimator')
        if row['status']=='SURVIVED':
            if row['relative_event_time_s'] or row['channel_id']:raise ValueError('survivor cannot carry a fabricated first-removal event')
        elif row['status']=='REACTED':
            time=finite(row['relative_event_time_s'],'actual first-removal time')
            if time<0 or time>horizon or row['channel_id'] not in allowed:raise ValueError('first-removal event outside physical horizon/source family')
        else:raise ValueError('homogeneous cohort cannot silently censor escape/surface/unknown status')
    return rows


def score_cohort(path,plan):
    """Independent endpoint/branch statistics; NEVER a standalone formal PASS."""
    if plan.get('test_kind')=='trajectory-integrated-hazard':
        raise ValueError('a spatially varying exposure cannot be scored as a constant lifetime')
    rows=read_original_parents(path,plan)
    n=len(rows);survivors=sum(r['status']=='SURVIVED' for r in rows);removed=n-survivors
    horizon=float(plan['horizon_s']);reference=float(plan['reference_total_rate_s'])
    if plan.get('test_kind')=='zero-hazard-invariant':
        # Complete point-Sun shadow/zero-density controls have no finite
        # lifetime to fit. Evaluate their exact no-removal identity instead
        # of inventing a finite error or treating all survivors as a failure.
        return {'schema':'srcMoon-M3-independent-cohort-score-v1',
            'status':'PASS' if removed==0 else 'FAIL',
            'scope':'exact zero-hazard original-parent invariant only',
            'formal_status':'SKIPPED / NOT QUALIFIED BY SCORER ALONE',
            'metrics':{'original_parents':n,'survivors':survivors,'first_removals':removed,
                       'reference_total_rate_s':0.,'primary_estimator':'zero-hazard'},
            'branches':[{'channel_id':c['id'],'count':sum(r['channel_id']==c['id'] for r in rows),
                         'source_probability':None,'status':'PASS' if not any(r['channel_id']==c['id'] for r in rows) else 'FAIL'}
                        for c in plan['source_channel_rates_s']],
            'parent_records_sha256':sha(path),'primary_status':'PASS' if removed==0 else 'FAIL',
            'observational_scoring':False}
    rate=None if survivors==0 else -math.log(survivors/n)/horizon
    lifetime=None if rate is None or rate==0 else 1/rate
    actual=rate if plan['primary_estimator']=='rate' else lifetime
    truth=reference if plan['primary_estimator']=='rate' else 1/reference
    relative=None if actual is None else abs(actual/truth-1)
    metrics={'original_parents':n,'survivors':survivors,'first_removals':removed,
        'endpoint_rate_s':rate,'endpoint_lifetime_s':lifetime,
        'primary_estimator':plan['primary_estimator'],'relative_primary_error':relative,
        'relative_primary_limit':ENSEMBLE_LIMIT[plan['test_id']]}
    primary='PASS' if relative is not None and relative<=metrics['relative_primary_limit'] else 'FAIL'
    branches=[];alpha=.05/plan['branch_claims']
    for source in plan['source_channel_rates_s']:
        p=float(source['rate_s'])/reference;count=sum(row['channel_id']==source['id'] for row in rows)
        if removed==0:
            branches.append({'channel_id':source['id'],'status':'SKIPPED','reason':'no independent removal trials available','count':count});continue
        low,high=clopper_pearson(count,removed,alpha)
        half=(high-low)/(2*p) if p>0 else None
        ok=(count==0) if p==0 else low<=p<=high and half<=.10
        branches.append({'channel_id':source['id'],'status':'PASS' if ok else 'FAIL','count':count,
            'source_probability':p,'frequency':count/removed,'allocated_alpha':alpha,'interval':[low,high],
            'relative_halfwidth':half,'relative_halfwidth_limit':.10 if p>0 else None,
            'endpoint_partial_rate_s':None if rate is None else rate*count/removed})
    exposure=math.fsum(horizon if row['status']=='SURVIVED' else float(row['relative_event_time_s']) for row in rows)
    metrics['event_time_exposure_s']=exposure
    metrics['event_time_rate_diagnostic_s']=removed/exposure if exposure>0 else None
    # Diagnostic estimator is never substituted when the fixed primary fails.
    statuses=[primary]+[r['status'] for r in branches]
    status='FAIL' if 'FAIL' in statuses else 'SKIPPED' if 'SKIPPED' in statuses else 'PASS'
    return {'schema':'srcMoon-M3-independent-cohort-score-v1','status':status,
        'scope':'independent statistics only; native/source/provider/numerical/budget/predecessor binding required',
        'formal_status':'SKIPPED / NOT QUALIFIED BY SCORER ALONE','metrics':metrics,'branches':branches,
        'parent_records_sha256':sha(path),'primary_status':primary,'observational_scoring':False}


def score_refinement(test_id,levels):
    """Three predeclared independent-reference errors, fixed approved budget.

    Each row supplies a positive explicit metric scale and its independently
    derived reference value. No observational/model-output baseline is invented.
    Exact zero at every level is allowed; positive errors must decrease on each
    refinement and the finest must meet the approved numerical bias allocation.
    """
    if test_id not in NUMERICAL_LIMIT or not isinstance(levels,list) or len(levels)<3:raise ValueError('existing ID and at least three refinement levels required')
    count=[];errors=[]
    def exact(value):
        # Native binary64 values retain their exact binary representation;
        # high-precision independent-reference strings retain their original
        # decimal value. Converting both to float BEFORE subtraction could
        # erase an actual numerical error at a large coordinate/offset.
        return Decimal.from_float(value) if isinstance(value,float) else Decimal(str(value))
    for row in levels:
        count.append(integer(row.get('steps'),'predeclared common steps'))
        finite(row.get('scale'),'independent metric scale',True)
        finite(row.get('actual'),'numerical metric');finite(row.get('reference'),'independent reference')
        with localcontext() as context:
            context.prec=80
            errors.append(abs(exact(row['actual'])-exact(row['reference']))/exact(row['scale']))
    if count!=sorted(set(count)):raise ValueError('common refinements must have distinct increasing step counts')
    decreasing=all(b<a or (a==0 and b==0) for a,b in zip(errors,errors[1:]))
    passed=decreasing and errors[-1]<=Decimal(str(NUMERICAL_LIMIT[test_id]))
    numeric=[float(value) for value in errors]
    return {'status':'PASS' if passed else 'FAIL','relative_errors':[v if math.isfinite(v) else None for v in numeric],
        'relative_errors_decimal':[str(value) for value in errors],
        'approved_relative_bias_limit':NUMERICAL_LIMIT[test_id],'decreasing':decreasing,
        'formal_status':'independent numerical metric only; other qualification gates remain'}
