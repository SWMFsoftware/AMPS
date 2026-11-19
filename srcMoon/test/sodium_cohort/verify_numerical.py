#!/usr/bin/env python3
"""Independent first-hit timing audit of actual homogeneous native cohorts.

No production RNG, rate, chemistry or exposure implementation is called. Real
D04 native data and QA remain authoritative; semantic HMAC and Decimal logs
are independent numerical references, never generated observations.
"""
from __future__ import annotations
import csv
from decimal import Decimal,localcontext,InvalidOperation
import hashlib
import hmac
import json
from pathlib import Path
import struct
import sys
import math
from fractions import Fraction

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import m3_qualification as q

D04=Path('/data/vtenishe/moon_validation_data/photochemistry')
# Frozen I21 acceptance.json: inverse-square rate relative error <=0.1%.
# This is separate from the approved0.2% total numerical-bias allocation and
# every caller-supplied arithmetic/ledger bound. It is not configurable here.
I21_INVERSE_SQUARE_RELATIVE_LIMIT=Decimal('0.001')


def adjacent_binary64(value,upward):
    """Exact adjacent IEEE binary64, including signed zero/subnormals.

    Python3.8 has no math.nextafter. Bit stepping is the sole implementation,
    rather than a numeric epsilon or a platform-dependent fallback. Native
    transport already requires IEEE binary64; unsupported host floats are an
    explicit missing prerequisite. Finite endpoints may return an infinity,
    which the owning epoch-representation check must refuse explicitly.
    """
    if sys.float_info.radix!=2 or sys.float_info.mant_dig!=53:
        raise q.MissingPrerequisite('independent clock proof requires IEEE binary64 Python floats')
    if not isinstance(value,float) or not math.isfinite(value) or type(upward) is not bool:
        raise ValueError('adjacent binary64 requires a finite float and explicit direction')
    if value==0:
        bits=1 if upward else (1<<63)|1
    else:
        bits=struct.unpack('>Q',struct.pack('>d',value))[0]
        bits+=1 if (value>0)==upward else -1
    return struct.unpack('>d',struct.pack('>Q',bits))[0]


def event_time_links(events,parents,origin_text):
    """Reconcile first-parent times with actual native physical event epochs.

    The CSV stores binary64 first-event times, while exact hex diagnostics
    retain the production relative long-double value. Their required binary64
    projection is exact. Absolute epochs may differ by the established TWO
    adjacent binary64 additions, independently implemented by exact bit steps;
    this representational envelope is not a configurable physical time error.
    Geometry basis epoch is deliberately not substituted for transport origin.
    """
    origin=q.finite(origin_text,'actual common transport origin')
    checks=[];records={row['original_parent_id']:row for row in parents if row['status']=='REACTED'}
    if set(records)!=set(events):raise ValueError('first-parent event set contradicts exact native event metrics')
    for parent,event in events.items():
        if 'relative_event_time_s' not in event or 'event_epoch_et_s' not in event:
            raise ValueError('actual relative/absolute event chronology diagnostics unavailable')
        relative=event['relative_event_time_s'];epoch=event['event_epoch_et_s']
        if relative<0:raise ValueError('native relative first-event time is negative')
        displayed=q.finite(records[parent]['relative_event_time_s'],'reported first-event time')
        projected=float(Fraction.from_float(origin)+relative)
        if not math.isfinite(projected):raise ValueError('native absolute time projection overflows')
        low=high=projected
        for _ in range(2):
            low=adjacent_binary64(low,False);high=adjacent_binary64(high,True)
        if not math.isfinite(low) or not math.isfinite(high):raise ValueError('absolute epoch neighbour envelope exceeds finite binary64')
        checks.append({'original_parent_id':parent,'criterion':'CSV-exact-native-relative-time-projection',
            'status':'PASS' if displayed==float(relative) else 'FAIL'})
        checks.append({'original_parent_id':parent,'criterion':'native-physical-epoch-from-actual-transport-origin',
            'status':'PASS' if Fraction.from_float(low)<=epoch<=Fraction.from_float(high) else 'FAIL'})
    return {'status':'FAIL' if any(r['status']=='FAIL' for r in checks) else 'PASS','checks':checks,
        'events_evaluated':len(events),'scope':'exact clock-receipt reconciliation; independent first-hit timing remains separate'}


def inverse_square_checks(events,reference_rate):
    """Evaluate both actual source-law and accepted fixed-point event rates.

    Input is the exact hexadecimal parser's Fraction mapping. A native flag,
    a precomputed residual or a loose user ledger tolerance is not a reference.
    Whole-cohort event completeness is verified by the owning audit separately.
    """
    checks=[]
    with localcontext() as context:
        context.prec=80
        for parent,event in events.items():
            for key in ('source_rate_event_s','accepted_rate_s'):
                if key not in event:raise ValueError('actual native point-source rate metric unavailable: '+key)
                value=event[key]
                actual=Decimal(value.numerator)/Decimal(value.denominator)
                relative=None if reference_rate==0 else abs(actual-reference_rate)/reference_rate
                passed=(actual==0) if relative is None else relative<=I21_INVERSE_SQUARE_RELATIVE_LIMIT
                checks.append({'original_parent_id':parent,'quantity':key,'actual_s':str(actual),
                    'independent_s':str(reference_rate),'relative_error':None if relative is None else str(relative),
                    'relative_limit':str(I21_INVERSE_SQUARE_RELATIVE_LIMIT),'status':'PASS' if passed else 'FAIL'})
    return {'status':'FAIL' if any(r['status']=='FAIL' for r in checks) else 'PASS','checks':checks,
        'events_evaluated':len(events),'scope':'existing fixed I21 inverse-square criterion; empty set is only a bookkeeping identity'}


def reference_receipt(receipt):
    """Read the native schema's explicit nested geometry, without guessing.

    Top-level geometry also supports trace diagnostics. If duplicated there,
    it must agree exactly with the numerical reference block; a missing point
    visibility cannot be inferred from a solar cosine or selected event rate.
    """
    if receipt.get('schema')!='srcMoon-native-sodium-input-receipt-v1':raise ValueError('unsupported actual native input receipt schema')
    fields=('solar_endpoint','total_acceleration_macro','point_sun_visible','sun_position_lso_m',
        'source_position_lso_m','source_reference_distance_m','source_reference_distance_identity','geometry_basis_epoch_et_s')
    reference=receipt['homogeneous_reference']
    if not isinstance(reference,dict) or set(reference)!=set(fields):raise ValueError('complete explicit native homogeneous reference required')
    for key in fields:
        top='geometry_basis_et_s' if key=='geometry_basis_epoch_et_s' else key
        if top in receipt and receipt[top]!=reference[key]:raise ValueError('contradictory duplicated native geometry: '+key)
    return reference


def source_uniform(seed,cohort,parent):
    """Independently reproduce the approved source-clock purpose1/draw0.

    Standard Python hmac implements RFC2104, independent of the production
    C++ SHA service. A source parent has an empty genealogy. Conversion of the
    first52 digest bits to the midpoint binary64 grid is exact.
    """
    u64=lambda n:struct.pack('>Q',n)
    key=u64(seed);label=cohort.encode('utf-8')
    message=(b'AMPS-MOON-M3-HMAC-SHA256-COUNTER-v1'+u64(len(label))+label+
             u64(parent)+u64(0)+u64(1)+u64(0))
    bits=int.from_bytes(hmac.new(key,message,hashlib.sha256).digest()[:7],'big')>>4
    return (bits+.5)*2.**-52


def source_rate(reference,identities):
    """Read complete real Na endpoint, retaining original binary64 parse.

    Hash every declared processed output before interpreting source records.
    Missing data/QA is SKIPPED; changed bytes or conflicting native records are
    ERROR. Native Sun coordinates are a geometry-service input; I15 owns their
    physical qualification. No independent ephemeris is invented here.
    """
    for name in ('README.md','provenance.json','qa_report.json'):
        path=D04/name
        if not path.is_file():raise q.MissingPrerequisite('D04 metadata unavailable: '+str(path))
        identities[str(path)]=q.sha(path)
    (D04/'README.md').read_text()
    if q.strict_object(D04/'qa_report.json').get('status')!='PASS':
        raise q.MissingPrerequisite('D04 extraction QA is not PASS')
    for item in q.strict_object(D04/'provenance.json')['outputs']:
        path=D04/item['path']
        if not path.is_file():raise q.MissingPrerequisite('D04 source unavailable: '+str(path))
        if q.sha(path)!=item['sha256']:raise ValueError('D04 native source identity mismatch: '+str(path))
        identities[str(path)]=item['sha256']
    endpoint=reference['solar_endpoint']
    if endpoint not in ('quiet','active'):raise ValueError('explicit separate native solar endpoint required')
    table=D04/'processed/phidrates/photo_rates_1au.csv'
    with table.open(newline='') as stream:
        family=[r for r in csv.DictReader(stream) if r['species']=='Na']
    if len(family)!=2 or {r['solar_activity'] for r in family}!={'0.0','1.0'}:
        raise ValueError('complete native Na endpoint family required')
    rows=[r for r in family if r['solar_activity']==('0.0' if endpoint=='quiet' else '1.0')]
    row=rows[0]
    if row['optical_depth']!='0.0' or row['distance_au']!='1.0':raise ValueError('D04 native normalization changed')
    macro=reference['total_acceleration_macro']
    if not isinstance(macro,str):raise ValueError('actual acceleration macro expansion required')
    macro=''.join(macro.split())
    zero_bindings={prefix+'TotalParticleAcceleration_default(a,s,p,x,v,n)'+suffix
        for prefix in ('','PIC::Mover::') for suffix in ('',';')}
    if macro not in zero_bindings:
        raise ValueError('homogeneous proof requires the actual existing zero-acceleration production binding')
    visible=reference['point_sun_visible']
    if type(visible) is not int or visible not in (0,1):raise ValueError('actual point-Sun0/1 visibility required')
    sun=reference['sun_position_lso_m'];point=reference['source_position_lso_m']
    if len(sun)!=3 or len(point)!=3:raise ValueError('actual LSO Sun/source vectors required')
    def D(x):
        if isinstance(x,bool):raise ValueError('native geometry must be explicit SI numbers, not booleans')
        try:value=Decimal(str(x))
        except (InvalidOperation,ValueError,TypeError) as error:raise ValueError('invalid native SI geometry number') from error
        if not value.is_finite():raise ValueError('native geometry must be finite SI values')
        return value
    with localcontext() as context:
        context.prec=80
        squared=sum((D(a)-D(b))**2 for a,b in zip(sun,point))
        unit=D(reference['source_reference_distance_m'])
        if squared<=0 or unit<=0 or not reference.get('source_reference_distance_identity'):
            raise ValueError('source distance and its normalization/provenance required')
        # The maintained source registry parses each archived scalar to double.
        # Preserve that actual representation; no table value is tuned/rounded
        # to agree with the native event timing under evaluation.
        rate=Decimal.from_float(float(row['rate_1au_s-1']))*unit*unit/squared*visible
    return rate,row['reaction_id']


def verify(runs,plan):
    """Audit native-output directories for frozen increasing common steps.

    `runs` is a list of (output-directory, mpi-ranks) actual launches. This is
    a per-cohort numerical proof, not global I21 qualification. All membership,
    source and horizon checks precede timing scores; corrupt CSV is ERROR.
    """
    q.validate_plan(plan);identities={};comparisons=[];levels={1:[],8:[]}
    baseline=None
    for directory,ranks in runs:
        if type(ranks) is not int or ranks not in levels:raise ValueError('explicit1/8MPI layouts required')
        directory=Path(directory);receipt_path=directory/'input_receipt.json'
        receipt=q.strict_object(receipt_path);identities[str(receipt_path)]=q.sha(receipt_path)
        if receipt['case_kind']!='homogeneous-fixed-position':raise ValueError('variable exposure cannot enter constant-hazard cohort scoring')
        for actual,expected,label in ((receipt['original_parents'],plan['original_parents'],'parent count'),
                (receipt['seed'],plan['seed'],'seed'),(receipt['cohort_id'],plan['cohort_id'],'semantic cohort'),
                (float(receipt['physical_horizon_s']),float(plan['horizon_s']),'common horizon'),
                (float(receipt['initial_represented_number']),float(plan['initial_represented_number']),'represented number')):
            if actual!=expected:raise ValueError('actual native input contradicts frozen '+label)
        reference=reference_receipt(receipt);rate,channel=source_rate(reference,identities)
        geometry=json.dumps(reference,sort_keys=True)
        if baseline is None:baseline=geometry
        elif baseline!=geometry:raise ValueError('geometry/source/force convention differs between matched refinements')
        if len(plan['source_channel_rates_s'])!=1 or plan['source_channel_rates_s'][0]['id']!=channel:
            raise ValueError('frozen cohort omits or changes complete native Na family')
        declared=Decimal.from_float(float(plan['reference_total_rate_s']))
        # Only the known binary64 parse/geometry arithmetic projection is
        # allowed here. This is a fixed representation envelope, not the0.2%
        # numerical gate or an adjustable physical/source uncertainty.
        if abs(declared-rate)>abs(rate)*Decimal(64)*Decimal(2)**-52:
            raise ValueError('declared source rate contradicts independently scaled D04 rate')
        path=directory/'parent_records.csv';rows=q.read_original_parents(path,plan)
        identities[str(path)]=q.sha(path)
        if {int(r['original_parent_id']) for r in rows}!=set(range(plan['original_parents'])):
            raise ValueError('actual cohort must retain every frozen global original-parent ID')
        maximum=Decimal(0);horizon=Decimal.from_float(float(plan['horizon_s']))
        with localcontext() as context:
            context.prec=80
            for row in rows:
                parent=int(row['original_parent_id']);u=Decimal.from_float(source_uniform(plan['seed'],plan['cohort_id'],parent))
                event=None if rate==0 else -u.ln()/rate
                expected='REACTED' if event is not None and event<=horizon else 'SURVIVED'
                if row['status']!=expected:
                    raise ArithmeticError('native first-removal membership differs from independent clock: parent'+str(parent))
                error=Decimal(0) if event is None or expected=='SURVIVED' else abs(Decimal(row['relative_event_time_s'])-event)/horizon
                maximum=max(maximum,error)
                comparisons.append({'mpi_ranks':ranks,'steps':receipt['physical_steps'],'original_parent_id':parent,
                    'expected_status':expected,'actual_status':row['status'],'reference_event_s':None if event is None else str(event),
                    'actual_event_s':row['relative_event_time_s'],'normalized_error':str(error)})
        levels[ranks].append({'steps':receipt['physical_steps'],'actual':str(maximum),'reference':'0','scale':'1'})
    scores={}
    for ranks,rows in levels.items():
        rows.sort(key=lambda r:r['steps'])
        scores[str(ranks)]=q.score_refinement('I21',rows)
    if [r['steps'] for r in levels[1]]!=[r['steps'] for r in levels[8]]:
        raise ValueError('matched1/8MPI refinements must use the same frozen steps')
    status='PASS' if all(r['status']=='PASS' for r in scores.values()) else 'FAIL'
    # The installed single-step Na topology has at most ONE accepted parent
    # removal: source clock1, channel1, recoil2 and daughter clock1. Count all
    # launches conservatively even though matched layouts reuse semantic keys.
    # This retained support/tail diagnostic is not a substitute for numerical
    # refinement, PRF/statistical checks or a physical model-error bound.
    finite_resolution=q.finite_grid_tail_bound(5*plan['original_parents']*len(runs))
    return {'schema':'srcMoon-M3-independent-homogeneous-numerics-v1','status':status,'scores':scores,
        'comparisons':comparisons,'input_sha256':identities,'finite_uniform_resolution':finite_resolution,
        'scope':'actual per-cohort source-clock timing; no observational or full I21 qualification'}
