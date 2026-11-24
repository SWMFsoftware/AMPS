#!/usr/bin/env python3
"""Independent first-hit timing audit of actual native Na cohorts.

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
import re
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
            # For a moving case the verifier supplies the independent
            # inverse-square reference at this *actual* event position.
            point_rate=reference_rate(event) if callable(reference_rate) else reference_rate
            for key in ('source_rate_event_s','accepted_rate_s'):
                if key not in event:raise ValueError('actual native point-source rate metric unavailable: '+key)
                value=event[key]
                actual=Decimal(value.numerator)/Decimal(value.denominator)
                relative=None if point_rate==0 else abs(actual-point_rate)/point_rate
                passed=(actual==0) if relative is None else relative<=I21_INVERSE_SQUARE_RELATIVE_LIMIT
                checks.append({'original_parent_id':parent,'quantity':key,'actual_s':str(actual),
                    'independent_s':str(point_rate),'relative_error':None if relative is None else str(relative),
                    'relative_limit':str(I21_INVERSE_SQUARE_RELATIVE_LIMIT),'status':'PASS' if passed else 'FAIL'})
    return {'status':'FAIL' if any(r['status']=='FAIL' for r in checks) else 'PASS','checks':checks,
        'events_evaluated':len(events),'scope':'existing fixed I21 inverse-square criterion; empty set is only a bookkeeping identity'}


def lunar_segment_geometry(source, sun, radius_m):
    """Project the origin onto the *finite* source-to-Sun segment in LSO metres.

    Fraction parses the printed SI coordinates exactly. Keeping the decision
    in squared distance avoids introducing rounding in the grazing case. The
    square root is used only for the human-readable distance in the artifact.
    """
    if len(source) != 3 or len(sun) != 3:
        raise ValueError('source and Sun need exactly three LSO coordinates')
    point = [Fraction(value) for value in source]
    ray = [Fraction(value)-initial for value,initial in zip(sun,point)]
    radius = Fraction(radius_m)
    norm2 = sum(value*value for value in ray)
    if radius <= 0 or norm2 == 0:
        raise ValueError('positive lunar radius and nonzero finite solar segment required')
    projection = -sum(initial*delta for initial,delta in zip(point,ray))/norm2
    position = min(Fraction(1),max(Fraction(0),projection))
    distance2 = sum((initial+position*delta)**2 for initial,delta in zip(point,ray))
    return {'intersection': 'interior' if distance2 < radius*radius else
            'exterior' if distance2 > radius*radius else 'tangent',
            'radius_m': str(radius), 'closest_segment_fraction': str(position),
            'minimum_distance_squared_m2': str(distance2),
            'minimum_distance_m': math.sqrt(float(distance2))}


def fixed_moon_shadow(plan,launch,receipt):
    """Certify a lunar-only dark case independently of the native rate and flag.

    An authenticated fixed-position launch enforces both zero velocity and the
    zero-force binding. With the frozen environment one segment is sufficient
    for the entire physical horizon. Other geometries need a different oracle.
    """
    from m2_initial_state import assignments
    if receipt['case_kind'] != 'homogeneous-fixed-position':
        raise q.MissingPrerequisite('fixed lunar shadow requires a stationary source')
    # prepare_cohort freezes the generated native input followed by the actual
    # application root and its recursive includes. The owner authenticates
    # the complete graph before invoking this check.
    native_input=Path(launch['input']['path'])
    documents=launch['input_documents']
    if len(documents)<2 or Path(documents[0]['path']).resolve()!=native_input.resolve():
        raise ValueError('missing native or application root in frozen input graph')
    application=Path(documents[1]['path'])
    settings=assignments(application)
    earth=re.findall(r'^\s*earth_occultation\s*=\s*(off|on)\s*$',native_input.read_text(),re.M)
    if (settings.get('surface_geometry',{}).get('value')!='sphere'
            or settings.get('simulation_environment_mode',{}).get('value')!='frozen'
            or earth!=['off']):
        raise q.MissingPrerequisite('fixed lunar-only oracle requires a frozen sphere and Earth occultation off')
    descriptor=plan['native_adapter']['producer_manifest']
    path=Path(descriptor['path'])
    if q.sha(path)!=descriptor['sha256']:
        raise ValueError('frozen native producer manifest changed')
    producer=q.strict_object(path)
    key='src/general/constants.PlanetaryData.h'
    header=Path(producer['source_root'])/key
    if not header.is_file():
        raise q.MissingPrerequisite('production lunar-radius source unavailable')
    if q.sha(header)!=producer['production_source_hashes'][key]:
        raise ValueError('lunar-radius source differs from the compiled producer')
    radii=re.findall(r'^\s*#define\s+_MOON__RADIUS_\s+([0-9.eE+-]+)\s*$',header.read_text(),re.M)
    if len(radii)!=1:
        raise ValueError('expected one explicit lunar radius in metres in the compiled planetary constants')
    reference=reference_receipt(receipt)
    geometry=lunar_segment_geometry(reference['source_position_lso_m'],reference['sun_position_lso_m'],radii[0])
    if geometry['intersection']=='tangent':
        raise q.MissingPrerequisite('exact grazing contact has no strict fixed lunar-shadow classification')
    geometry.update(status='PASS' if geometry['intersection']=='interior' and reference['point_sun_visible']==0 else 'FAIL',
        production_point_sun_visible=reference['point_sun_visible'],
        lunar_radius_source={'path':str(header),'sha256':q.sha(header)},
        application_input={'path':str(application),'sha256':q.sha(application)},
        frame='Moon-centred LSO',geometry_basis_epoch_et_s=reference['geometry_basis_epoch_et_s'])
    return geometry


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


def frozen_reference(plan):
    """Authenticate the optional reference snapshot produced by plan preparation.

    Older explicit plans carry the independently calculated rate without a
    reference snapshot. New prepared plans also freeze the geometry basis and
    every D04 dependency, so equal rates cannot hide a changed source point.
    """
    binding=plan.get('reference_binding')
    if binding is None:return None,{}
    receipt=binding['receipt'];path=Path(receipt['path'])
    if not path.is_file():raise q.MissingPrerequisite('frozen geometry receipt unavailable: '+str(path))
    digest=q.sha(path)
    if digest!=receipt['sha256']:raise ValueError('frozen geometry receipt identity changed')
    reference=reference_receipt(q.strict_object(path))
    identities={}
    source_rate(reference,identities)
    if identities!=binding['data_sha256']:
        raise ValueError('D04 identities differ from the frozen cohort reference')
    for record in binding.get('spice_files',[]):
        kernel=Path(record['path'])
        if not kernel.is_file():raise q.MissingPrerequisite('frozen SPICE input unavailable: '+str(kernel))
        digest=q.sha(kernel)
        if digest!=record['sha256']:raise ValueError('frozen SPICE input identity changed: '+str(kernel))
        identities[str(kernel)]=digest
    identities[str(path)]=receipt['sha256']
    return reference,identities


class StraightLineExposure:
    """Independent frozen-geometry oracle with an analytic transverse bound.

    The declared particle moves along x. Its distance squared from the Sun
    is d(t)^2+b^2 with d(t)>0 and fixed transverse b. Integrating the axial
    problem gives an upper bound on H; dividing by 1+b^2/min(d)^2 gives a
    lower bound. The same two functions bracket the optical-depth root.
    """
    def __init__(self,plan,receipt,initial_rate):
        wanted=plan.get('trajectory_reference',{})
        actual=receipt.get('trajectory_reference',{})
        if wanted.get('oracle')!='illuminated-zero-force-lso-x-axis-v1':
            raise ValueError('explicit frozen straight-line oracle identity required')
        if actual.get('zero_force') is not True or actual.get('terrain_faces')!=0:
            raise ValueError('numerical oracle requires actual zero-force motion and no terrain faces')
        self.start=self.vector(receipt['source_position_lso_m'])
        self.sun=self.vector(receipt['sun_position_lso_m'])
        self.velocity=self.vector(actual['source_velocity_lso_ms'])
        if self.velocity!=self.vector(wanted['source_velocity_lso_ms']) or not self.velocity[0] or any(self.velocity[1:]):
            raise ValueError('actual and declared nonzero x-only velocity must match')
        self.horizon=Decimal.from_float(float(plan['horizon_s']))
        with localcontext() as ctx:
            ctx.prec=80
            end=self.start[0]+self.velocity[0]*self.horizon
            self.d0=self.sun[0]-self.start[0]
            minimum=self.sun[0]-max(self.start[0],end)
            if minimum<=0 or initial_rate<=0:
                raise ValueError('straight-line reference requires a positive distance and rate throughout')
            # Independent sufficient proof for each reference ray segment:
            # all x coordinates of its endpoints exceed the sphere's xmax.
            # The whole segment is then on the same strictly separated side.
            ray_min=min(self.start[0],end,self.sun[0])
            spheres=actual['opaque_spheres']
            if not isinstance(spheres,list) or not spheres:
                raise ValueError('native Moon opacity is required in the linked spherical fixture')
            for sphere in spheres:
                centre=self.vector(sphere['centre_lso_m'])
                radius=Decimal.from_float(q.finite(sphere['radius_m'],'sphere radius',True))
                if ray_min<=centre[0]+radius:
                    raise q.MissingPrerequisite('straight-line oracle cannot prove full finite-ray visibility')
            transverse=sum((a-b)**2 for a,b in zip(self.sun[1:],self.start[1:]))
            self.relative_projection_bound=transverse/minimum**2
            self.strength=initial_rate*(self.d0**2+transverse)
            self.end_depth_upper=self.axial_depth(self.horizon)
            self.end_depth_lower=self.end_depth_upper/(1+self.relative_projection_bound)

    @staticmethod
    def vector(values):
        if not isinstance(values,list) or len(values)!=3:
            raise ValueError('explicit SI three-vector required')
        return tuple(Decimal.from_float(q.finite(v,'SI component')) for v in values)

    def point_rate(self,event):
        positions=(event[f'parent_event_x{i}'] for i in range(3))
        squared=sum((sun-Decimal(x.numerator)/Decimal(x.denominator))**2
                    for sun,x in zip(self.sun,positions))
        if squared<=0:
            raise ValueError('native event hit the singular point source')
        return self.strength/squared

    def axial_depth(self,t):
        return self.strength*t/(self.d0*(self.d0-self.velocity[0]*t))

    def first_hit(self,threshold):
        if threshold>self.end_depth_upper:
            return None
        if threshold>self.end_depth_lower:
            raise q.MissingPrerequisite('transverse analytical envelope cannot classify this parent')
        def root(depth):
            return depth*self.d0**2/(self.strength+depth*self.d0*self.velocity[0])
        return root(threshold),root(threshold*(1+self.relative_projection_bound))


def homogeneous_timing_breakdown(depth,reference_rate,native,display_text,precision='binary64-log'):
    """Signed time differences in seconds, independent of the acceptance test.

    The historical source clock stores binary64 -log(U). Its archived
    receipts retain the native event time and instantaneous rate, but do not
    serialize the original libm result. Calculate the correctly rounded clock
    independently from Decimal's log, and attribute any discrepancy between
    libm and that reference to the *combined* transport remainder. In a
    compensated-extended-log case the receipt does not separately retain the
    evaluated libm result. Use the unrounded Decimal depth for the accepted
    rate comparison, and leave the extended log contribution inside the
    combined native clock/transport remainder rather than estimating it.
    In a homogeneous case a diagnostic first-hit root at the actual instantaneous
    rate does not require a model of the integration/clock update history.
    """
    def decimal_fraction(value):
        return Decimal(value.numerator)/Decimal(value.denominator)

    with localcontext() as context:
        context.prec=80
        if not {'accepted_rate_s','relative_event_time_s'}<=native.keys():
            raise ValueError('homogeneous first hit requires the exact native time and accepted rate')
        rate=decimal_fraction(native['accepted_rate_s'])
        if depth<=0 or reference_rate<=0 or rate<=0:
            raise ValueError('a homogeneous first hit requires positive optical depth and rates')
        if precision not in ('binary64-log','compensated-extended-log'):
            raise ValueError('unsupported optical clock precision in native receipt')
        clock_depth=Decimal.from_float(float(depth)) if precision=='binary64-log' else depth
        independent=depth/reference_rate
        clock_at_reference_rate=clock_depth/reference_rate
        root_at_native_rate=clock_depth/rate
        exact_native=decimal_fraction(native['relative_event_time_s'])
        projected=Decimal.from_float(float(native['relative_event_time_s']))
        displayed=Decimal(display_text)
        # Preserve a missing component in extended mode. Reporting zero
        # would imply exact knowledge of the native libm rounding error.
        quantified_clock = ({'source_clock_rounding_s':clock_at_reference_rate-independent}
                            if precision=='binary64-log' else {})
        return {
            **quantified_clock,
            'accepted_rate_effect_s':root_at_native_rate-clock_at_reference_rate,
            'native_clock_and_transport_remainder_s':exact_native-root_at_native_rate,
            'parent_binary64_projection_s':projected-exact_native,
            'parent_csv_text_rendering_s':displayed-projected,
            'native_minus_reference_s':exact_native-independent,
            'csv_text_minus_reference_s':displayed-independent,
        }


def verify(runs,plan):
    """Audit native-output directories for frozen increasing common steps.

    `runs` supplies (output-directory, mpi-ranks, parsed exact native events)
    for every actual launch. This is
    a per-cohort numerical proof, not global I21 qualification. All membership,
    source and horizon checks precede timing scores; corrupt CSV is ERROR.
    """
    q.validate_plan(plan);reference_basis,identities=frozen_reference(plan);comparisons=[];levels={1:[],8:[]}
    baseline=None
    timing_breakdown=[]
    varying=plan.get('test_kind')=='trajectory-integrated-hazard'
    expected_kind='trajectory-integrated-hazard' if varying else 'homogeneous-fixed-position'
    for directory,ranks,events in runs:
        if type(ranks) is not int or ranks not in levels:raise ValueError('explicit1/8MPI layouts required')
        directory=Path(directory);receipt_path=directory/'input_receipt.json'
        receipt=q.strict_object(receipt_path);identities[str(receipt_path)]=q.sha(receipt_path)
        if receipt['case_kind']!=expected_kind:raise ValueError('actual transport case does not match frozen plan')
        for actual,expected,label in ((receipt['original_parents'],plan['original_parents'],'parent count'),
                (receipt['seed'],plan['seed'],'seed'),(receipt['cohort_id'],plan['cohort_id'],'semantic cohort'),
                (float(receipt['physical_horizon_s']),float(plan['horizon_s']),'common horizon'),
                (float(receipt['initial_represented_number']),float(plan['initial_represented_number']),'represented number')):
            if actual!=expected:raise ValueError('actual native input contradicts frozen '+label)
        reference=reference_receipt(receipt);rate,channel=source_rate(reference,identities)
        if reference_basis is not None and reference!=reference_basis:
            raise ValueError('actual geometry, endpoint or zero-force binding differs from the frozen reference')
        precision=receipt.get('optical_clock_precision','binary64-log')
        if precision not in ('binary64-log','compensated-extended-log'):
            raise ValueError('unsupported native source-clock precision')
        if 'optical_clock_precision' in plan and precision!=plan['optical_clock_precision']:
            raise ValueError('native source-clock precision differs from frozen plan')
        geometry=json.dumps([reference,receipt.get('trajectory_reference'),precision],sort_keys=True)
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
        motion=StraightLineExposure(plan,receipt,rate) if varying else None
        maximum=Decimal(0);horizon=Decimal.from_float(float(plan['horizon_s']))
        maxima={}
        with localcontext() as context:
            context.prec=80
            for row in rows:
                parent=int(row['original_parent_id']);u=Decimal.from_float(source_uniform(plan['seed'],plan['cohort_id'],parent))
                depth=-u.ln()
                root_interval=motion.first_hit(depth) if motion else None
                event=(sum(root_interval)/2 if root_interval is not None else None) if motion else (None if rate==0 else depth/rate)
                expected='REACTED' if event is not None and event<=horizon else 'SURVIVED'
                if row['status']!=expected:
                    raise ArithmeticError('native first-removal membership differs from independent clock: parent'+str(parent))
                actual_time=row['relative_event_time_s']
                position_error=None
                if expected=='REACTED':
                    native=events[row['original_parent_id']]
                    exact_time=native['relative_event_time_s']
                    actual_time=Decimal(exact_time.numerator)/Decimal(exact_time.denominator)
                    if motion is None:
                        components=homogeneous_timing_breakdown(depth,rate,native,row['relative_event_time_s'],precision)
                        # A full cohort can contain many thousands of events.
                        # Keep one actual signed example at the largest
                        # magnitude for each cause and level. Every parent is
                        # still scored against its independent physical root.
                        for name,value in components.items():
                            if name not in maxima or abs(value)>abs(maxima[name][1]):
                                maxima[name]=(parent,value)
                if motion is not None and expected=='REACTED':
                    positions=[native[f'parent_event_x{i}'] for i in range(3)]
                    actual_position=[Decimal(x.numerator)/Decimal(x.denominator) for x in positions]
                    position_error=max(abs(x-x0-v*event) for x,x0,v in zip(actual_position,motion.start,motion.velocity))
                    if position_error>Decimal(str(q.NUMERICAL_LIMIT['I21']))*abs(motion.velocity[0])*horizon:
                        raise ArithmeticError('native straight-line event position exceeds I21 numerical allocation')
                error=Decimal(0) if event is None or expected=='SURVIVED' else abs(Decimal(actual_time)-event)/horizon
                maximum=max(maximum,error)
                comparisons.append({'mpi_ranks':ranks,'steps':receipt['physical_steps'],'original_parent_id':parent,
                    'expected_status':expected,'actual_status':row['status'],'reference_event_s':None if event is None else str(event),
                    'actual_event_s':str(actual_time),'normalized_error':str(error),
                    'reference_root_width_s':None if root_interval is None else str(root_interval[1]-root_interval[0]),
                    'position_error_m':None if position_error is None else str(position_error)})
        levels[ranks].append({'steps':receipt['physical_steps'],'actual':str(maximum),'reference':'0','scale':'1'})
        if motion is None:
            timing_breakdown.append({'mpi_ranks':ranks,'steps':receipt['physical_steps'],
                'events_evaluated':len(events),'optical_clock_precision':precision,'component_extrema':{
                    name:{'original_parent_id':parent,'signed_seconds':str(value),
                          'absolute_over_horizon':str(abs(value)/horizon)}
                    for name,(parent,value) in maxima.items()}})
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
    return {'schema':'srcMoon-M3-independent-'+('trajectory' if varying else 'homogeneous')+'-numerics-v1',
        'status':status,'scores':scores,
        'comparisons':comparisons,'input_sha256':identities,'finite_uniform_resolution':finite_resolution,
        'homogeneous_timing_breakdown':timing_breakdown,
        'scope':'actual per-cohort source-clock timing; no observational or full I21 qualification'}
