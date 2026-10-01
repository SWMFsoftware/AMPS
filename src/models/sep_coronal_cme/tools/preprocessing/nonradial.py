"""Steady corotating characteristic/Piola reference producer, with bounded support.

This explicitly solves the Section-6.5 spatial characteristic. An arbitrary
time-dependent twist is not accepted. Callback authorities are supplied by the
offline producer, never searched for in sibling application directories.
"""
from __future__ import annotations
import math
from .core import require, finite

def dot(a,b): return math.fsum(x*y for x,y in zip(a,b))
def norm(a): return math.sqrt(dot(a,a))
def cross(a,b): return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def plus(a,b,factor=1): return [x+factor*y for x,y in zip(a,b)]
def scaled(a,factor): return [factor*x for x in a]
def determinant(m): return dot(m[0],cross(m[1],m[2]))
def inverse(m):
    d=determinant(m); require(d>0,"Piola fold/orientation failure")
    rows=[scaled(cross(m[1],m[2]),1/d),scaled(cross(m[2],m[0]),1/d),scaled(cross(m[0],m[1]),1/d)]
    return [list(row) for row in zip(*rows)]
def matvec(m,v): return [dot(row,v) for row in m]

class SteadyPiolaProvider:
    def __init__(self, start_radius_m, outer_radius_m, radial_step_m, difference_step_m,
                 corotating_velocity, reference_field, authority_fingerprint, boundary_tolerance,
                 physical_boundary_field=None):
        finite([start_radius_m,outer_radius_m,radial_step_m,difference_step_m,boundary_tolerance])
        require(0<start_radius_m<outer_radius_m and radial_step_m>0 and difference_step_m>0 and boundary_tolerance>0,
                "invalid steady winding support/controls")
        require(authority_fingerprint and callable(corotating_velocity) and callable(reference_field), "missing immutable flow/field authority")
        self.start,self.outer,self.step,self.difference=start_radius_m,outer_radius_m,radial_step_m,difference_step_m
        self.velocity,self.reference=corotating_velocity,reference_field
        self.authority,self.tolerance=authority_fingerprint,boundary_tolerance
        self.boundary_field=physical_boundary_field or reference_field

    def boundary_trace(self, direction):
        direction=scaled(direction,1/norm(direction));x=scaled(direction,self.start)
        velocity=self.velocity(x);projection=dot(velocity,direction)
        require(projection>0,"turning boundary characteristic")
        radial_derivative=scaled(velocity,1/projection)
        # Tangential derivatives are the identity sphere; the radial column
        # comes from the flow. Reference and physical one-sided vector fields
        # are distinct authorities when the entering field is nonradial.
        matrix=[[(1.0 if i==j else 0.0)+(radial_derivative[i]-direction[i])*direction[j]
                 for j in range(3)] for i in range(3)]
        actual=scaled(matvec(matrix,self.reference(x)),1/determinant(matrix))
        expected=self.boundary_field(x)
        require(norm(plus(actual,expected,-1))<=self.tolerance*norm(expected),"complete one-sided boundary vector/normal-flux trace mismatch")
        return actual

    def forward(self, reference):
        finite(reference); radius=norm(reference)
        require(self.start<=radius<=self.outer,"winding extrapolation outside support")
        point=scaled(reference,self.start/radius)
        def rate(x):
            velocity=self.velocity(x); finite(velocity)
            projection=dot(velocity,x)/norm(x)
            require(projection>0,"nonpositive radial projection/characteristic turning")
            return scaled(velocity,1/projection)
        n=max(1,math.ceil((radius-self.start)/self.step)); step=(radius-self.start)/n
        for _ in range(n):
            a=rate(point);b=rate(plus(point,a,step/2));c=rate(plus(point,b,step/2));d=rate(plus(point,c,step))
            point=plus(point,[a[i]+2*b[i]+2*c[i]+d[i] for i in range(3)],step/6)
        require(abs(norm(point)-radius)<=self.tolerance*radius,"characteristic radius constraint not converged")
        return point

    def deformation(self, reference):
        radius=norm(reference);require(radius>self.start+self.difference and radius<self.outer-self.difference,
                                       "one-sided boundary deformation needs its separate trace authority")
        columns=[]
        for j in range(3):
            step=[0.,0.,0.];step[j]=self.difference
            a,b=self.forward(plus(reference,step)),self.forward(plus(reference,step,-1))
            columns.append(scaled(plus(a,b,-1),1/(2*self.difference)))
        matrix=[list(row) for row in zip(*columns)]
        require(determinant(matrix)>0,"Piola fold/Jacobian failure")
        return matrix

    def evaluate(self, reference, mass_per_flux, sector, interface_normal=None):
        require(mass_per_flux>0 and sector in {-1,1},"invalid tube mass/sector authority")
        self.boundary_trace(reference)
        point=self.forward(reference);matrix=self.deformation(reference);jacobian=determinant(matrix)
        field=scaled(matvec(matrix,self.reference(reference)),1/jacobian)
        velocity=self.velocity(point)
        require(norm(cross(velocity,field))<=self.tolerance*norm(velocity)*norm(field),
                "steady ideal field-aligned branch/complete boundary field mismatch")
        density=mass_per_flux*norm(field)/norm(velocity)
        normal=None
        if interface_normal is not None:
            inverse_transpose=[list(row) for row in zip(*inverse(matrix))]
            pushed=matvec(inverse_transpose,interface_normal);normal=scaled(pushed,1/norm(pushed))
        return dict(position_m=point,magnetic_field_t=field,velocity_m_per_s=velocity,mass_density_kg_m3=density,
                    mass_per_flux=mass_per_flux,sector=sector,interface_normal=normal,deformation=matrix,jacobian=jacobian)

    def inverse_position(self, physical, seed, tolerance_m, iterations=25):
        # Diffeomorphism follows from a smooth single-valued flow with positive
        # radial projection and an identity boundary. Newton is bounded and any
        # failed inverse refuses provider publication; it does not pick a tube.
        x=list(seed)
        for _ in range(iterations):
            residual=plus(self.forward(x),physical,-1)
            if norm(residual)<=tolerance_m:return x
            x=plus(x,matvec(inverse(self.deformation(x)),residual),-1)
        raise ValueError("winding inverse uniqueness/convergence not established")
