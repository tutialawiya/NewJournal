"""Kinematic cold-pool forward models and conditional source inversions.

Coordinates km, observation times minutes, velocities m/s, orientation degrees
counterclockwise from east. Background is a prescribed translation velocity,
not a fitted or independently validated outflow-propagation velocity.
"""
import numpy as np
from scipy.optimize import least_squares


def arrival_times(stations,source,t_start,speed,wind=(0.,0.),ratio=1.,angle=0.):
    """Expanding advected ellipse; ratio=1 gives an advected circle.

    Semi-axis growth rates: a_dot=c*sqrt(r), b_dot=c/sqrt(r).
    Return first positive arrival; NaN if a front cannot reach a station.
    """
    stations=np.asarray(stations,float);source=np.asarray(source,float)
    if speed<=0 or ratio<1:raise ValueError('Positive speed and axis ratio >=1 required')
    theta=np.deg2rad(angle)
    rotation=np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]])
    displacement=(stations-source)@rotation.T
    velocity=(np.asarray(wind,float)*.06)@rotation.T
    axis=np.array([speed*.06*np.sqrt(ratio),speed*.06/np.sqrt(ratio)])
    d=displacement/axis;v=velocity/axis
    A=np.dot(v,v)-1;B=-2*(d@v);C=np.sum(d*d,axis=1)
    if A< -1e-10:
        elapsed=2*C/(-B+np.sqrt(B*B-4*A*C))
        elapsed=np.where(C<1e-14,0.,elapsed)
    else:
        elapsed=np.full(len(stations),np.nan)
        for i,(b,c) in enumerate(zip(B,C)):
            roots=np.roots([A,b,c]) if abs(A)>1e-10 else np.array([-c/b]) if abs(b)>1e-10 else np.array([])
            positive=[float(r.real) for r in roots if abs(r.imag)<1e-7 and r.real>=0]
            if positive:elapsed[i]=min(positive)
    return float(t_start)+elapsed


def front_points(source,t_start,time,speed,wind=(0.,0.),ratio=1.,angle=0.,points=181):
    elapsed=time-t_start
    if elapsed<0:return np.empty((0,2))
    theta=np.linspace(0,2*np.pi,points)
    local=np.column_stack([np.cos(theta)*speed*.06*np.sqrt(ratio)*elapsed,
                           np.sin(theta)*speed*.06/np.sqrt(ratio)*elapsed])
    a=np.deg2rad(angle);rotation=np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])
    return local@rotation.T+np.asarray(source)+np.asarray(wind)*.06*elapsed


def invert(stations,times,wind=(0.,0.),fixed_speed=None,ratio=1.,angle=0.,free_shape=False,
           extent=100.,starts=16,seed=8,initial=None):
    stations=np.asarray(stations);times=np.asarray(times);rng=np.random.default_rng(seed)
    # Subcritical advection is imposed so the front can reach every azimuth.
    rotation=np.array([[np.cos(np.deg2rad(angle)),np.sin(np.deg2rad(angle))],[-np.sin(np.deg2rad(angle)),np.cos(np.deg2rad(angle))]])
    critical=np.linalg.norm((np.asarray(wind)@rotation.T)/np.array([np.sqrt(ratio),1/np.sqrt(ratio)]))
    minimum=max(1.,1.05*(np.linalg.norm(wind)*np.sqrt(3.) if free_shape else critical))
    if fixed_speed is not None and fixed_speed<=critical:
        return None,[]
    lower=[-extent,-extent,float(times.min()-360)];upper=[extent,extent,float(times.min()-1e-5)]
    if fixed_speed is None:lower.append(minimum);upper.append(30.)
    if free_shape:lower.extend([1.,-89.9]);upper.extend([3.,89.9])
    low=np.array(lower);high=np.array(upper)
    def decode(p):
        c=fixed_speed if fixed_speed is not None else p[3]
        r,theta=(p[-2],p[-1]) if free_shape else (ratio,angle)
        return c,r,theta
    def residual(p):
        c,r,theta=decode(p)
        pred=arrival_times(stations,p[:2],p[2],c,wind,r,theta)
        return np.nan_to_num(pred-times,nan=1e5,posinf=1e5,neginf=-1e5)
    solutions=[]
    candidates=[]
    if initial is not None:candidates.append(np.asarray(initial))
    for i in range(starts):
        p=np.array([rng.uniform(-30,30),rng.uniform(-30,30),times.min()-rng.uniform(10,120)])
        if fixed_speed is None:p=np.r_[p,rng.uniform(minimum,min(14.,29.))]
        if free_shape:p=np.r_[p,rng.uniform(1.1,2.5),rng.uniform(-80,80)]
        candidates.append(p)
    for p in candidates:
        result=least_squares(residual,np.clip(p,low+1e-7,high-1e-7),bounds=(low,high),x_scale='jac',max_nfev=500,ftol=1e-9,xtol=1e-9,gtol=1e-9)
        c,r,theta=decode(result.x)
        values=np.linalg.svd(result.jac,compute_uv=False)
        rank=int(np.sum(values>max(values[0],1e-20)*1e-7))
        solutions.append({'x':float(result.x[0]),'y':float(result.x[1]),'t_start':float(result.x[2]),'speed':float(c),'ratio':float(r),'angle':float(theta),
                          'rmse_minutes':float(np.sqrt(np.mean(result.fun**2))),'jacobian_rank':rank,'parameters':len(result.x),
                          'condition_number':float(values[0]/values[-1]) if values[-1]>1e-12 and rank==len(result.x) else np.inf,
                          'near_domain_boundary':bool(max(abs(result.x[:2]))>.98*extent),
                          'any_optimizer_bound_active':bool(np.any(result.active_mask!=0)),
                          'origin_time_lower_bound_active':bool(result.x[2]<low[2]+1e-4),
                          'parameters_vector':result.x.tolist(),'success':bool(result.success)})
    solutions.sort(key=lambda x:x['rmse_minutes'])
    return solutions[0],solutions
