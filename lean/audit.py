"""Independent numerical audit of the attached production solver.
This is not a formal semantics or a validated interval calculation.
Run: python audit.py
"""
from pathlib import Path
import os, sys, json, hashlib, importlib.util, ast, time
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import numpy as np
import mpmath as mp
mp.mp.dps=80
ROOT=Path(__file__).resolve().parent
def exact(x):
    """Interpret each input binary64 exactly before 80-digit evaluation."""
    p, q = float(x).as_integer_ratio()
    return mp.mpf(p) / q


def reference_residual(s):
    z = list(map(exact, s.z))
    n = s.n
    f = [mp.mpf(0)] * (2*n+11)
    f[1] = z[2] - z[1]*exact(s.Hoverd)
    f[2] = z[2] - exact(s.height)*z[3]**2
    f[3] = z[4]*z[3] - 2*mp.pi
    f[4] = z[5]+z[7]-z[4]
    f[5] = z[1]*(z[6]+z[7]-z[4])-z[8]
    f[6] = z[s.Current_criterion+4]-exact(s.Current)*mp.sqrt(z[1])
    f[7] = z[10]+z[n+10]+2*sum(z[11:n+10])
    f[8] = z[10]-z[n+10]-z[2]
    for m in range(n+1):
        X = m*mp.pi/n
        Y = z[10+m]
        psi = u = v = mp.mpf(0)
        for j in range(1,n+1):
            b = z[n+10+j]
            # Direct hyperbolic basis, independent of the exponential implementation.
            S = mp.sinh(j*(Y+z[1]))/mp.cosh(j*z[1])
            C = mp.cosh(j*(Y+z[1]))/mp.cosh(j*z[1])
            psi += b*S*mp.cos(j*X)
            u += j*b*C*mp.cos(j*X)
            v += j*b*S*mp.sin(j*X)
        f[9+m] = psi-z[8]-z[7]*Y
        f[n+10+m] = ((-z[7]+u)**2+v**2)/2+Y-z[9]
    return f[1:]


def reference_point(s, Xf, Yf):
    z = list(map(exact, s.z))
    X, Y, kd = exact(Xf), exact(Yf), z[1]
    u=v=ux=vx=mp.mpf(0)
    for j in range(1,s.n+1):
        b=exact(s.B[j])
        C=mp.cosh(j*(Y+kd))/mp.cosh(j*kd)
        S=mp.sinh(j*(Y+kd))/mp.cosh(j*kd)
        u+=j*b*C*mp.cos(j*X)
        v+=j*b*S*mp.sin(j*X)
        ux-=j*j*b*C*mp.sin(j*X)
        vx+=j*j*b*S*mp.cos(j*X)
    u=(z[5]+u)/mp.sqrt(kd)
    v=v/mp.sqrt(kd)
    ux*=mp.sqrt(kd)
    vx*=mp.sqrt(kd)
    c=z[4]/mp.sqrt(kd)
    return u, v, (u-c)*ux+v*vx


def delta(a,b):
    return abs(float(a)-float(b))


def direct_integrals(s,rho,nphase=256,nvertical=64):
    """Independent integration of the defining physical integrals.
    No call to _point, get_kinematics or _calc_integral_props_cpp.
    Periodic phase quadrature and Gauss-Legendre integration in depth.
    """
    ph=np.arange(nphase)*2*np.pi/nphase
    j=np.arange(1,s.n+1,dtype=float)
    coeff=s.Y[1:s.n+1].copy();coeff[-1]*=0.5
    eta=np.cos(ph[:,None]*j)@coeff/s.k
    nodes,weights=np.polynomial.legendre.leggauss(nvertical)
    dy=(s.d+eta)[:,None]/2
    y=-s.d+dy*(nodes[None,:]+1)
    a=s.k*y[:,:,None]*j; b=s.k*s.d*j
    den=1+np.exp(-2*b)
    C=(np.exp(a)+np.exp(-a-2*b))/den
    S=(np.exp(a)-np.exp(-a-2*b))/den
    amp=j*s.B[1:s.n+1]*np.sqrt(s.g/s.k)
    u=s.EulerianCurrent+np.sum(C*amp*np.cos(ph[:,None,None]*j),axis=2)
    v=np.sum(S*amp*np.sin(ph[:,None,None]*j),axis=2)
    r=s.z[9]*s.g/s.k
    pressure=rho*(r-s.g*y-0.5*((u-s.c)**2+v*v))
    def integral(a): return float(np.mean(np.sum(a*weights[None,:]*dy,axis=1)))
    return {'Impulse':integral(rho*u),
            'KineticEnergy':integral(0.5*rho*(u*u+v*v)),
            'PotentialEnergy':float(0.5*rho*s.g*np.mean(eta*eta)),
            'Sxx':integral(pressure+rho*u*u)-0.5*rho*s.g*s.d*s.d,
            'Power':integral((pressure+rho*s.g*y+0.5*rho*(u*u+v*v))*u)}


def surface_diagnostics(s, points=2048):
    X=np.arange(points)*2*np.pi/points
    j=np.arange(1,s.n+1,dtype=float)
    coef=s.Y[1:s.n+1].copy();coef[-1]*=.5
    eta=np.cos(X[:,None]*j)@coef
    a=eta[:,None]*j;b=s.z[1]*j
    C=(np.exp(a)+np.exp(-a-2*b))/(1+np.exp(-2*b))
    S=(np.exp(a)-np.exp(-a-2*b))/(1+np.exp(-2*b))
    B=s.B[1:s.n+1]
    psi=(S*np.cos(X[:,None]*j))@B
    u=(C*np.cos(X[:,None]*j))@(j*B)
    v=(S*np.sin(X[:,None]*j))@(j*B)
    stream=psi-s.z[8]-s.z[7]*eta
    dynamic=.5*((u-s.z[7])**2+v**2)+eta-s.z[9]
    ux=-(C*np.sin(X[:,None]*j))@(j*j*B)*np.sqrt(s.z[1])
    vx=(S*np.cos(X[:,None]*j))@(j*j*B)*np.sqrt(s.z[1])
    urel=(u-s.z[7])/np.sqrt(s.z[1]);vd=v/np.sqrt(s.z[1])
    ax=s.g*(urel*ux+vd*vx)
    absolute_u=(u+s.z[5])/np.sqrt(s.z[1])*np.sqrt(s.g*s.d)
    sech=2*np.exp(-s.z[1]*j)/(1+np.exp(-2*s.z[1]*j))
    bed_u=(s.z[5]+np.cos(X[:,None]*j)@(j*B*sech))/np.sqrt(s.z[1])*np.sqrt(s.g*s.d)
    return {'points':points,'streamline_max':float(max(abs(stream))),
            'bernoulli_max':float(max(abs(dynamic))),
            'surface_acceleration_dense':float(max(abs(ax))),
            'surface_vertical_velocity_dense':float(max(abs(vd))*np.sqrt(s.g*s.d)),
            'surface_acceleration_reported':float(s.acc_max),
            'surface_vertical_velocity_reported':float(s.w_max),
            'surface_horizontal_crest':float(absolute_u[0]),
            'surface_horizontal_max_abs':float(max(abs(absolute_u))),
            'bed_horizontal_crest':float(bed_u[0]),
            'bed_horizontal_max_abs':float(max(abs(bed_u)))}


def jacobian_check(s):
    z=s.z.copy();A=s._jacobian();num=s.num;maximum=0.
    # Central differences of the direct 80-digit equations, at 30 digits.
    with mp.workdps(35):
        for j in range(1,num+1):
            step=1e-5*max(1.,abs(z[j]))
            s.z[j]=z[j]+step;fp=reference_residual(s)
            s.z[j]=z[j]-step;fm=reference_residual(s)
            actual=exact(z[j]+step)-exact(z[j]-step)
            for i,(p,m) in enumerate(zip(fp,fm)):
                numerical=float((p-m)/actual)
                maximum=max(maximum,abs(numerical-A[i,j-1])/max(1.,abs(A[i,j-1])))
            s.z[j]=z[j]
    s.z[:]=z
    return maximum


def interval_crosscheck(s):
    certificate=s.finite_system_certificate
    f=reference_residual(s)
    z=list(map(exact,s.z));n=s.n
    h,t,d,g,u,rho=map(exact,(s.H_target,s.T_target,s.d,s.g,s.Uc,1025.0))
    f[0]=z[2]-z[1]*h/d
    f[1]=z[2]-h/(g*t*t)*z[3]**2
    f[5]=z[s.Current_criterion+4]-u/mp.sqrt(g*d)*mp.sqrt(z[1])
    result={'all_residuals_enclosed':all(mp.mpf(lo)-mp.mpf('1e-70')<=v<=mp.mpf(hi)+mp.mpf('1e-70') for v,(lo,hi) in zip(f,certificate['residual_intervals']))}
    k=z[1]/d;velocity=mp.sqrt(g/k);scale=rho*g/k**2
    pulse=z[8]+z[1]*z[5]
    ke=(z[4]*pulse-z[5]*(z[7]*z[1]-z[8]))/2
    pe=sum(exact(s.Y[j])**2/4 for j in range(1,n))+exact(s.Y[n])**2/16
    beta=2*z[9]-z[4]**2+2*z[4]*z[5]
    outputs={'L':2*mp.pi/k,'c':z[4]*velocity,'Impulse':rho*velocity/k*pulse,
        'KineticEnergy':scale*ke,'PotentialEnergy':scale*pe,'EnergyDensity':scale*(ke+pe),
        'Sxx':scale*(4*ke-3*pe+beta*z[1]-2*z[5]*pulse),
        'Power':scale*velocity*(z[4]*(3*ke-2*pe)+beta*(pulse+z[4]*z[1])/2-2*z[4]*z[5]*pulse),
        'MeanSquareBedVelocity':g/k*sum((j*z[n+10+j]/mp.cosh(j*z[1]))**2/2 for j in range(1,n+1))}
    for name,value in outputs.items():
        lo,hi=certificate['outputs'][name]['interval']
        result[name+'_enclosed']=mp.mpf(lo)-mp.mpf('1e-70')<=value<=mp.mpf(hi)+mp.mpf('1e-70')
    return result


def run_case(mod, case, backend):
    H,T,d,U,criterion=case
    mod.NUMBA_AVAILABLE=backend=='numba'
    s=mod.FentonStreamFunction(H,T,d,U)
    s.Current_criterion=criterion;s.solve()
    result={'case':case,'backend':backend,'converged':bool(s.converged),'message':s.last_error,
            'accuracy_verified':bool(s.accuracy_verified),'status':s.solution_status,
            'accuracy_metrics':s.accuracy_metrics,
            'finite_system_verified':s.finite_system_verified,
            'finite_system_certificate':s.finite_system_certificate}
    if not s.converged and s.solution_status != 'resolution_failure':return result
    result["interval_crosscheck"] = interval_crosscheck(s)
    assert all(result["interval_crosscheck"].values()), result["interval_crosscheck"]
    f=reference_residual(s);rhs=np.zeros_like(s.z);s._eqns(rhs)
    point_error=0.
    for X in (0.,np.pi/3,np.pi):
        for Y in (-s.z[1],-s.z[1]/2,s._surface_keta(X)):
            point_error=max(point_error,max(delta(a,b) for a,b in zip(s._point(X,Y),reference_point(s,X,Y))))
    rho,g,d,c,U,Q=map(exact,(mod.RHO,s.g,s.d,s.c,s.EulerianCurrent,s.VolumeFluxQ))
    k=exact(s.k);ys=list(map(exact,s.Y))
    pe=rho*g/k**2*(sum(y*y/4 for y in ys[1:s.n])+ys[s.n]**2/16)
    impulse=rho*(c*d-Q);ke=(c*impulse-Q*rho*U)/2
    bed2=g/k*sum((j*exact(s.B[j])/mp.cosh(j*exact(s.z[1])))**2/2 for j in range(1,s.n+1))
    absolute=bed2+U**2
    power=c*(3*ke-2*pe)+absolute*(impulse+rho*c*d)/2-2*c*U*impulse
    stress=4*ke-3*pe+rho*absolute*d-2*U*impulse
    momentum=stress-2*c*impulse+rho*d*(c*c+g*d/2)
    reference={'Impulse':impulse,'KineticEnergy':ke,'PotentialEnergy':pe,
               'MeanSquareBedVelocity':bed2,'MeanSquareBedAbsoluteVelocity':absolute,
               'Power':power,'Sxx':stress,'EnergyDensity':ke+pe,'MomentumFlux':momentum}
    errors={name:delta(getattr(s,name),v) for name,v in reference.items()}
    checks={name:errors[name]<=1e-7+1e-8*abs(float(v)) for name,v in reference.items()}
    # Single-phase trapezoidal momentum flux has a different discretisation.
    checks['MomentumFlux']=errors['MomentumFlux']<=1e-6+1e-5*abs(float(momentum))
    quad=direct_integrals(s,mod.RHO)
    quad_errors={key:delta(getattr(s,key),v) for key,v in quad.items()}
    node_error=max(abs(s._surface_keta(m*np.pi/s.n)-s.z[10+m]) for m in range(s.n+1))
    result.update(N=s.n,condition=s.jacobian_condition,fraction=s.continuation_fraction,
                  residual_80digits=float(max(map(abs,f))),
                  residual_code_difference=max(delta(a,b) for a,b in zip(rhs[1:],f)),
                  point_difference=point_error,node_difference=float(node_error),
                  outputs={key:float(getattr(s,key)) for key in list(reference)+['L','c','k','Cg','EulerianCurrent','StokesCurrent','VolumeFluxQ','acc_max','w_max','ExcursionBed','u_bed','u_surf']},
                  reference={key:float(v) for key,v in reference.items()},errors=errors,checks=checks,
                  quadrature=quad,quadrature_errors=quad_errors,
                  surface=surface_diagnostics(s),z=s.z.tolist(),Y=s.Y[:s.n+1].tolist())
    result['checks'].update(residual=result['residual_80digits']<=1e-10,
                           residual_implementation=result['residual_code_difference']<=1e-10,
                           kinematics=point_error<=1e-10,surface_nodes=node_error<=1e-10,
                           conditioning=result['condition']<=1e11,target_height=result['fraction']==1.)
    return result


def guard_checks(mod):
    mod.NUMBA_AVAILABLE=False;results={}
    for i,case in enumerate([(0,9,10,0),(1,-1,10,0),(1,9,float('inf'),0),(1,9,10,float('nan'))]):
        s=mod.FentonStreamFunction(*case);s.solve()
        results['invalid_input_'+str(i)]=not s.converged and bool(s.last_error)
    s=mod.FentonStreamFunction(1,9,10,0);s.Current_criterion=3;s.solve()
    results['invalid_current_selector']=not s.converged
    s=mod.FentonStreamFunction(1,9,10,0);s.nstep=1;s.number=3;s._newton=lambda _:0.
    s.solve();results['stagnation_rejected']=not s.converged and s.residual_max>s.residual_tolerance_final
    s=mod.FentonStreamFunction(1,9,10,0);original=s._compute_Y_and_B
    def perturb_after_step():original();s.z[9]+=.1
    s._compute_Y_and_B=perturb_after_step;s.nstep=1;s.solve()
    results['postprocess_residual_guard']=not s.converged
    s=mod.FentonStreamFunction(1,9,10,0)
    s.z[1]=1.;s.z[4]=1.;s.z[7]=1.;s.z[9]=.5;s.Y[s.n]=.2
    s._calc_integral_props_cpp()
    results['terminal_mode_energy_weight']=abs(s.PotentialEnergy-mod.RHO*s.g*s.d**2*.2**2/16)<1e-10
    return results


def refinement_checks(mod, report):
    mod.NUMBA_AVAILABLE=False
    result=[]
    for case in [(2,3,5,0),(3,4,5,1)]:
        s=mod.FentonStreamFunction(*case);s.solve()
        coarse=direct_integrals(s,mod.RHO,256,64)
        fine=direct_integrals(s,mod.RHO,512,96)
        finer=direct_integrals(s,mod.RHO,1024,128)
        result.append({'case':case,'quadrature_256_64':coarse,
                       'quadrature_512_96':fine,'quadrature_1024_128':finer,
                       'fine_to_finer_difference':{k:abs(fine[k]-finer[k]) for k in fine},
                       'surface_8192':surface_diagnostics(s,8192)})
    report['refinement_checks']=result


def main():
    sources=json.loads((ROOT/'sources.json').read_text())
    for key in ('python','readme'):
        path=ROOT/sources[key]['filename']
        assert hashlib.sha256(path.read_bytes()).hexdigest()==sources[key]['sha256'],f'{key}: source hash mismatch'
    path=ROOT/sources['python']['filename']
    spec=importlib.util.spec_from_file_location('fenton_production_audit',path)
    mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
    available=mod.NUMBA_AVAILABLE
    cases=[(1,9,10,0,1),(1,9,10,.5,1),(1,9,10,-.5,1),(1,9,5,0,1),(.2,4,60,0,1),
           (2,9,10,0,1),(1,9,10,.2,2),(3,9,5,0,1),(3,9,5,1,1),
           *[(h,t,5,u,1) for h in (1,2) for t in (3,4) for u in (0,1)],
           (3,4,5,1,1),(1,3,20,-.5,1),(2,4,5,.2,2),(1,3,5,.5,2)]
    report={'scope':'Empirical numerical tests of the acceptance contract, including explicit accuracy rejections; not a universal program proof.',
            'sources':sources,'numpy':np.__version__,'mpmath':mp.__version__,
            'precision_digits':80,'cases':[]}
    (ROOT/'results').mkdir(exist_ok=True)
    for backend in ['numpy']+(['numba'] if available else []):
        for case in cases:
            r=run_case(mod,case,backend);report['cases'].append(r)
            failed=[k for k,v in r.get('checks',{}).items() if not v]
            print(backend,case,'converged',r['converged'],'failed',failed,flush=True)
            (ROOT/'results/numerical.json').write_text(json.dumps(report,indent=2,allow_nan=False,default=lambda x: x.item()))
    report['guards']=guard_checks(mod)
    mod.NUMBA_AVAILABLE=False
    report['jacobian_checks']=[]
    for case in [(1.,9.,10.,.5,1),(1.,3.,5.,.5,2)]:
        s=mod.FentonStreamFunction(*case[:4],fourier_order=12);s.Current_criterion=case[4];s.solve()
        err=jacobian_check(s)
        report['jacobian_checks'].append({'case':case,'N':s.n,'max_scaled_error':err,'pass':err<1e-7})
        print('Jacobian',case,err,flush=True)
    report['rejections']=[]
    for case in [(3,3,5,0),(3,3,5,1),(3,4,5,0)]:
        s=mod.FentonStreamFunction(*case);s.solve()
        report['rejections'].append({'case':case,'rejected':not s.converged,'message':s.last_error})
        print('Target rejection',case,not s.converged,flush=True)
    (ROOT/'results/numerical.json').write_text(json.dumps(report,indent=2,allow_nan=False,default=lambda x: x.item()))
    refinement_checks(mod,report)
    (ROOT/'results/numerical.json').write_text(json.dumps(report,indent=2,allow_nan=False,default=lambda x: x.item()))
    print('Completed audit; inspect the report, including discretisation diagnostics.',flush=True)
    budgets_pass=all((c['converged'] and c['accuracy_verified'] and c['finite_system_verified'] and all(c.get('checks',{}).values()))
                    or (not c['converged'] and not c['accuracy_verified']
                        and c['status']=='resolution_failure' and not all(c.get('checks',{}).values()))
                    for c in report['cases'])
    extrema_pass=all(c['outputs']['u_bed'] >= c['surface']['bed_horizontal_max_abs']-1e-8
                     and c['outputs']['u_surf'] >= c['surface']['surface_horizontal_max_abs']-1e-8
                     and c['outputs']['u_bed'] >= 0 and c['outputs']['u_surf'] >= 0
                     for c in report['cases'] if c['converged'])
    all_guards=all(report['guards'].values()) and all(j['pass'] for j in report['jacobian_checks'])
    if not (budgets_pass and extrema_pass and all_guards):
        print('AUDIT STATUS: acceptance-contract check failed.',flush=True)
        return 2
    print('AUDIT STATUS: all finite-matrix acceptance and rejection checks passed; no universal program certificate is claimed.',flush=True)
    return 0

if __name__=='__main__':sys.exit(main())


