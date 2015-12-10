#q,kappa,dkappadr,dPsi_T/dr
from toroidal_flux_q import toroidal_flux_derivative,func_polyfit,func_polyfit_derivative,func_q
#psi(r),dpsiNdr
from r_to_psi import func_psi,psi_pedestal_width,drdPsiN_at_psiN_one
import numpy

def drdpsiN_of_r(a,simul):
    q95=simul.inputs.Miller_q
    kappa90=simul.inputs.Miller_kappa
    RBar=simul.RBar
    BBar=simul.BBar
    #rescale q_profile to fit whatever q95 we have in the simulation
    scale_q=q95/3.5
    scale_kappa=kappa90/1.58
    #actual profile has q95~3.5
    
    kappay=scale_kappa*numpy.array([1.23,1.37,1.58])
    kappax=numpy.array([0,0.5,0.9])
    deg=2
    #in what follows, rho=r/a
    kappa_rho=func_polyfit(kappax,kappay,deg) #kappa(rho)
    dkappadrho=func_polyfit_derivative(kappax,kappay,deg) #dkappadrho(rho) 

    #generate the integrated B_T as function of minor radius of eliptical cross-section
    BT_r=lambda x: 1/simul.inputs.RHat(0.0) #should be B_T(r), here taken as constant and at theta=0
    dBTdr=lambda x: 0 #dBTdr(r)
    #calculate radial derivative of toroidal flux Psi_T (assumes small theta dependence in B)
    dPsiTdr=toroidal_flux_derivative(a,kappa_rho,dkappadrho,BT_r,dBTdr)

    #calculate q profile
    q_rho=func_q(1.0,3.5*scale_q,0.6,1.6*scale_q) #q(rho)
    q_r=lambda r:q_rho(r/a) #q(rho=r/a)

    #calculate psi(r)
    psi_r=func_psi(q_r,dPsiTdr) #psi(r), will be in the same units as BBar and RBar
    psiAHat=psi_r(a)/(BBar*RBar**2) #translates to/from unitless psi to psi_N

    simul.inputs.changevar("physicsParameters","psiAHat",psiAHat)
    simul.inputs.read(simul.input_filename)
    #print  drdPsiN_at_psiN_one(q_r,dPsiTdr,psi_r,a)
    return drdPsiN_at_psiN_one(q_r,dPsiTdr,psi_r,a)[0]
