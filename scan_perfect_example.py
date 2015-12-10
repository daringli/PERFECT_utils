from perfect_simulation import perfect_simulation,normalized_perfect_simulation #to get info about simulation
from generate_compatible_profiles import generate_compatible_profiles
from generate_psi_of_r import drdpsiN_of_r


def simul_manipulation(simul):
    #increase impurity Z and mass (roughly)
    z_index=1
    simul.inputs.charges[z_index]=simul.inputs.charges[z_index]+1
    simul.inputs.masses[z_index]=simul.inputs.masses[z_index]+1
    #str_of_array=' '.join(map(str,l))
    
    simul.inputs.changevar("speciesParameters","charges",' '.join(map(str,simul.inputs.charges)))
    simul.inputs.changevar("speciesParameters","masses",' '.join(map(str,simul.inputs.masses)))


#simulation initialization
dir="test"
input_filename="input.namelist"
norm_filename="norms.namelist"
output_filename="perfectOutput.h5"



Tped=0.9
mI=0
zI=1
eI=2
a=0.7 #minor radius in meters
xwidth=0.03 #pedestal width in r

Tped_e=0.9
dTCoredx_e=-0.1*Tped_e/xwidth
dTpeddx_e=-Tped_e/xwidth
dTSOLdx_e=-0.05*Tped_e/xwidth

Tped_d=0.9
dTCoredx_d=-0.1*Tped_d/xwidth
dTpeddx_d=-0.1*Tped_d/xwidth
dTSOLdx_d=-0.1*Tped_d/xwidth

Tped_N=0.9
dTCoredx_N=-0.1*Tped_N/xwidth
dTpeddx_N=-0.1*Tped_N/xwidth
dTSOLdx_N=-0.1*Tped_N/xwidth

nped_d=0.4
dnCoredx_d=-0.6
dnpeddx_d=-nped_d/xwidth
dnSOLdx_d=-0.6

cmd="perfect -ksp_type gmres -ksp_gmres_restart 100 -ksp_max_it 1000"

num_simuls=5
simuls=[0]*num_simuls
simul_dir_names=[dir+str(x+1) for x in range(num_simuls)]


original_simul=normalized_perfect_simulation(input_filename,norm_filename,None)
original_simul.species=["d","N","e"]
dxdpsiN=drdpsiN_of_r(a,original_simul) #generate dr/dpsi_N. Returns a function

for i in range(0,num_simuls):
    if i==0:
        simuls[0]=original_simul.copy_simulation_norm_to_dir(simul_dir_names[0])
    else:
        simuls[i]=simuls[i-1].copy_simulation_norm_to_dir(simul_dir_names[i])
        simul_manipulation(simuls[i])
        #raise pedestal gradient for the ion density (makes them less negative)
        dnpeddx_d=dnpeddx_d+0.2
    generate_compatible_profiles(simuls[i],mI=mI,zI=zI,eI=eI,xwidth=xwidth,a=a,dxdpsiN=dxdpsiN,Tped_e=Tped_e,dTCoredx_e=dTCoredx_e,dTpeddx_e=dTpeddx_e,dTSOLdx_e=dTSOLdx_e,Tped_d=Tped_d,dTCoredx_d=dTCoredx_d,dTpeddx_d=dTpeddx_d,dTSOLdx_d=dTSOLdx_d,Tped_N=Tped_N,dTCoredx_N=dTCoredx_N,dTpeddx_N=dTpeddx_N,dTSOLdx_N=dTSOLdx_N,nped_d=nped_d,dnCoredx_d=dnCoredx_d,dnpeddx_d=dnpeddx_d,dnSOLdx_d=dnSOLdx_d)
    simuls[i].always_run_simulation(cmd)
