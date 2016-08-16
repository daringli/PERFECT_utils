#this class generates n_i, n_z and Phi for a given n_e, T_e, T_i and T_z
#so that that the orderings in PERFECT are satisfied

from perfect_simulation import perfect_simulation,normalized_perfect_simulation #to get info about simulation
import h5py #to write profiles
import numpy #matrix things, etc
import scipy.optimize #solve horrible system
from bezier_transition import bezier_transition, derivative_bezier_transition #for smooth transition between 2 functions
from perfectInputFile import perfectInput
from perfectPsiAHatProfileFile import create_psiAHat_of_Npsi

#for main
import matplotlib.pyplot as plt
from sys import argv

def generate_nonuniform_grid_bezier(inputFilename,psi1,psi2,slope12,smoothness):
    #Calculates psi(psiN), where psiN is uniform and psi is nonuniform

    #psi1, psi2: psi points between which to use denser grid
    #slope12: slope between corresponding psiN1, psiN2 points
    
    
    inputs = perfectInput(inputFilename)
    Npsi = inputs.Npsi
    psiAHat = inputs.psiAHat #psiAHat for a would-be uniform map
    psiN = inputs.psi #array of uniform grid points
    psiNMin = inputs.psiMin
    psiNMax = inputs.psiMax

    #slope01, slope23: slope outside tighter grid section. Assumed equal
    slope01 = slope12*(psiAHat*(psiNMax - psiNMin) - (psi2-psi1))/(slope12*(psiNMax - psiNMin) - (psi2-psi1))
    slope23=slope01

    # psiN points corresponding to psi1, psi2
    psiN1 = psiNMin + (psi1 - psiAHat*psiNMin)/slope01
    psiN2 = psiNMax - (psiAHat*psiNMax - psi2)/slope23
    
    #create linear functions to transition between
    inner = (lambda psiN: (psiAHat*psiNMin + slope01*(psiN-psiNMin)))
    middle = (lambda psiN: inner(psiN1) + slope12*(psiN-psiN1))
    outer = (lambda psiN: middle(psiN2) + slope23*(psiN-psiN2))

    derivativeInner = (lambda psiN: slope01)
    derivativeMiddle = (lambda psiN: slope12)
    derivativeOuter = (lambda psiN: slope23)
    
    
    #generate end points for bezier curve from smoothness parameter
    #smoothness of 1 gives endpoints of transitions in the middle of the pedestal
    #so that it just has time to transition to middle linear function
    distance = smoothness*(psiN2-psiN1)/2
    psiNList = [psiN1,psiN2]
    pair1=[distance, distance]
    pair2=[distance, distance]
    pairList=[pair1,pair2]
    splineList=[inner,middle,outer]
    derivativeSplineList = [derivativeInner,derivativeMiddle,derivativeOuter]

    #calculate nonuniform grid points array
    #psi=bezier_transition(splineList,psiNList,pairList,psiN) 
    psi,dpsidpsiN=derivative_bezier_transition(splineList,derivativeSplineList,psiNList,pairList,psiN) 

    create_psiAHat_of_Npsi("psiAHat.h5",Npsi,dpsidpsiN(psiN),psi(psiN))
    
    return (psi,dpsidpsiN)

def generate_nonuniform_grid_arctanh(inputFilename,slope,psi0):
    #Calculates psi(psiN), where psiN is uniform and psi is nonuniform
    # psi = h(psiN)
    # h = Atanh(P(psiN)), where P = a2 s**psiN + a1*psiN + a0
    # polynomial determined to make h(psiMin) = psiMin, h(psiMax) = psiMax
    # and have the slope slope at psi0
    inputs = perfectInput(inputFilename)
    Npsi = inputs.Npsi
    psiAHat = inputs.psiAHat #psiAHat for a would-be uniform map
    psiN = inputs.psi #array of uniform grid points
    psiNMin = inputs.psiMin
    psiNMax = inputs.psiMax
    t1 = numpy.tanh(psiNMin)
    t2 = numpy.tanh(psiNMax)
    k = slope
    A = numpy.array([[psiNMin**3,psiNMin**2,psiNMin,1],[psiNMax**3,psiNMax**2,psiNMax,1],[psi0**3,psi0**2,psi0,1],[3*psi0**2,2*psi0,1,0]])
    b = numpy.array([t1,t2,0,k])
    coef=numpy.linalg.solve(A,b)

    h = (lambda psiN: psiAHat*numpy.arctanh(coef[0]*psiN**3 + coef[1]*psiN**2 + coef[2]*psiN+coef[3]))
    dhdpsiN = (lambda psiN: psiAHat*(3*coef[0]*psiN**2 + 2*coef[1]*psiN + coef[2])/(1+(coef[0]*psiN**3 + coef[1]*psiN**2 + coef[2]*psiN+coef[3])**2))

    psi=[h(pN) for pN in psiN] 
    dpsidpsiN= [dhdpsiN(pN) for pN in psiN]

    create_psiAHat_of_Npsi("psiAHat.h5",Npsi,dpsidpsiN,psi)
    
    return (psi,dpsidpsiN)

def generate_nonuniform_grid_Int_arctan(inputFilename,a,b,c,psiN1,psiN2):
    #Calculates psi(psiN), where psiN is uniform and psi is nonuniform
    # psi = h(psiN)
    # dh/dpsiN = c - b[1 + (2/pi)atan((x-s_1)/a)] - [1 + (2/pi)atan((x-s_2)/a)]
    # parameters determined to set transition speed and pedestal value, 0<b<c, a<0
    # psiN0: start of simulation region
    # psiN1: start of pedestal region
    # psiN2: end of pedestal region
    # psiN3: end of simulation region
    # s0,s1,s2,s3: corresponding uniform grid points
    # psi0,psi3: corresponding psiN*psiAHat points
    pi = numpy.pi
    log=numpy.log
    arctan=numpy.arctan

    inputs = perfectInput(inputFilename)
    psiAHat=inputs.psiAHat
    psiN=inputs.psi

    psiN0=psiN[0]
    psiN3=psiN[-1]
    psi0=psiAHat*psiN0
    psi3=psiAHat*psiN3
    dpsiN=psiN3-psiN0

    a=a*dpsiN
    c=c*psiAHat
    b=b*psiAHat
    
    # set start of simulation region at psi0 by convention
    s0=psiN0
    #need to solve for s1, s2 before we specify h
    # we do this numerically since h is hellish

    #solve for pedestal width in s
    dpsi_ped=psiAHat*(psiN2-psiN1)
    f1 = (lambda ds_ped: c*ds_ped + (2*a*b/pi)*(0.5*log((ds_ped/a)**2 + 1) - (ds_ped/a)*arctan(ds_ped/a)))
    ds_ped = scipy.optimize.fsolve((lambda ds_ped: f1(ds_ped) - dpsi_ped),psiAHat*dpsi_ped)
    #solve for s1, expressing s2=s1 + ds_ped
    dpsi_10 = psiAHat*(psiN1 - psiN0)
    f2 = (lambda s1: c*(s1-s0) -  (a*b/pi)*(
        0.5*log(((s0-s1)/a)**2 + 1) - 0.5*log(((s0-(s1 + ds_ped))/a)**2 + 1)
        + ((s0-(s1 + ds_ped))/a)*arctan((s0-(s1 + ds_ped))/a) - ((s0-s1)/a)*arctan((s0-s1)/a)
        - (ds_ped/a)*arctan(ds_ped/a) + 0.5*log((ds_ped/a)**2+1)
                                          )
    )
    s1 = scipy.optimize.fsolve((lambda s1: f2(s1) - dpsi_10),psiAHat*dpsi_10)
    s2 = s1 + ds_ped
    
    #we can now get functions for h and dh/ds
    # smooth step functions and characteristic function of pedestal region
    H1 = (lambda s : 0.5 + (1/pi)*arctan((s-s1)/a))
    H2 = (lambda s : 0.5 + (1/pi)*arctan((s-s2)/a))
    chi =  (lambda s : H1(s) - H2(s))
    dhds = (lambda s: c - b*chi(s))
    h = (lambda s: psi0 + c*(s-s0) - (a*b/pi)*(
        0.5*log(((s0-s1)/a)**2 + 1) - 0.5*log(((s0-s2)/a)**2 + 1)
        + ((s0-s2)/a)*arctan((s0-s2)/a) - ((s0-s1)/a)*arctan((s0-s1)/a)
        + ((s-s1)/a)*arctan((s-s1)/a) - ((s-s2)/a)*arctan((s-s2)/a)
        + 0.5*log(((s-s2)/a)**2 + 1) - 0.5*log(((s-s1)/a)**2 + 1)
    )
    )

    #solve for s3
    s3 = scipy.optimize.fsolve((lambda x:h(x) - psi3),psiN3)

    # since we'd rather have s3 at psiN3, we rescale s linearly
    rescale=(psiN3-psiN0)/(s3-s0)
    s1=s0+(s1-s0)*rescale
    s2=s0+(s2-s0)*rescale
    s3=s0+(s3-s0)*rescale
    b=b/rescale
    c=c/rescale
    a=a*rescale
    print "s1: " + str(s1)
    print "s2: " + str(s2)
    print "s3: " + str(s3)

    #uniform grid goes between s0 and s3
    Npsi = inputs.Npsi
    s = numpy.linspace(s0,s3,Npsi)
    psi=h(s)
    dpsids= dhds(s)

    create_psiAHat_of_Npsi("psiAHat.h5",Npsi,dpsids,psi)

    #plt.plot(s,H1(s))
    #plt.plot(s,-H2(s))
    #plt.show()
    
    #return (s,psi,dpsids,s1,s2)
    return h,dhds

##################################
#    For testing purposes        #
##################################


    
if __name__ == "__main__":
    inputFilename = argv[1]
    bezier=False
    arctanh=False
    Int_arctan=True

    #get psiAHat to translate our old psiN coordinate to a psi coordinate
    inputs = perfectInput(inputFilename)
    
    if bezier:
        psiAHat = inputs.psiAHat

        #pedestal start and stop in psi
        psi1 = 0.94927395957025573*psiAHat
        psi2 = 0.97463697978512787*psiAHat

        psiN = inputs.psi #we use the "old psi grid" as our new uniform psiN grid

        #psiN(psi) slope in area with denser grid


        f, axarr = plt.subplots(2,2)

        slopeOverPsiAHatList = [1, 0.5, 0.25]
        smoothnessList = [0.5,0.5, 0.5]
        
        for slopeOverPsiAHat,smoothness in zip(slopeOverPsiAHatList,smoothnessList): 
            slope = psiAHat*slopeOverPsiAHat

            psi,dpsidpsiN=generate_nonuniform_grid_bezier(inputFilename,psi1,psi2,slope,smoothness)

            axarr[0,0].plot(psiN,psi,label="slope " + str(slopeOverPsiAHat))
            axarr[1,0].plot(psiN,dpsidpsiN,label="slope " + str(slopeOverPsiAHat))


        slopeOverPsiAHatList = [0.25, 0.25, 0.25]
        smoothnessList = [1,0.5, 0.25]

        for slopeOverPsiAHat,smoothness in zip(slopeOverPsiAHatList,smoothnessList): 
            slope = psiAHat*slopeOverPsiAHat

            psi,dpsidpsiN=generate_nonuniform_grid_bezier(inputFilename,psi1,psi2,slope,smoothness)

            axarr[0,1].plot(psiN,psi)

            axarr[1,1].plot(psiN,dpsidpsiN)

            #plt.legend(bbox_to_anchor=(0.03, 0.95), loc=2, borderaxespad=0.)
        axarr[0,0].set_title("Varying slope")
        axarr[0,0].set_ylabel(r"$\psi/(\bar{R}^2\bar{B})$")
        axarr[1,0].set_ylabel(r"$d/d\psi_N (\psi/(\bar{R}^2\bar{B}))$")

        axarr[0,1].set_title("Varying smoothness")
        axarr[1,0].set_xlabel(r"$\psi_N$")
        axarr[1,1].set_xlabel(r"$\psi_N$")



        plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
        plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
        plt.show()
    elif arctanh:
        f, axarr = plt.subplots(2,1)
        psiAHat = inputs.psiAHat
        psi1 = 0.94927395957025573*psiAHat
        psi2 = 0.97463697978512787*psiAHat
        psi0 = (psi2 + psi1)/(2*psiAHat)
        print "psi0: " + str(psi0)
        psiN = inputs.psi
        
        slopeOverPsiAHatList = [0.01,0.1,0.2]
        for slopeOverPsiAHat in slopeOverPsiAHatList: 
            slope = slopeOverPsiAHat/psiAHat

            psi,dpsidpsiN=generate_nonuniform_grid_arctanh(inputFilename,slope,psi0)

            axarr[0].plot(psiN,psi)

            axarr[1].plot(psiN,dpsidpsiN)

            #plt.legend(bbox_to_anchor=(0.03, 0.95), loc=2, borderaxespad=0.)
        axarr[0].set_title("Varying slope")
        axarr[0].set_ylabel(r"$\psi/(\bar{R}^2\bar{B})$")
        axarr[1].set_ylabel(r"$d/d\psi_N (\psi/(\bar{R}^2\bar{B}))$")

        
        plt.setp([axarr[0].get_xticklabels()], visible=False)
        #plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
        plt.show()
    elif Int_arctan:
        f, axarr = plt.subplots(2,2)
        psiAHat = inputs.psiAHat
        psiN1 = 0.94927395957025573
        psiN2 = 0.97463697978512787
        psiN0= inputs.psi[0]
        psiN3= inputs.psi[-1]
        
        dpsiN=(psiN3-psiN0)
        print "psiAHat " + str(psiAHat)
        
        colorlist=['b','g','r','k']
        alist=[0.5,0.1,0.025,0.01]
        for a,col in zip(alist,colorlist):
            a=a
            b=0.5
            c=1

            s,psi,dpsidpsiN,s1,s2=generate_nonuniform_grid_Int_arctan(inputFilename,a,b,c,psiN1,psiN2)
            #print s
            #print psi
            axarr[0,0].plot(s,psi/psiAHat,color=col)
            axarr[1,0].plot(s,dpsidpsiN/psiAHat,color=col)
            
            axarr[0,0].axvline(x=s1,color=col,linestyle=':')
            axarr[0,0].axhline(y=psiN1,color=col,linestyle=':')
            axarr[1,0].axvline(x=s1,color=col,linestyle=':')

            axarr[0,0].axvline(x=s2,color=col,linestyle=':')
            axarr[0,0].axhline(y=psiN2,color=col,linestyle=':')
            axarr[1,0].axvline(x=s2,color=col,linestyle=':')

        blist=[0.125,0.25,0.5,0.9]
        for b,col in zip(blist,colorlist): 
            a=0.025
            b=b
            c=1

            s,psi,dpsidpsiN,s1,s2=generate_nonuniform_grid_Int_arctan(inputFilename,a,b,c,psiN1,psiN2)
            axarr[0,1].plot(s,psi/psiAHat,color=col)
            axarr[1,1].plot(s,dpsidpsiN/psiAHat,color=col)

            axarr[0,1].plot(s,psi/psiAHat,color=col)
            axarr[1,1].plot(s,dpsidpsiN/psiAHat,color=col)
            
            axarr[0,1].axvline(x=s1,color=col,linestyle=':')
            axarr[0,1].axhline(y=psiN1,color=col,linestyle=':')
            axarr[1,1].axvline(x=s1,color=col,linestyle=':')

            axarr[0,1].axvline(x=s2,color=col,linestyle=':')
            axarr[0,1].axhline(y=psiN2,color=col,linestyle=':')
            axarr[1,1].axvline(x=s2,color=col,linestyle=':')

           

            #plt.legend(bbox_to_anchor=(0.03, 0.95), loc=2, borderaxespad=0.)
        axarr[0,0].set_title("Varying transition")
        axarr[0,1].set_title("Varying density")

        axarr[0,0].set_ylabel(r"$h/\hat{\psi}_A$")
        axarr[1,0].set_ylabel(r"$d_s h/\hat{\psi}_A$")

        axarr[1,0].set_xlabel(r"$s$")
        axarr[1,1].set_xlabel(r"$s$")
        axarr[1,0].set_ylim([0,1.5])
        axarr[1,1].set_ylim([0,1.5])

        axarr[0,0].set_xlim([s[0],s[-1]])
        axarr[0,1].set_xlim([s[0],s[-1]])
        axarr[1,0].set_xlim([s[0],s[-1]])
        axarr[1,1].set_xlim([s[0],s[-1]])

        axarr[0,0].set_ylim([s[0],s[-1]])
        axarr[0,1].set_ylim([s[0],s[-1]])



        plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
        plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
        plt.show()
