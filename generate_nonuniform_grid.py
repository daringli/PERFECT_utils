#this class generates n_i, n_z and Phi for a given n_e, T_e, T_i and T_z
#so that that the orderings in PERFECT are satisfied

from perfect_simulation import perfect_simulation,normalized_perfect_simulation #to get info about simulation
import h5py #to write profiles
import numpy #matrix things, etc
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
    psi=bezier_transition(splineList,psiNList,pairList,psiN) 
    dpsidpsiN=derivative_bezier_transition(splineList,derivativeSplineList,psiNList,pairList,psiN) 

    create_psiAHat_of_Npsi("psiAHat.h5",Npsi,dpsidpsiN,psi)
    
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


##################################
#    For testing purposes        #
##################################


    
if __name__ == "__main__":
    inputFilename = argv[1]
    bezier=False

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
    else:
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
