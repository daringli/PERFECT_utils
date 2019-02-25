from __future__ import division

import numpy as np
from perfectGeometryFile import perfectGeometry

from diff_matrix import diff_matrix

from perfectInputFile import perfectInput

from scipy.interpolate import interp1d, interp2d

import h5py
import f90nml # You may need to install this module if you do not have it already, e.g. "pip install --user f90nml"


def read_file(filename):
    with open(filename,'r') as f:
        l = f.readline()
        ls = l.split()
        array = np.array([float(x) for x in ls])
    return array
        
class EXPRO_geometry(object):
        
    def read_psi(self):
        self.psi = read_file(self.dirname + "/" + self.psi_filename)
        self.psia = self.psi[-1]
        self.psiN = self.psi/self.psia
        self.Npsi = len(self.psi)

    def read_theta(self):
        self.theta = read_file(self.dirname + "/" + self.theta_filename)
        self.Ntheta = len(self.theta)

    def read_I(self):
        self.I = read_file(self.dirname + "/" + self.I_filename)
        
    def read_B(self):
        array = read_file(self.dirname + "/" + self.B_filename)
        self.B=array.reshape(self.shape,order=self.order)

    def read_R(self):
        array = read_file(self.dirname + "/" + self.R_filename)
        self.R=array.reshape(self.shape,order=self.order)
        
    def read_J(self):
        array = 1/read_file(self.dirname + "/" + self.J_filename)
        self.J=array.reshape(self.shape,order=self.order)
        
    def __init__(self,dirname):
        self.order='F'
        
        self.dirname = dirname
        self.B_filename = "B"
        self.R_filename = "bigR"
        self.J_filename = "jacob"
        self.I_filename = "I"
        self.psi_filename = "psi"
        self.theta_filename = "theta"
        
        self.read_psi()
        self.read_theta()
        self.shape = (self.Npsi,self.Ntheta)
        self.psiN2D,self.theta2D = np.meshgrid(self.psiN,self.theta,indexing='ij')

        self.read_I()
        self.read_B()
        self.read_R()
        self.read_J()

        # differentiation matrices
        self.ddtheta=diff_matrix(self.theta[0],self.theta[0]+2*np.pi,self.Ntheta,order=4,periodic=True,endpoint=False)
        self.ddpsiN=diff_matrix(self.psiN[0],self.psiN[-1],self.Npsi,order=4,periodic=False,endpoint=True)

        self.ddtheta_B = np.einsum('kj,ij',self.ddtheta,self.B)
        self.ddpsiN_B = np.tensordot(self.ddpsiN,self.B,([1],[0]))
        self.ddpsiN_I = np.tensordot(self.ddpsiN,self.I,([1],[0]))

    def create_perfect_geometry(self,p_dirname,p_filename="input.namelist",norm_filename = "norms.namelist"):
        pinput = perfectInput(p_dirname + "/" + p_filename)
        
        normfile = f90nml.read(p_dirname + "/" + norm_filename)
        norm_group_name="normalizationParameters"
        BBar=normfile[norm_group_name]["BBar"]
        RBar=normfile[norm_group_name]["RBar"]
        TBar=normfile[norm_group_name]["TBar"]
        mBar=normfile[norm_group_name]["mBar"]
        nBar=normfile[norm_group_name]["nBar"]
        eBar=normfile[norm_group_name]["eBar"]
        ePhiBar=normfile[norm_group_name]["ePhiBar"]
        vBar=np.sqrt(2*TBar/mBar)
    

        pNtheta = pinput.Ntheta
        ptheta = pinput.theta

        pNpsi = pinput.Npsi        
        pinput.psiAHat = self.psia/(BBar*RBar**2)
        psiMin = pinput.psiMin
        psiMax = pinput.psiMax
        if pinput.psiGridType == 0:
            ppsiN = pinput.psi
        else:
            psiAHatFilename = pinput.psiAHatFilename
            psiN_to_psi_profiles=h5py.File(p_dirname + "/" + psiAHatFilename,'r')
            psiN_to_psi_groupname="/Npsi"+str(pNpsi)+"/"
            psiAHatArray = psiN_to_psi_profiles[psiN_to_psi_groupname+"psiAHatArray"][()]
            # these are the actual (nonuniform) psi and psiN
            ppsi = psiN_to_psi_profiles[psiN_to_psi_groupname+"psiArray"][()]
            ppsiN = psiN_to_psi_profiles[psiN_to_psi_groupname+"psiArray"][()]/self.psia
            internal_gradient_conversion_factor = psiAHatArray/self.psia
        
        IHat_interp = interp1d(self.psiN,self.I/(RBar*BBar))
        ddpsiN_IHat_interp = interp1d(self.psiN,self.ddpsiN_I/(RBar*BBar))
        
        # use %2pi on theta to go from [-pi,pi) to [0,2pi)
        # BHat_interp = interp2d(self.psiN2D,self.theta2D%(2*np.pi),self.B/(BBar))
        # ddtheta_BHat_interp = interp2d(self.psiN2D,self.theta2D%(2*np.pi),self.ddtheta_B/(BBar))
        # ddpsiN_BHat_interp = interp2d(self.psiN2D,self.theta2D%(2*np.pi),self.ddpsiN_B/(BBar))
        # RHat_interp = interp2d(self.psiN2D,self.theta2D%(2*np.pi),self.R/(RBar))
        # JHat_interp = interp2d(self.psiN2D,self.theta2D%(2*np.pi),self.J/(BBar/RBar))

        BHat_interp = interp2d(self.theta%(2*np.pi),self.psiN,self.B/(BBar))
        ddtheta_BHat_interp = interp2d(self.theta%(2*np.pi),self.psiN,self.ddtheta_B/(BBar))
        ddpsiN_BHat_interp = interp2d(self.theta%(2*np.pi),self.psiN,self.ddpsiN_B/(BBar))
        RHat_interp = interp2d(self.theta%(2*np.pi),self.psiN,self.R/(RBar))
        JHat_interp = interp2d(self.theta%(2*np.pi),self.psiN,self.J/(BBar/RBar))
        
        BHat = BHat_interp(ptheta,ppsiN)
        dBHatdpsi = ddpsiN_BHat_interp(ptheta,ppsiN)
        dBHatdtheta = ddtheta_BHat_interp(ptheta,ppsiN)
        RHat = RHat_interp(ptheta,ppsiN)
        JHat = JHat_interp(ptheta,ppsiN)
        IHat = IHat_interp(ppsiN)
        dIHatdpsi = ddpsiN_IHat_interp(ppsiN)

        pg_filename = pinput.geometryFilename
        pg = perfectGeometry(p_dirname + "/" + pg_filename)

        if pNpsi==1:
            BHat = BHat[np.newaxis,:]
            RHat = RHat[np.newaxis,:]
            JHat = JHat[np.newaxis,:]
            dBHatdpsi = dBHatdpsi[np.newaxis,:]
            dBHatdtheta = dBHatdtheta[np.newaxis,:]
            
        pg.create_geometry_for_Npsi_Ntheta(pNpsi, pNtheta, psiMin, psiMax, BHat, dBHatdpsi, dBHatdtheta, RHat, JHat, IHat, dIHatdpsi)

        
if __name__=="__main__":
    import matplotlib.pyplot as plt
    import sys

    argc = len(sys.argv)
    if argc == 3:
        neo_dir = sys.argv[1]
        perfect_dir  = sys.argv[2]
    elif argc == 2:
        neo_dir = sys.argv[1]
        perfect_dir = '.'
    else:
        raise ValueError("Need an eqdsk filename as an input argument")
    
    #neo_dir ='/home/bstefan/neoSimuls/t'
    #perfect_dir = "/home/bstefan/perfectSimuls/T17-09/JET89454_NU_geo_test"
    eg = EXPRO_geometry(neo_dir)
    # print "Printing B..."
    # plt.pcolor(eg.psiN2D,eg.theta2D,eg.B,rasterized=True)
    # plt.colorbar()
    # plt.xlabel(r"$\psi_N$")
    # plt.ylabel(r"$\theta$")
    # plt.savefig('EXPRO_B.pdf')
    print "Creating perfect geometry..."
    eg.create_perfect_geometry(perfect_dir)
