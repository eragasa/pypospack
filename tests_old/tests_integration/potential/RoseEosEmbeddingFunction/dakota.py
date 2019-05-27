#DO NOT CALL DIRECTLY, load as module.
# Script: dakota.py
# Purpose: Module providing support scripts for Dakota Fitting
# Syntax:  
# Example: 
# Author:  cjo 2/13/15
from numpy import *
from scipy.optimize import brentq
globalcutoff = 5.75

def psi(r):
    ''' Implements cutoff function for smoothly bringing function to zero '''
    if type(r) == ndarray:
        s = empty_like(r)
        for i in xrange(len(r)): 
            x = r[i]
            if x > 1.0:
                s[i]= 0.
            elif ((x > 0.0) and (x <= 1.0)):
                s[i] = ( -6.*x**5 + 15.*x**4 - 10.*x**3 + 1.)
            else:
                s[i] = 1.
    else:
        if r > 1.0:
            s = 0.
        elif ((r > 0.0) and (r <= 1.0)):
            s = ( -6.*r**5 + 15.*r**4 - 10.*r**3 + 1.)
        else:
            s = 1.

    return s

def fpair(r,params):
    ''' implemnts Morse functional form of pair potential '''
    De,am,r0,rp = params
    pairvals = De * ( ( 1.-exp( -am*(r-r0) ) )**2. - 1. ) * psi( (r-rp)/(globalcutoff-rp) )
    return pairvals

def fdens(r,params):
    ''' implements Exponential e- density functional form '''
    rho0,lambda0,r0,rd = params
    densvals = rho0 * exp( -1.*( r - r0 )/lambda0 )* psi( (r-rd)/(globalcutoff-rd))
    return densvals

def fembedbjs(r,params):
    ''' implements BJS functional form of embedding function '''
    F0,gamma,F1 = params
    embedvals = F0 * ( 1. - gamma*log(r+1.e-30) ) * r**gamma + F1 * r
    return embedvals

def rhofxn(a,rho0,r0,lambda0,rd,rhostar):
    ''' calculates ideal e- density based on exponential functional form 
        data input format:  rho0  r0    lambda0  rd   rhostar ''' 
    return rho0*(12.*exp(-(a/sqrt(2.)-r0)/lambda0) * psi( (a/sqrt(2.)-rd) / (globalcutoff-rd) )
            + 6.*exp(-(a-r0)/lambda0) * psi( (a-rd) / (globalcutoff-rd) )
            + 24.*exp(-(a*sqrt(1.5)-r0)/lambda0) * psi( (a*sqrt(1.5)-rd) / (globalcutoff-rd) )
            + 12.*exp(-(a*sqrt(2.)-r0)/lambda0) * psi( (a*sqrt(2.)-rd) / (globalcutoff-rd) )
            + 24.*exp(-(a*sqrt(2.5)-r0)/lambda0) * psi( (a*sqrt(2.5)-rd) / (globalcutoff-rd) )
            + 8.*exp(-(a*sqrt(3.)-r0)/lambda0) * psi( (a*sqrt(3.)-rd) / (globalcutoff-rd) )
            ) - rhostar

def fembedFoiles(rho,params):
    ''' implements Foiles-style embeding function 
        (i.e. density and pair potential forced to match Rose EOS)
        parameter list:
        p[0]   p[1]     p[2] p[3] p[4]    p[5]       p[6]   p[7] p[8] 
        E_coh  a(morse) r0  rho0 lambda0 lambdarose De      rp   rd  '''
    Ecoh,am,r0,rho0,lambda0,lambdarose,De,rp,rd = params
    embedvals = empty_like(rho)
    k=0
    for rhostar in rho:
        #solve density for lattice constant (a) where density (rhostar) is found
        rhop = (rho0,r0,lambda0,rd,rhostar)
        a = brentq(rhofxn,0.,10000.,rhop,xtol=1.0e-8) #lattice constant where rhostar is found
        #find E_Rose for lattice constant
        astar = (a-r0*sqrt(2.)) / (lambdarose*sqrt(2.)*r0) 
        Erose = Ecoh*(1+astar)*exp(-astar) 
        #find pair potential for lattice constant
        pp = (De,am,r0,rp)
        Epot = 12.*fpair(a/sqrt(2.),pp)+6.*fpair(a,pp)+24.*fpair(sqrt(1.5)*a,pp)+12.*fpair(sqrt(2.)*a,pp)+24.*fpair(sqrt(2.5)*a,pp)+8.*fpair(sqrt(3.)*a,pp)
        #calculate value of embedding fxn
        embedvals[k] = Erose - 0.5*Epot
        k += 1

    return embedvals

def write_eamhead(eam,params,title,*elist):
    ''' write header for eam file'''
    Nrho,drho,Nr,dr = params
    eam.write('Author: cjo & sf Citation: Sandia National Laboratories, Albuquerque NM\n')
    eam.write(title+'\n')
    eam.write('***********************\n')
    #write num atoms
    if len(elist)==1:
        eam.write('    1 {0:s}\n'.format(elist[0]))
    elif len(elist)==2:
        eam.write('    2 {0:s} {1:s}\n'.format(elist[0],elist[1]))
    else:
        raise StandardError, "Only two elements can be used at this time!!!"

    eam.write('{0:5d}{1:24.16e}{2:5d}{3:24.16e}{4:24.16f}\n'.format(Nrho,drho,Nr,dr,globalcutoff))

def write_singleeam(eam,el):
    ''' Write Single atom section of  EAM/ALLOY style (funcfl) Potential for LAMMPS'''
    ###################### Write Single Element Section #####################################
    #for el in [Pt,Au,Pd]:
    eam.write('{0:5d}{1:15.5g}{2:15.5g} {3:s}\n'.format(el.atnum,el.mass,el.alat,el.lattype))
    #write embedding fxn for Nrho values
    for i in xrange(0,len(el.embed),5):
        eam.write('{0:24.16E}{1:24.16E}{2:24.16E}{3:24.16E}{4:24.16E}\n'.format(
            el.embed[i],el.embed[i+1],el.embed[i+2],el.embed[i+3],el.embed[i+4]))
    #write density for Nr values
    for i in xrange(0,len(el.rho),5):
        eam.write('{0:24.16E}{1:24.16E}{2:24.16E}{3:24.16E}{4:24.16E}\n'.format(
            el.rho[i],el.rho[i+1],el.rho[i+2],el.rho[i+3],el.rho[i+4]))

def write_paiream(eam,r,pair):
    ''' Write pair section of EAM/ALLOY style (funcfl) Potential for LAMMPS'''
    ###################### Write Pair Potential Section #####################################
    #for pair in [phiPtPt,phiAuPt,phiAuAu,phiPdPt,phiPdAu,phiPdPd]:
    #pair = phiAuAu
    pairr=r[:]*pair[:]
    #write pair fxn for Nr values * r
    for i in xrange(0,len(pairr),5):
        eam.write('{0:24.16E}{1:24.16E}{2:24.16E}{3:24.16E}{4:24.16E}\n'.format(
            pairr[i],pairr[i+1],pairr[i+2],pairr[i+3],pairr[i+4]))


def getparams(dakotain):
    ''' Creates dictionary of parameters provided by Dakota input to analysis script:
        Is agnostic to number and name of parameters '''
    infile  = open(dakotain,'r')
    plist=[] #dictionary values

    for line in infile: 
        line.strip()
        if line.split()[1] == 'functions':
            break
        elif line.split()[1] =='variables':
            continue
        else:
            plist.append( (line.split()[1],float(line.split()[0])) )

    pdict = dict(plist)
    
    infile.close()
    return pdict
