# -*- coding: utf-8 -*-

import sys
import numpy as np
import scipy.linalg as sc

#=============================================================================
# 1. Assembling the stiffness matrix, with option for 2nd order effect
#-----------------------------------------------------------------------------

def stiffness(L, EI, P=None):
    ''' Global stiffness matrix for a discretized column.
        All input parameters are numpy vectors, one element per 
        column discretization element.
        
        L:   lengths
        EI:  stifnesses
        P:   axial loading (negative for compression)
    '''

# Preliminary settings
    
    N   =  len(L) + 1
    KG  =  np.zeros((2*N,2*N))

    ke11 = 12*EI/L/L/L
    ke12 =  6*EI/L/L
    ke22 =  4*EI/L
    ke24 =  2*EI/L

    if (P is None):
        kg11 = 36*P/L/30
        kg12 =  3*P  /30*np.ones(L.shape)
        kg22 =  4*P*L/30
        kg24 =    P*L/30

# Iterate for all discretization elements
        
    for k in range(N-1):
    
        i1 = 2* k           #0
        i2 = 2* k + 1       #1
        i3 = 2*(k + 1)      #2
        i4 = 2*(k + 1) + 1  #3

# Linear elastic stifness matrix
    
        k11k = ke11[k] # 12
        k12k = ke12[k] # 6L
        k22k = ke22[k] # 4L^2
        k24k = ke24[k] # 2L^2
        
        KG[i1,i1] += k11k
        KG[i1,i2] += k12k;    KG[i2,i1] += k12k
        KG[i1,i3] -= k11k;    KG[i3,i1] -= k11k
        KG[i1,i4] += k12k;    KG[i4,i1] += k12k
        
        KG[i2,i2] += k22k;
        KG[i2,i3] -= k12k;    KG[i3,i2] -= k12k
        KG[i2,i4] += k24k;    KG[i4,i2] += k24k
        
        KG[i3,i3] += k11k
        KG[i3,i4] -= k12k;    KG[i4,i3] -= k12k
        KG[i4,i4] += k22k

# Linear elastic stifness matrix

        if (P is None):    
            k11k = kg11[k] # 36
            k12k = kg12[k] # 3L
            k22k = kg22[k] # 4L^2
            k24k = kg24[k] # L^2
            
            KG[i1,i1] += k11k
            KG[i1,i2] += k12k;    KG[i2,i1] += k12k
            KG[i1,i3] -= k11k;    KG[i3,i1] -= k11k
            KG[i1,i4] += k12k;    KG[i4,i1] += k12k
            
            KG[i2,i2] += k22k
            KG[i2,i3] -= k12k;    KG[i3,i2] -= k12k
            KG[i2,i4] -= k24k;    KG[i4,i2] -= k24k
            
            KG[i3,i3] += k11k
            KG[i3,i4] -= k12k;    KG[i4,i3] -= k12k
            KG[i4,i4] += k22k

    return KG

#=============================================================================
# 2. Assembling consistent mass matrix
#-----------------------------------------------------------------------------

def consistMass(L, mu):
    ''' Global consistent mass matrix for a discretized column.
        All input parameters are numpy vectors, one element per 
        column discretization element.
        
        L:   lengths
        mu:  mass per unit lengths
    '''

# Preliminary settings

    N   =  len(L) + 1
    MG  =  np.zeros((2*N,2*N))
    
    m11 = 156*mu*L    /420
    m12 =  22*mu*L*L  /420
    m13 =  54*mu*L    /420
    m14 =  13*mu*L*L  /420
    m22 =   4*mu*L*L*L/420
    m24 =   3*mu*L*L*L/420

# Iterate for all discretization elements
        
    for k in range(N-1):
    
        i1 = 2* k           #0
        i2 = 2* k + 1       #1
        i3 = 2*(k + 1)      #2
        i4 = 2*(k + 1) + 1  #3
    
        m11k = m11[k]
        m12k = m12[k]
        m13k = m13[k]
        m14k = m14[k]
        m22k = m22[k]
        m24k = m24[k]
 
        MG[i1,i1] += m11k
        MG[i1,i2] += m12k;    MG[i2,i1] += m12k
        MG[i1,i3] += m13k;    MG[i3,i1] += m13k
        MG[i1,i4] -= m14k;    MG[i4,i1] -= m14k
        
        MG[i2,i2] += m22k;
        MG[i2,i3] += m14k;    MG[i3,i2] += m14k
        MG[i2,i4] -= m24k;    MG[i4,i2] -= m24k
        
        MG[i3,i3] += m11k;    
        MG[i3,i4] -= m12k;    MG[i4,i3] -= m12k
        MG[i4,i4] += m22k
     
    return MG

#=============================================================================
# 3. Complete case analysis
#-----------------------------------------------------------------------------

def analyseCase(H, EI, mu, P1, n):

    L = (H/n)*np.ones(n)
    P =  P1 + np.cumsum(L*mu)
    
    K =  stiffness(L, EI, P)
    M =  consistMass(L, mu)

    K       =  K[:-2,:-2]
    M       =  M[:-2,:-2]
    M[0,0] +=  P1/9.81

    w2, Phi = sc.eig(K, M)
    iw      = w2.argsort()
    w2      = w2[iw]
    Phi     = Phi[:,iw]
    wk      = np.sqrt(np.real(w2)) 
    fk      = wk/2/np.pi
    
    return fk[0]

#=============================================================================
