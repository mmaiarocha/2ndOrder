function K = Tower2DStiff(K,L,n1,n2,I,E,P);
%TOWER2DSTIFF Add a 2D beam segment to tower model stiffness matrix.
%        K = Tower2DStiff(K,L,n1,n2,I,E,P) sum up a 2D beam stiffness
%        matrix to a global stiffness matrix. The required input is:
%
%        K:  the global stiffness matrix.
%        L:  the beam length.
%        n1: number of starting structural node.
%        n2: number of ending structural node.
%        I:  moment of inertia of beam cross section.
%        E:  Young's modulus of material.
%        P:  axial load for considering geometric stiffness
%           (where negative values mean compression).
%
%   Programmed by:
%   Marcelo M. Rocha
%   LDEC/PPGEC/UFRGS
%   10-Jan-2002.
%
%========================================================================

%------------------------------------------------------------------------
% 1. Address in the global matrix

i1  = 2*n1 - 1;
i2  = 2*n1;
i3  = 2*n2 - 1;
i4  = 2*n2;

% Linear elastic stifness matrix

k11 = 12*E*I/L/L/L;
k12 =  6*E*I/L/L;
k22 =  4*E*I/L;

K(i1,i1) = K(i1,i1) + k11;
K(i1,i2) = K(i1,i2) + k12;    K(i2,i1) = K(i2,i1) + k12;
K(i1,i3) = K(i1,i3) - k11;    K(i3,i1) = K(i3,i1) - k11; 
K(i1,i4) = K(i1,i4) + k12;    K(i4,i1) = K(i4,i1) + k12; 
K(i2,i2) = K(i2,i2) + k22;
K(i2,i3) = K(i2,i3) - k12;    K(i3,i2) = K(i3,i2) - k12;
K(i2,i4) = K(i2,i4) + k22/2;  K(i4,i2) = K(i4,i2) + k22/2;
K(i3,i3) = K(i3,i3) + k11;    
K(i3,i4) = K(i3,i4) - k12;    K(i4,i3) = K(i4,i3) - k12;
K(i4,i4) = K(i4,i4) + k22;

% Geometric stifness matrix

kg11 = 36*P/30/L;
kg12 =  3*P/30;
kg22 =  4*P*L/30;

K(i1,i1) = K(i1,i1) + kg11;
K(i1,i2) = K(i1,i2) + kg12;    K(i2,i1) = K(i2,i1) + kg12;
K(i1,i3) = K(i1,i3) - kg11;    K(i3,i1) = K(i3,i1) - kg11; 
K(i1,i4) = K(i1,i4) + kg12;    K(i4,i1) = K(i4,i1) + kg12; 
K(i2,i2) = K(i2,i2) + kg22;
K(i2,i3) = K(i2,i3) - kg12;    K(i3,i2) = K(i3,i2) - kg12;
K(i2,i4) = K(i2,i4) - kg22/4;  K(i4,i2) = K(i4,i2) - kg22/4;
K(i3,i3) = K(i3,i3) + kg11;
K(i3,i4) = K(i3,i4) - kg12;    K(i4,i3) = K(i4,i3) - kg12;
K(i4,i4) = K(i4,i4) + kg22;

%======================================================================
% End of function m-file TOWER2DSTIFF.


