function M = Tower2DMass(M,L,n1,n2,mi);
%TOWER2DMASS Add a 2D beam segment to tower model mass matrix.
%        M = Tower2DMass(M,L,n1,n2,mi) sum up a 2D beam consistent mass
%        matrix to a global mass matrix. The required input is:
%
%        M:  the global mass matrix.
%        L:  the beam length.
%        n1: number of starting structural node.
%        n2: number of ending structural node.
%        mi: the element mass per unit length.
%
%   Programmed by:
%   Marcelo M. Rocha
%   LDEC/PPGEC/UFRGS
%   12-Apr-2002.
%
%========================================================================

%------------------------------------------------------------------------
% 1. Address in the global matrix

i1  = 2*n1 - 1;
i2  = 2*n1;
i3  = 2*n2 - 1;
i4  = 2*n2;

% Consistent mass matrix

m11 = 156*mi*L/420;
m12 =  22*mi*L*L/420;
m13 =  54*mi*L/420;
m14 =  13*mi*L*L/420;
m22 =   4*mi*L*L*L/420;
m24 =   3*mi*L*L*L/420;

M(i1,i1) = M(i1,i1) + m11;
M(i1,i2) = M(i1,i2) + m12;    M(i2,i1) = M(i2,i1) + m12;
M(i1,i3) = M(i1,i3) + m13;    M(i3,i1) = M(i3,i1) + m13; 
M(i1,i4) = M(i1,i4) - m14;    M(i4,i1) = M(i4,i1) - m14; 
M(i2,i2) = M(i2,i2) + m22;
M(i2,i3) = M(i2,i3) + m14;    M(i3,i2) = M(i3,i2) + m14;
M(i2,i4) = M(i2,i4) - m24;    M(i4,i2) = M(i4,i2) - m24;
M(i3,i3) = M(i3,i3) + m11;    
M(i3,i4) = M(i3,i4) - m12;    M(i4,i3) = M(i4,i3) - m12;
M(i4,i4) = M(i4,i4) + m22;

%======================================================================
% End of function m-file TOWER2DMASS.


