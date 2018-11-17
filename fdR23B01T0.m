% fdR23B01T0 function
% Author: Satish Nallapaneni
% Revision June 27, 2014
% INPUTS:
% R: radius ratio R2/R1
% rv: dimensionless location starting at rd=r/R1=1 and ending at rd=R2/R1=R
% tv: dimensionless time starting at td=0
% A = desired accuracy (1E-A =10^-A); A=2,3, ..., 15
% Bi: Biot number
% OUTPUTS:
% Td: dimensionless temperature calculated at (xd,td) to desired accuracy A
% qd: dimensionless heat flux calculated at (xd,td) to desired accuracy A
function [Td,qd]=fdR23B01T0(rv,tv,R,Bi,A)
pkg load mapping
  srv=length(rv);
stv=length(tv);
Temp=zeros(stv,srv);
% Preallocating arrays for speed
flux=zeros(stv,srv);
% Preallocating arrays for speed
% calculate number of eigenvalues required to obtain solution with accuracy A
mmax1=floor(2*(0.5+(R-1)/pi*sqrt(A*log(10)/min(tv))));
bet=feigR23(mmax1,R,Bi);% Call the function to get eigenvalues
for ir = 1:srv
% begin space loop
r=rv(ir);
term1=1; % steady-state temperature solution
term2=0;
% steady-state heat flux solution
for it=1:stv
% begin time loop
t=tv(it);
% calculate m_max for every timestep as it reduces computational cost
mmax=floor(1.2*(0.5+(R-1)/pi*sqrt(A*log(10)/t)));
term3 = 0;
term4 = 0;
for ii=1:mmax
bt=bet(ii);
V0=-bt*besselj(1,bt*R)+Bi*besselj(0,bt*R);
S0=-bt*bessely(1,bt*R)+Bi*bessely(0,bt*R);
Nr = (S0*besselj(0,bt*r)-V0*bessely(0,bt*r));
L1 = (S0*besselj(0,bt*R)-V0*bessely(0,bt*R));
Denominator = (Bi^2+bt^2)*besselj(1,bt)*besselj(1,bt)-V0^2;
Norm = Denominator/(besselj(1,bt)*besselj(1,bt));
term3 = term3 - (pi*pi/2)*Bi*R*(exp(-bt*bt*t))*Nr*L1/Norm;
Nm = (S0*besselj(1,bt*r)-V0*bessely(1,bt*r));
term4 = term4 - (pi*pi/2)*Bi*R*exp(-bt*bt*t)*bt*Nm*L1/Norm;
end
T(ir)=term1;
Q(ir)=term2;
Temp(it,ir) = term1 + term3; % total temperature solution
flux(it,ir) = term2 + term4; % total heat flux solution
end
end
Td=abs(Temp); qd=abs(flux);
