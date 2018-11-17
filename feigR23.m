% feigR23.m
% Author: Satish Nallapaneni
% Revision April 21, 2014
% R must be greater than 1
% Bi must be greater than or equal to 0
% Eigen condition:
% f = S0*besselj(1,beta) - V0*bessely(1,beta) where,
% V0 = -beta*besselj(1,beta*R)+Bi*besselj(0,beta*R) and
% S0 = -beta*bessely(1,beta*R)+Bi*bessely(0,beta*R)
function BB=feigR23(num,R,Bi) %
x=R-1;
if Bi>0
for j= 1:num
eig21 = (j-1/2)*pi/x + ((0.267 - 0.007*x)/(1 + 0.5*x) + (0.05 + 0.04*x*j^(2/3))/(x + 1))/j^1.667 + (j - 1)*0.06/j^2;
n=j-1;
if n>0
eig22=n*pi/x + 0.058*x^0.25/n;
else
eig22=Bi/(100+Bi);
end
eigi=(eig22+1.25*Bi^0.6*eig21/j^1.5)/(1+1.25*Bi^0.6/j^1.5);
eigi=roundn(((eigi*10^10)/10^10),-30);
x0 = eigi;
dx=pi/x/10;
for pp=1:2
UL = (x0 + dx);
LL = (x0 - dx);
f0=eigenfunction(x0,R,Bi);
f1=eigenfunction(LL,R,Bi);
f2=eigenfunction(UL,R,Bi);
fp = ((f2 - f1)/2)/dx;
% fp: first derivative using first order central difference
fpp = ((f1 + f2 - 2*f0)/dx^2);
% fpp: second derivative using first order central difference
h = (-f0/fp);
eps=-fpp*h^2/2/(fp + h*fpp);
x0 = x0 + h + eps; % Newton Raphson method for finding eigenvalues
dx = h/5;
end
for n=1:3
V0x=-x0*besselj(1,x0*R)+Bi*besselj(0,x0*R);
S0x=-x0*bessely(1,x0*R)+Bi*bessely(0,x0*R);
fx=(S0x)*besselj(1,x0)-(V0x)*bessely(1,x0);
DV0=-R*(x0*besselj(0,x0*R)+Bi*besselj(1,x0*R));
DS0=-R*(x0*bessely(0,x0*R)+Bi*bessely(1,x0*R));
fpx=DS0*besselj(1,x0)+S0x*(besselj(0,x0)-besselj(1,x0)/x0)-DV0*bessely(1,x0)-V0x*(bessely(0,x0)-bessely(1,x0)/x0);
% fpx: actual first derivative of the eigencondition w.r.t. x0
x0=x0-(fx/fpx); % Newton Raphson method for finding eigenvalues
end
eig(j)=x0;
end
for ii=1:num
index(ii)=ii; x0=eig(ii); Bv(ii)=Bi; Ra(ii)=R;
fz(ii)=eigenfunction(x0,R,Bi);
n(ii)=x0/((ii-0.5)*pi/x);
end
else
for ii= 1:num
x0=1/(R-1)*(ii+.01)*pi;
for it=1:5
fxn=besselj(1,x0)*bessely(1,x0*R)-besselj(1,x0*R)*bessely(1,x0);% =0 for eigencondition
f1=(besselj(0,x0)-besselj(2,x0))*bessely(1,x0*R);
f2=R*besselj(1,x0)*(bessely(0,x0*R)-bessely(2,x0*R));
f2=f2-R*bessely(1,x0)*(besselj(0,x0*R)-besselj(2,x0*R));
f3=-(bessely(0,x0)-bessely(2,x0))*besselj(1,x0*R);% First derivative
fpxn=(f1+f2+f3)/2;
delt=fxn/fpxn; x0=x0-delt;%Newton Raphson method for finding eigenvalues
end%it
eig(ii)=x0;
fz(ii)=besselj(1,x0)*bessely(1,x0*R)-besselj(1,x0*R)*bessely(1,x0);
index(ii)=ii;Bv(ii)=0; Ra(ii)=R;
n(ii)=eig(ii)/pi*(R-1)/ii;
end%ii
end %if Bi>0
%sprintf(' index    R    Biot#    beta    m    Root Value')
%BBB=[index' Ra' Bv' eig' n' fz'];
%fprintf('%5.0f %10.5f %5.3f %12.10f %2.5f %12.5e\n',BBB')
BB=eig;
