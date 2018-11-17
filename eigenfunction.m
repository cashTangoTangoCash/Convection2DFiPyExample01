#authors: Satish Nallapaneni & James V. Beck – June 29, 2014
% function that calculates eigencondition for a given eigenvalue
function f=eigenfunction(beta,R,Bi)
S0 = -beta*bessely(1, beta*R) + Bi*bessely(0, beta*R);
V0 = -beta*besselj(1, beta*R) + Bi*besselj(0, beta*R);
f = (S0*besselj(1, beta) - V0*bessely(1, beta));
