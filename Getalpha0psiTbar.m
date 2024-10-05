function [alpha0v]=Getalpha0psiTbar(M0,x,t1v,Phibarst)
%
% Obtain alpha0v = M0*Phibarstp(:,:,locTp1,is)*x;.
%
% input variables:
%
%   M0 - Matrix definining variable in nonlinearity
%   x - state value
%   t1v - (horizontal) vector of times
%   Phibarstp - same for \Abar'
%
% output variables:
%  alpha0v - time-indexed (hor) vector of alpha^T_\cdot values at s times.
%
% Recover data.
%
nfullp1=size(t1v,2);
%
% Get $\alpha_0(s;T) for "s"\in ]t_0,\Tbar[.
%
alpha0v=zeros(1,nfullp1);
for is=1:nfullp1
  alpha0v(1,is)=M0*squeeze(Phibarst(:,:,1,is))'*x;
end;

