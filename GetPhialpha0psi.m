function [alpha0v]=GetPhialpha0psi(M0,c1,x,t1v,Phibarst,Phibarstp)
%
% Obtain Phibarst, Phibarstp and alpha0v.
%
% First, recover data.
%
%M0p=M0';
ntp1=size(t1v,2);
%nt=ntp1-1;
%t0=t1v(1,1);
%t1=t1v(1,ntp1);
%t0v=t1v;
%del=t1v(1,2)-t1v(1,1);
%
% Get $\alpha_{t_0}^T (indep of T)
%
alpha0v=zeros(1,ntp1);
for it0=1:ntp1
  alpha0v(1,it0)=M0*Phibarstp(:,:,ntp1,it0)*x;
end;
%'In GetPhialpha0'
%alpha0v

