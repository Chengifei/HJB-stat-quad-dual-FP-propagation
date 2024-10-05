function [alphav,Phibarst,Phibarstp,Dssig,Dpsigs,alpha0v]=...
    GetalphaPhifuncnlpsi(M0,c1,x,t1v,Psiv,first,alphastartv,...
    funcparams,whichmethod)
%
% N.B.: ASSUMES $\Theta(\alpha)=\atan(\alpha)$!!!!!!
%
%
% Obtain Phibarst, Phibarstp, alpha0v and alphav.
%
% input variables:
%
%  M0 - Matrix definining variable in nonlinearity
%  c1 - stat-quad coef.
%  x - state value
%  t1v - (horizontal) vector of times
%  Abarv - array of time-dependent linear dynamics asociated to q
%  Psiv - Fundamental matrix function (vector) associated Abarv.
%  first - 1 if first time-step; 0 otherwise
%  alphastartv - alpha function vector for start of fixed-point iteration.
%               This is ignored if first==1.
%   funcparams - the parameters ([ktemp,eps]).
%   whichmethod - 1 for fixed-point method 1 and 2 for method 2.
%
% output variables:
%  alphav - time-indexed (hor) vector of alpha^T(\cdot) values.
%  Phibarst - state-trans matrix for $\Abar$ from t to s - obtained for
%             all times (up through the terminal terminal time).
%  Phibarstp - same for \Abar'
%  alpha0v - time-indexed (hor) vector of alpha values at init times.
%
% functions called:
%   getDs,
%   Getalphafuncnl[gradtheta,newtonforeta]
%
%
% First, recover data.
%
M0p=M0';
t0v=t1v;
ntp1=size(t1v,2);
nt=ntp1-1;
t0=t1v(1,1);
t1=t1v(1,ntp1);
t0v=t1v;
del=t1v(1,2)-t1v(1,1);
invdel=1/del;
%
% Get $\Phibar_\tau^s$ for all generic \tau, s.
% Also get $(\Phibar_\tau^s)'$.
% NOTE the ordering matches the subscript ordering of state transition
% matrices.
%
% Also, getting these for the entire duration, rather than only up to the
% initial terminal time, t10.
%
nfullp1=size(Psiv,3);
Phibarst=zeros(2,2,nfullp1,nfullp1);
Phibarstp=zeros(2,2,nfullp1,nfullp1);
Psiinvv=zeros(2,2,nfullp1);
for it=1:nfullp1
    Psiinvv(:,:,it)=inv(Psiv(:,:,it));
end
for it0=1:nfullp1
  for it1=1:nfullp1
      Phibarst(:,:,it1,it0)=Psiv(:,:,it1)*Psiinvv(:,:,it0);
      Phibarstp(:,:,it1,it0)=(Phibarst(:,:,it1,it0))';
  end
end
%
% Get Dssig and DPsigs
%
[Dssig,Dpsigs]=getDs(x,t1v,Phibarst,Phibarstp);
%
% Get $\alpha_{t_0}^T (indep of T)
%
alpha0v=zeros(1,ntp1);
for it0=1:ntp1
  alpha0v(1,it0)=M0*Phibarstp(:,:,ntp1,it0)*x;
end;
%
% Get alpha^T_s for all terminal times, T,
% and intermediary times, s, using (8.5).
% Currently uses a fixed-point method.
%
if first==1
  [alphav]=Getalphafuncnl(M0,c1,x,t1v,Dssig,Dpsigs,...
      alpha0v,alphastartv,funcparams,whichmethod);
else
  [alphav]=Getalphafuncnl(M0,c1,x,t1v,Dssig,Dpsigs,...
      alpha0v,alphastartv,funcparams,whichmethod);
end
%
% Plot the resulting alpha function, if desired.
%
plotalp=0;
if plotalp==1
  'Plotting from within GetalphaPhifuncnlpsi.'
  figure;
  plot(t1v,alphav);
end;