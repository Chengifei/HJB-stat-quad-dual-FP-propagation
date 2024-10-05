function [alphafullv]=...
    GetalpPhinlpsiTbar(M0,c1,t1v,Phibarst,DMz,...
    nt0,alpha0v,alphastartv,funcparams)
%
% N.B.: ASSUMES $\Theta(\alpha)=\atan(\alpha)$!!!!!!
%
%
% Obtain Phibarst, Phibarstp, alpha0v and alphav.
%
% input variables:
%
%   M0 - Matrix definining variable in nonlinearity
%   c1 - stat-quad coef.
%   t1v - (horizontal) vector of times
%         Note: T="last element of t1v".
%         \alpha^T_s for t\le s,T \le \Tbar. Here, typically we
%         take t=0. Also, we're really only interested for s\in[t,T];
%         the use of \Tbar with s\in (T,\Tbar] is a mathematical
%         artifact that is useful for obtaing the function over [t,T].
%   Abarv - array of time-dependent linear dynamics asociated to q
%   Phibarst - state-trans matrix for $\Abar$ from t to s - obtained for
%             all times (up through the terminal terminal time).
%   Phibarstp - same for \Abar'
%   nt0 - The step at which we hit t10 (t10=t0+del*nt0).
%   alpha0v - time-indexed (hor) vector of alpha values at init times.
%   alphastartv - alpha function vector for start of fixed-point iteration.
%               This is ignored if first==1.
%   funcparams - the parameters ([ktemp,eps]).
%
% output variables:
%  alphav - time-indexed (hor) vector of alpha^T(\cdot) values.
%
% functions called:
%   getDs,
%   Getalphafuncnl[gradtheta,newtonforeta],
%   extendfixedptalp[etainverse]
%
%
% First, recover data.
%
M0p=M0';
nfullp1=size(t1v,2);
nfull=nfullp1-1;
ntp1=nt0+1;
locT=nt0;
t1vshort=t1v(1,1:ntp1); %short version only goes up to ntp1 (i.e., to t10).
t0=t1v(1,1);
t10=t1v(1,ntp1);
del=t1v(1,2)-t1v(1,1); % step size
invdel=1/del;
%
% Get $\alpha_{t_0}^T (shorter vector)
%
alpha0vshort=alpha0v(1,1:ntp1);
%
% Get alpha^T_s for all terminal times, T,
% and intermediary times, s, using (8.5).
% Currently uses a fixed-point method.
% Method 2 is the one [that's guaranteed to work for T sufficiently small.
%
whichmethod=2;  % This is the only option!
[alphavshort]=Getalphafuncnl(M0,c1,t1vshort,...
    DMz,alpha0vshort,alphastartv,funcparams);
%
% Now, we need to extend the fixed-point solution to the
% entire interval (up to \Tbar).
%
[alphafullv]=extendfixedptalp(alphavshort,alpha0v,M0,c1,t1v,nt0,...
    Phibarst,DMz,funcparams);
%
% Plot the resulting alpha function, if desired.
%
plotalp=0;
if plotalp==1
  'Plotting from within GetalphaPhifuncnlpsi.'
  figure;
  plot(t1v,alphav);
end;