function [backsubserr,relbacksuberr]=backsubsalpnlTbar(alphachkv,...
    M0,c1,dim,t1v,DMz,alpha0v,locTp1,funcparams)
%
% N.B.: ASSUMES $\Theta(\alpha)=\atan(\alpha)$!!!!!!
%
% Backsubstitution check of our alpha vector
% (representing the alpha function of time).
% This is the nonlinear version, i.e., for $\grad\Theta\not=0$.
%
%
% inut variables:
%   alphachkv - the alpha vector to be checked.
%   alpha0v - current initial alpha (at t=0).
%   M0,c1,dim,t1v - parameters
%   Dssig,Dpsigs - internal variables.
%   T - We have \alpha^T_s for t\le s,T \le \Tbar. Here, typically we
%       take t=0. Also, we're really only interested for s\in[t,T];
%       the use of \Tbar with s\in (T,\Tbar] is a mathematical
%       artifact that is useful for obtaing the function over [t,T].
%   locTp1 - the position of T in t1v, i.e., where T=t1v(1,locT+1).
%   funcparams - the parameters ([ktemp,eps]).
%
%
% functions called:
%   gradtheta
%
% Recover data.
%
M0p=M0';
invc1=1/c1;
t0v=t1v;
nfullp1=size(t1v,2);
nfull=nfullp1-1;
%t0=t1v(1,1);
%t1=t1v(1,ntp1);
locT=locTp1-1;
locTm1=locT-1;
del=t1v(1,2)-t1v(1,1);
delo2=del/2;
%dim=size(x,1);
%
% First, we need to get the integral term. This integral only
% goes up to T, but should be for DMz with s,\sigma values in [1,nfullp1].
%
mc1delM0=-c1*del*M0;
intforalp=zeros(1,nfullp1);
for is=1:nfullp1
    if locTp1==1
        intforalp(1,is)=0;
    elseif locTp1==2
      intforalp(1,is)=0.5*mc1delM0*(DMz(:,1,is)*alphachkv(1,1)+...
        DMz(:,locTp1,is)*alphachkv(1,locTp1));
    else
      intforalp(1,is)=0.5*mc1delM0*(DMz(:,1,is)*alphachkv(1,1)+...
        DMz(:,locTp1,is)*alphachkv(1,locTp1));
      for isig=2:locT
        intforalp(1,is)=intforalp(1,is)+mc1delM0*DMz(:,isig,is)*...
            alphachkv(1,isig);
      end;
    end;  % end of if locTp1\in...
end; % end of for is=...
%'In backsubsalpnlTbar:'
%intforalpold
%intforalp
%pause;
%
% Get $\grad\Theta(\alpha)$ for alpha in alphachkv.
%
[gradthetav]=gradtheta(alphachkv,funcparams);
%
% Obtain the backsubstitution error.
%
backsubserr=alphachkv-(invc1*gradthetav+alpha0v+intforalp);
%'componennts of check in back...nl'
%[alphachkv;invc1*gradthetav;alpha0v;intforalp]
normbse=norm(backsubserr,2);
%denomalp=norm(alphachkv,2);
denomalp=norm(alphachkv,2)+norm(invc1*gradthetav,2)+norm(alpha0v,2)+...
    norm(intforalp,2);
%
% Avoid near-zero (or zero!) denonminator issue.
%
denomtoolowbnd=1e-10;
relbacksubserrreplacement=1e-5;
if denomalp<denomtoolowbnd
  'Denominator in rel diff comp in backs...Tbar near zero; using alternate.'
  relbacksuberr=relbacksubserrreplacement;
else
  relbacksuberr=normbse/denomalp;
end;
