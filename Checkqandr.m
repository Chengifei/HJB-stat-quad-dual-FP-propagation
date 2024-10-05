function [qchkv,rchkv,relqchkv,relrchkv]=Checkqandr(qv,rv,tv,Pv,alphav,...
    Amat,Cmat,M0,Gammatilde,c1,Phibarst,Abarv,funcparams)
%
% Check that q,r actualy satisfy their ODEs.
%
% input variables:
%   qv - 2 x nfull vector of q_\cdot for times in tv
%   rv - 1 x nfull vector of r_\cdot for times in tv
%   tv - 1 x nfull vector of times
%   Pv - 2 x 2 x nfull solution of DRE
%   alphav - 1 x nfull alpha control values
%   x1 - x_1
%   Amat,Cmat,M0,Ga0mmatilde, c1 - parameters
%   funcparams - theta function parameters ([ktemp,eps]).
%
% output variables:
%   qchkv - qdot minus the RHS of the ODE- should be zeros
%   rchkv - qdot minus the RHS of the ODE- should be zeros
%
% functions called:
%   gettheta
%
%
ntp1=size(tv,2);
nt=ntp1-1;
del=tv(1,2)-tv(1,1);
twodel=2*del;
%
% Get theta values.
%
[thetav]=gettheta(alphav,funcparams);
%
% Get time derivatives.
%
qdotv=zeros(2,ntp1);
rdotv=zeros(1,ntp1);
qdotv(:,1)=(qv(:,2)-qv(:,1))/del;
qdotv(:,ntp1)=(qv(:,ntp1)-qv(:,nt))/del;
rdotv(1,1)=(rv(1,2)-rv(1,1))/del;
rdotv(1,ntp1)=(rv(1,ntp1)-rv(1,nt))/del;
for it=2:nt
  qdotv(:,it)=(qv(:,it+1)-qv(:,it-1))/twodel;
  rdotv(1,it)=(rv(1,it+1)-rv(1,it-1))/twodel;
end;
%
% Do the checks.
%
Amatp=Amat';
M0p=M0';
qchkv=zeros(2,ntp1);
relqchkv=zeros(1,ntp1);
rchkv=zeros(1,ntp1);
relrchkv=zeros(1,ntp1);
for it=1:ntp1
  Abar=Pv(:,:,it)*Gammatilde-Amatp;
  Rhsq=Abar*qv(:,it)-c1*M0p*alphav(1,it);
  %Rhsq2=-1*Phibarst(:,:,it,it)*c1*M0p*alphav(1,it)+...
  %    Abarv(:,:,it)*qv(:,it);
  qchkv(:,it)=qdotv(:,it)-Rhsq;
  relqchkv(1,it)=norm(qchkv(:,it),2)/(norm(qdotv(:,it),2)+norm(Rhsq,2));
  %
  qGq=qv(:,it)'*Gammatilde*qv(:,it);
  aCa=c1*alphav(1,it)^2;
  Rhsr=qGq+aCa-2*thetav(1,it);
  rchkv(1,it)=rdotv(1,it)-Rhsr;
  %rdot11=rdotv(:,it);
  %[rdot11,Rhsr]
  relrchkv(1,it)=abs(rchkv(1,it))/(abs(rdotv(1,it))+abs(Rhsr));
  %rdot=rdotv(1,it);
  %if (it==1)||(it==2)
  %  [qGq,aCa,2*thetav(1,it),rdot,Rhsr*del]
  %  'Pausing in Checkqandr'
  %  pause;
  %end;
end;
