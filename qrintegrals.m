function [qs,rs,qv,rv]=qrintegrals(t1v,alphav,Phibarst,...
    M0,c1,Gammatilde,funcparams)
%
% Obtain q_s and r_s by integration.
%
% input variables:
%   as usual...
%   funcparams - theta function parameters ([ktemp,eps]).
%
% output variables:
%   q - q at initial time (t1v(1,1)), i.e., q_s in paper.
%   r - r at initial time (t1v(1,1)), i.e., r_s in paper.
%   qv - q_s for s in t1v(1,\cdot)
%   rv - r_s for s in t1v(1,\cdot)
%
% functions called:
%   gettheta
%
% Retrieve data
%
ntp1=size(t1v,2);
nt=ntp1-1;
del=t1v(1,2)-t1v(1,1);
invdel=1/del;
dim=size(M0,2);
M0p=M0';
%
% Get q by trapezoid rule.
%
% Horribly inefficient code in qrintegrals? Fix this?
%
% Horribly inefficient code!
%pause;
qv=zeros(dim,ntp1);
%qdotv4chk=zeros(2,ntp1);
%
% NEXT PART MAY NOT BE CORRECT?????
% SHOULDN"T "INTEGRATE SOLUTION" RATHER THAN ODE?
%
% Abarv(:,:,it)*qv(:,it)-c1+M0p*alphav4chk(1,it);
% q=0
% for is=nt:-1:1
%   q=q-del*(Abarv(:,:,is)*q-c1+M0p*alphav4chk(1,is))
% end;
%
for is=nt:-1:1
  if is==nt
    qv(:,is)=0.5*(Phibarst(:,:,is,is)*M0p*alphav(1,is)+...
        Phibarst(:,:,is,ntp1)*M0p*alphav(1,ntp1));
  else
    q=0.5*(Phibarst(:,:,is,is)*M0p*alphav(1,is)+...
        Phibarst(:,:,is,ntp1)*M0p*alphav(1,ntp1));
    for it=is+1:nt
      q=q+Phibarst(:,:,is,it)*M0p*alphav(1,it);
    end;
    qv(:,is)=q;
  end
end
%
% Be careful about signs,
%
qv=c1*del*qv;
qs=qv(:,1);
%
% Get Theta(alpha(t)).
%
[Thetav]=gettheta(alphav,funcparams);
%
% Get r by trapezoid rule.
%
rv=zeros(1,ntp1);
for is=nt:-1:1
  if is==nt
    rv(1,is)=0.5*(qv(:,is)'*Gammatilde*qv(:,is)+c1*alphav(1,is)^2-2*Thetav(1,is)+...
        qv(:,ntp1)'*Gammatilde*qv(:,ntp1)+c1*alphav(1,ntp1)^2-2*Thetav(1,ntp1));
  else
    r=0.5*(qv(:,is)'*Gammatilde*qv(:,is)+c1*alphav(1,is)^2-2*Thetav(1,is)+...
        qv(:,ntp1)'*Gammatilde*qv(:,ntp1)+c1*alphav(1,ntp1)^2-2*Thetav(1,ntp1));
    for it=is+1:nt
      r=r+qv(:,it)'*Gammatilde*qv(:,it)+c1*alphav(1,it)^2-2*Thetav(1,it);
    end;
    rv(1,is)=r;
    %if (is==1)||(is==2)
    %  [qv(:,is)'*Gammatilde*qv(:,is),c1*alphav(1,is)^2,2*Thetav(1,is)]
    %  'pausing in qrintegrals'
    %  pause;
    %end;
  end
end
%
% Be careful about signs.
%
%rv=del*rv;
rv=-del*rv;
rs=rv(1,1);
%'rdot in qrintegrals'
%[(rv(1,2)-rv(1,1))/del,del]
%'pausing.'
%pause;
