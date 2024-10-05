function [alphafullv]=extendfixedptalp(alphav,alpha0v,M0,c1,t1v,nt0,...
    Phibarst,DMz,funcparams)
%
% Extend the alpha function to cover all of [t_0,\Tbar].
%
%
% functions called: etainverse
%
%
%
% Extract data.
%
M0p=M0';
%
nfullp1=size(t1v,2);
nfull=nfullp1-1;
ntp1=nt0+1;
ntp2=ntp1+1;
del=t1v(1,2)-t1v(1,1);
%
% CHeck that we actually need to extend.
%
if ntp1==nfullp1
  'There was no need to extend alpha vector in extendalpha.m'
  alphafullv=alphav;
  return;
end;
%
% We're only going to add on to the end of alphav.
%
alphafullv=alphav;
%
% Get the over-all-time-part of the derivative.
% This is a "time-by-time-dimensioned" matrix.
%
% First, we need to do the integral that's on the RHS.
% Integral could be done only one time to save caomputation, I think.
%
mc1delM0=-c1*del*M0;
intforalp=zeros(1,nfullp1);
for is=ntp1:nfullp1
  if ntp1<3
    intforalp(1,is)=0.5*mc1delM0*(DMz(:,is,1)*alphav(1,1)+...
        DMz(:,is,2)*alphav(1,2));
  else
    intforalp(1,is)=0.5*mc1delM0*(DMz(:,is,1)*alphav(1,1)+...
        DMz(:,is,ntp1)*alphav(1,ntp1));
    for isig=2:nt0
      intforalp(1,is)=intforalp(1,is)+mc1delM0*DMz(:,is,isig)*...
          alphav(1,isig);
    end;
  end;  % end of if ntp1<3
end; % end of for is=...
%
% Get RHS.
%
rhsv=alpha0v+intforalp;
%
% Use only extended part.
%
rhsvshort=rhsv(1,ntp2:nfullp1);
%
% Initialize f-p loop in eta^{-1} computation.
%
nnewalp=size(rhsvshort,2);
alphastartv=zeros(1,nnewalp);
%
% Obtain alpha from eta^{-1} of its argument.
%
[newalphav]=etainverse(rhsvshort,alphastartv,funcparams,c1);
%
% Concatenate them.
%
%'in extend...we have alphafullv given by'
%alphav
alphafullv=[alphav,newalphav];
%
% Check length
%
if size(alphafullv)~=size(t1v)
    'Error in sizing of extended alpha in extendfixedptalp.m'
    'Pausing.'
    pause;
end;
