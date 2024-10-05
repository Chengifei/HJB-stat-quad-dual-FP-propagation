function [effectofhalving,shouldbeahalf]=checkstat(alphav,x,t1v,Phibarst,...
    M0,c1,Gammatilde,funcparams)
%
% CHecks how close to being a stat, an alpha^T_\cdot actually is
% in the stat_{\alpha_\cdot} [x'q^T_t+0.5*r^T_t]
% Do this by varying the alpha vector.
% Doubling the variation should approximately cause a factor of 4 change.
%
% output variables:
%
%   effectofhalving - should be near to 1/4.
%   shouldbeahalf - should be near 1/2 if halved properly.
%
% First get the baseline.
%
[qsbase,rsbase,qbasev,rbasev]=qrintegrals(t1v,alphav,Phibarst,...
    M0,c1,Gammatilde,funcparams);
xqpdrbase=x'*qsbase+0.5*rsbase;
%
% Now perform variations.
%
nalpp1=size(alphav,2);
nfullp1=size(t1v,2);
nfp1=nfullp1;
if nalpp1~=nfullp1
    'Error in checkstat. Vector sizes. Pausing.'
    pause;
end;
invnfp1=1/nfp1;
del=t1v(1,2)-t1v(1,1);
avalpsize=invnfp1*sum(abs(alphav));
%
% SEt the number of random variations that we want to use in the test.
%
numvar=9;
dxqpdrvarv=zeros(1,9);
variationav=zeros(1,9);
dxqpdrvar2v=zeros(1,9);
variationa2v=zeros(1,9);
%
% Vary alpha, and compute the result.
% Also do this for double the size of the variation.
% NOTE: Then using the doubled as the basling variation, and the regular
% as the halved variation.....
%
for ivar=1:numvar
  variationa=avalpsize*randn(1,nfp1);
  alpvaried=alphav+variationa;
  [qsvar,rsvar,qvarv,rvarv]=qrintegrals(t1v,alpvaried,Phibarst,...
      M0,c1,Gammatilde,funcparams);
  xqpdrvar=x'*qsvar+0.5*rsvar;
  dxqpdrvarv(1,ivar)=xqpdrvar-xqpdrbase;
  variationav(1,ivar)=norm(variationa,2);
  % now double
  variationa2=2*variationa;
  alpvaried2=alphav+variationa2;
  [qsvar2,rsvar2,qvar2v,rvar2v]=qrintegrals(t1v,alpvaried2,Phibarst,...
      M0,c1,Gammatilde,funcparams);
  xqpdrvar2=x'*qsvar2+0.5*rsvar2;
  dxqpdrvar2v(1,ivar)=xqpdrvar2-xqpdrbase;
  variationa2v(1,ivar)=norm(variationa2,2);
end;
%
% If we're at a stat, the resulting variation should be o(variationa^2).
%
normdx=norm(dxqpdrvarv,2);
normdv=norm(variationav,2);
normdx2=norm(dxqpdrvar2v,2);
normdv2=norm(variationa2v,2);
%
effectofhalving=normdx/normdx2;
shouldbeahalf=normdv/normdv2;

