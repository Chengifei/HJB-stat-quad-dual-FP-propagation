% tetainverse
%
% Tests etainverse
%
c1=3; % coefficient related to stat-quad duality
c1=-3; % coefficient related to stat-quad duality
c1=-2.1; % coefficient related to stat-quad duality
invc1=1/c1;
%
% parameters that define the nonllinear function.
%
%ktemp=2;
ktemp=-2;
eps=1;
%eps=0.5;
funcparams=[ktemp,eps];
%
% eta data
%
neta=10;
etalo=-5;
etahi=5;
etav=linspace(etalo,etahi,neta)
%
% starting values for the fixed-point method in etainverse.
%
alphastartv=zeros(1,neta);
%
[alphav]=etainverse(etav,alphastartv,funcparams,c1)
%
% Check these results.
%
[gradthetav]=gradtheta(alphav,funcparams);
chkv=etav-(alphav-invc1*gradthetav)
