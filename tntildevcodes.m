% tntildevcodes
%
%
% Set x range.
%
widthx1=0.0632;
x1min=-widthx1;
x1max=widthx1;
nx1=20;
x1v=linspace(x1min,x1max,nx1);
%
% Parameters.
%
c1=2; % coefficient related to stat-quad duality
c1=-3; % coefficient related to stat-quad duality
%ktemp=2;
ktemp=-2;
eps=1;
%eps=0.5;
funcparams=[ktemp,eps];
%
% Call algorithms.
%
[Ntildevold]=Ntildecomp(x1v,c1,funcparams);
%
[Ntildev]=Ntildeviafxdpt(x1v,c1,funcparams);
%
'Ntildevold:'
Ntildevold
'Ntildev:'
Ntildev
%
Ntildev-Ntildevold
