function [Phist]=getPhiT(s,t)
%
% Generate state transition matrix $\Phi^T_{s,t}$ for 
% Abar^T=[ 0 -1; 1 0]
%
cts=cos(s-t);
sts=sin(s-t);
Phist=[cts -sts; sts cts];