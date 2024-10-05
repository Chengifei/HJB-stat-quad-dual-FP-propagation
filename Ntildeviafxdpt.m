function [Ntildev,alphav,thetav]=Ntildeviafxdpt(x,c1,funcparams)
%
% Computes Ntilde(x) for our particular \Theta function.
% Here, Ntilde(x) \doteq 
%     \stat_\alpha \{ \Theta(\alpha) -(c_1/2)(\alpha-x)^2 \}
% Obtains the argstat by using fixed point to solve
%  \alpha = (1/c_1) [ x-\grad\Theta(\alpha) ].
%
% FIX?
%  \alpha = (1/c_1) [ x+\grad\Theta(\alpha) ].
%
% Note: This code is nearly identical to that of etainverse, from which
% most of it has been adapted.
%
% Note: Convergence criterial for the f-p alg are hardwired.
%
%
% input variables:
%   x1v - array of x_1 values.
%   c1 - c_1.
%   funcparams - theta function parameters ([ktemp,eps]).
%
% output variables:
%   Ntildev - vector of Ntilde values.
%
% functions called:
%   gettheta, gradtheta, gradeta
%
%
% Extract data.
%
x1v = x(1, :);
[alphav] = etainverse(x1v, zeros(size(x1v)), funcparams, c1);
%
% Get the stat from the argstat.
%
[thetav]=gettheta(alphav,funcparams);
c1o2=c1/2;
%Ntildev=thetav+(c1/2)*(x1v-alphav).^2;
%
% FIX?
%
Ntildev=thetav-(c1/2)*(x1v-alphav).^2;