function [sol,res,status]=newtonforeta(omegav,c1,...
    funcparams)
%
% Using Newton's method to solve
% F(alpha;k,eps,om)=(alpha-om)+(k1/c1)*alpha/sqrt(eps+alpha^2),
% which is \eta(\alpha)=om with \eta=\alpha-(1/c_1)\grad\theta(\alpha)
% for \theta=-k1*sqrt(eps+alpha^2)
%
% This version works for vectors of inputs.
%
% input variables:
%   omegav - a vector of parameters (see above).
%   c1 - c_1.
%   funcparams - [k,eps] parameters defining \theta function.
%
% outputs:
%   sol - the argstat
%   res - residuals
%   status - 0 for fail 1 for success
grad_theta = funcparams{2};
sol = zeros(size(omegav));
res = zeros(size(omegav));
for n = 1:length(omegav)
    [sol(n), res(n), status_] = fzero(@(alpha) c1 * alpha - grad_theta(alpha) - c1 * omegav(n), omegav(n));
    if status_ < 0
        status = 0;
        return;
    end
end
status = 1;
end
