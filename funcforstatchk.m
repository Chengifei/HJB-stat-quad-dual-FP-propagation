function [f]=funcforstatchk(a);
%
% Just a function evaluator; looking at stat-quad duals to
% find good nonlinearities for examples.
%
eps=1.0;
ktemp=-2;
f=ktemp*sqrt(eps+a^2);