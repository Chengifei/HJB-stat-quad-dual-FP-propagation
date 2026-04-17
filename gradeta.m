function [gradetav]=gradeta(alphav,c1,funcparams)
% Evaluates gradient of \eta at \alphav,
% with dualizing coefficient c1, and \grad^2\theta specified by funcparams
%\eta = \alpha - c_1^{-1} \grad\theta(\alpha)
    gradetav = 1 - c1 \ (funcparams{3}(alphav));
end
