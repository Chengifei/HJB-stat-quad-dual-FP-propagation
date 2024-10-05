function [gradetav]=gradeta(alphav,c1,funcparams)
%\eta = \alpha - c_1^{-1} \grad\theta(\alpha)
gradetav = 1 - c1 \ (funcparams{3}(alphav));
end