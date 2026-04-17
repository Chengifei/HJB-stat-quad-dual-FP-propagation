function [gradthetav]=gradtheta(alphav,funcparams)
% Evaluates gradient of \Theta at alphav, as specified by funcparams
gradthetav = funcparams{2}(alphav);
end
