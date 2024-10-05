function [gradthetav]=gradtheta(alphav,funcparams)
gradthetav = funcparams{2}(alphav);
end