function [thetav]=gettheta(alphav,funcparams)
% Evaluates \Theta at alphav
    thetav = funcparams{1}(alphav);
end
