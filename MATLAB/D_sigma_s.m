function [Dssig] = D_sigma_s(t, Phi)
% Computes \breve{\mathcal{D}}_{\tau, sigma}
% Input:
%   t -- grid points for time
%   \Phi -- the state transition matrices (4-D tensor)
% Output:
%   Dssig -- \breve{\mathcal{D}}_{\tau, sigma} as defined in Section 9.1
dim=size(Phi,1);
len=length(t);
Dssig=zeros(dim,dim,len,len);
for is=1:len
    for isig=1:len
        bound = min(is, isig);
        Dssig(:, :, is, isig) = squeeze(trapz(t(1:bound), pagemtimes(Phi(:,:,1:bound,is), "transpose", Phi(:,:,1:bound,isig), "none"), 3));
    end
end
