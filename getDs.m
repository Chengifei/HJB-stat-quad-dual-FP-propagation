function [Dssig,Dpsigs]=getDs(t,Phibarst)
%
% Obtain Dssig and DPsigs,
% i.e., \int_0^s \Phibar_{\tau,s)' \Phibar_{\tau,\sigma}\,d\tau
% and \int_0^\sigma \Phibar_{\tau,s)' \Phibar_{\tau,\sigma}\,d\tau.
%
% input variables:
%   dim - space dimension.
%   t1v - vector on time steps.
%
dim=size(Phibarst,1);
len=length(t);
Dssig=zeros(dim,dim,len,len);
Dpsigs=zeros(dim,dim,len,len);
for is=1:len
    for isig=1:len
         Dssig(:, :, is, isig) = squeeze(trapz(t(1:is), pagemtimes(Phibarst(:,:,1:is,is), "transpose", Phibarst(:,:,1:is,isig), "none"), 3));
         Dpsigs(:, :, isig, is) = squeeze(trapz(t(1:isig), pagemtimes(Phibarst(:,:,1:isig,is), "transpose", Phibarst(:,:,1:isig,isig), "none"), 3));
    end
end
