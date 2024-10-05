function [intonrhsv]=intonrhsofetaofalp(alphav,t1v,nt0,c1,M0,DMz)
%
% Hastily written code to check
% -c_1 M^0  \int_t0^T \Dbar_{\sigma,s}(M^0)'\alpha_\sigma\,d\sigma.
% Here DMz=\Dbar_{\sigma,s}(M^0)'
%
ntp1=nt0+1;
nfullp1=size(t1v,2);
del=t1v(1,2)-t1v(1,1);
%
mc1delM0=-c1*del*M0;
intonrhsv=zeros(1,nfullp1);
for is=2:nfullp1
  if ntp1==1
      intonrhsv(1,is)=0;
  elseif ntp1==2
    intonrhsv(1,is)=0.5*mc1delM0*(DMz(:,1,is)*alphav(1,1)+...
        DMz(:,ntp1,is)*alphav(1,ntp1));
  else
    intonrhsv(1,is)=0.5*mc1delM0*(DMz(:,1,is)*alphav(1,1)+...
        DMz(:,ntp1,is)*alphav(1,ntp1));
    for isig=2:nt0
        intonrhsv(1,is)=intonrhsv(1,is)+mc1delM0*DMz(:,isig,is)*alphav(1,isig);
    end;
  end; % end of if ntp1
end; % end of for is=...
