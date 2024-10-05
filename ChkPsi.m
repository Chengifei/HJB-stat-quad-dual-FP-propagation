function [Abarv,Psivchk,relPsivchk]=ChkPsi(Psiv,Pv,tv,Gammatilde,Amatp)
%
% Check the Psi computations, and output Abar's.
%
% functions called:
%   none.
%
% Extract data.
%
ntp1=size(Psiv,3)
ntp1chk=size(Pv,3)
if ntp1chk~=ntp1
    'Error in ChkiPSi! sizes do not match. Pausing.'
    pause;
end;
nt=ntp1-1;
dim=size(Psiv,1);
del=tv(1,2)-tv(1,1);
invdel=1/del;
%
% Set array sizes.
%
Abarv=zeros(dim,dim,ntp1);
Psichkmat=zeros(dim,dim);
Psivchk=zeros(1,ntp1);
Psidotv=zeros(dim,dim,ntp1);
relPsivchk=zeros(1,ntp1);
%
% Get time derivative of Psiv.
%
for it=1:nt
  Psidotv(:,:,it)=invdel*(Psiv(:,:,it+1)-Psiv(:,:,it));
end;
Psidotv(:,:,ntp1)=invdel*(Psiv(:,:,ntp1)-Psiv(:,:,nt));
%
% Get Abarv and Psivchk.
%
for it=ntp1:-1:1
    Pmat=Pv(:,:,it);
    Abarv(:,:,it)=Pmat*Gammatilde-Amatp;
    RHS=Abarv(:,:,it)*Psiv(:,:,it);
    Psichkmat=Psidotv(:,:,it)-RHS;
  Psivchk(1,it)=norm(Psichkmat,'fro');
  relPsivchk(1,it)=Psivchk(1,it)/(norm(Psidotv(:,:,it),'fro')+norm(RHS,'fro'));
end;