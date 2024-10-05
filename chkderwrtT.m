function [swrelnormsdererr,relnormdererr]=chkderwrtT(alphav,dalpdtv,...
    Phibarstp,Dssig,Dpsigs,t1v,x,c1,Abarv,M0,funcparams)
%
% Checks the computation of derivative of \alpha^T_s wrt T.
%
%
% Recover data.
%
dim=size(x,1);
ntp1=size(alphav,2);
del=t1v(1,2)-t1v(1,1);
%
% ntp1=istep, I think.
%
Abarpv=zeros(dim,dim,ntp1);
for is=1:ntp1
    Abarpv(:,:,is)=Abarv(:,:,is)';
end;
%
% The - M^0 \Phibar'_{T,s} \Abar' x term.
%
alpzerochkv=zeros(1,ntp1);
for is=1:ntp1
  alpzerochkv(1,is)=-M0*Phibarstp(:,:,ntp1,is)*Abarpv(:,:,ntp1)*x;
end;
%
% \Dbar_{\sigma,s} def from of 7/6/23, page 3 multiplied by (M^0)'
%
M0p=M0';
DMz=zeros(dim,ntp1,ntp1);
for ik=1:ntp1
   for ij=1:ntp1
    if ij<ik
      DMz(:,ik,ij)=squeeze(Dpsigs(:,:,ij,ik)*M0p);
    else
      DMz(:,ik,ij)=squeeze(Dssig(:,:,ik,ij)*M0p);
    end
   end
end
%
% The c_1 M^0 \Dbar_{T,s} M^0^' \alpha^T_T term
%
alptermtermv=zeros(1,ntp1);
for is=1:ntp1
    %alptermtermv(1,is)=c1*M0*DMz(:,ntp1,is)*alphav(:,ntp1);
    %alptermtermv(1,is)=0.5*c1*M0*(DMz(:,ntp1,is)*alphav(:,ntp1)+...
    %    DMz(:,ntp1-1,is)*alphav(:,ntp1-1));
    alptermtermv(1,is)=0.5*c1*M0*(DMz(:,is,ntp1)*alphav(:,ntp1)+...
        DMz(:,is,ntp1-1)*alphav(:,ntp1-1));
end;
%
%    DMz(:,ntp1,ntp1-1)
%    squeeze(Dpsigs(:,:,ntp1,ntp1-1)*M0p)
%    DMz(:,ntp1,ntp1)
%    [alphav(:,ntp1-1),alphav(:,ntp1)]
%
%
% The
% c_1 M^0\int_t^T \Dbar_{\sigma,s} M^0' \frac{d\alpha^T_\sigma}{dT} \d\sigma
% term.
% Via rectangle rule for now (only a check).
%
nt=ntp1-1;
delconeM0=del*c1*M0;
inttermv=zeros(1,ntp1);
for is=1:ntp1
  inttermv(1,is)=0.5*delconeM0*(DMz(:,is,1)*dalpdtv(:,1)+...
      DMz(:,is,ntp1)*dalpdtv(:,ntp1));
  for isig=2:nt
    inttermv(1,is)=inttermv(1,is)+delconeM0*DMz(:,is,isig)*dalpdtv(:,isig);
  end;
end;
%
% TEMP CHECK STUFF
%
%Amattemp=zeros(ntp1,ntp1);
%for is=1:ntp1
%  Amattemp(is,1)=0.5*delconeM0*DMz(:,is,1);
%  Amattemp(is,ntp1)=0.5*delconeM0*DMz(:,is,ntp1);
%  for isig=2:ntp1-1
%    Amattemp(is,isig)=delconeM0*DMz(:,is,isig);
%  end;
%end;
%Amattemp(1,:)
%Amattemp(2,:)
%Amattemp(ntp1,:)
%
% End TEMP CHK STUFF
%
%
% The \grad(\alpha^T_s)\frac{d\alpha^T_s}{dT} term.
%
[gradetav]=gradeta(alphav,c1,funcparams);
%gradetaterm=zeros(1,ntp1);
%for is=1:ntp1
%    gradetaterm(:,is)=gradetav(:,is)*dalpdtv(:,is)
%end;
gradetaterm=gradetav.*dalpdtv;
%
% The should-be-zero term
%
%shouldbezerochkv=zeros(1,ntp1);
shouldbezerochkv=alpzerochkv+alptermtermv+inttermv+gradetaterm;
%
% Write and/or plot.
%
% Time-step-wise norms.
%
stepwisenormsbzc=abs(shouldbezerochkv);
stepwisesumnorms=abs(alpzerochkv)+abs(alptermtermv)+...
    abs(inttermv)+abs(gradetaterm);
swrelnormsdererr=stepwisenormsbzc(1,:)./stepwisesumnorms(1,:);
%
% Show terms
%
%[gradetaterm;inttermv;alpzerochkv;alptermtermv]
%
% Overall realtive norm error.
%
normsbzc=norm(shouldbezerochkv,2);
sumnorms=norm(alpzerochkv,2)+norm(alptermtermv,2)+...
    norm(inttermv,2)+norm(gradetaterm,2);
relnormdererr=normsbzc/sumnorms;
