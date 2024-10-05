function [backsubchkW,relbschkW,Wtarray]=backsubchkofW(Warray,x1v,x2v,tfull1v,...
    C,Amat,M0,c1,Gammatilde,funcparams,Pv,qvarray,alpvarray)
%
%
% Does a backsubstitution check into the original HJ PDE.
%
% input variables:
%   Warray - W(t,x1,x2) on t1v,x1v,x2v.
%   tfull1v,x1v,x2v - associated domain arrays
%   funcparams - the parameters ([ktemp,eps]).
%   Others are parameters; see paper.
%
% Extract data
% (Get \sigma\sigma' from $\Gammatilde=\sigma\sigma'+c_2\Lzero\Lzerop$.)
%
M0p=M0';
sigsigp=Gammatilde; % NOTE: Uses L^0=0!
Sc1=C-c1*M0p*M0;
invc1=1/c1;
ntfp1=size(tfull1v,2);
nx1p1=size(x1v,2);
nx2p1=size(x2v,2);
dt=tfull1v(1,2)-tfull1v(1,1);
dx1=x1v(1,2)-x1v(1,1);
dx2=x2v(1,2)-x2v(1,1);
idt=1/dt;
idx1=1/dx1;
idx2=1/dx2;
i2dt=idt/2;
i2dx1=idx1/2;
i2dx2=idx2/2;
%
% Get \tilde N(x_1) on x1v
%
[Ntildev,alphavfromx1v,thetafromx1v]=Ntildeviafxdpt(x1v,c1,funcparams);
% OLD VERSION: [Ntildev]=Ntildecomp(x1v,c1,funcparams);
%
%
% Loop over time and space variables.
%
sumterms=zeros(ntfp1,nx1p1,nx2p1);
Wtarray=zeros(ntfp1,nx1p1,nx2p1);
%Wx1a=zeros(ntfp1,nx1p1,nx2p1);
%Wx2a=zeros(ntfp1,nx1p1,nx2p1);
backsubchkW=zeros(ntfp1,nx1p1,nx2p1);
relbschkW=zeros(ntfp1,nx1p1,nx2p1);
%'FIX TEMP CHANGE IN backsub...W!!!!'
%'pausing'
%pause;
for it=ntfp1:-1:1
  t=tfull1v(1,it);
  for ix1=1:nx1p1
    x1=x1v(1,ix1);
    for ix2=1:nx2p1
      x2=x2v(1,ix2);
      %
      % Get partials.
      %
      if it==1
        Wt=idt*(Warray(it+1,ix1,ix2)-Warray(it,ix1,ix2));
      elseif it==ntfp1
        Wt=idt*(Warray(it,ix1,ix2)-Warray(it-1,ix1,ix2));
      else
        Wt=i2dt*(Warray(it+1,ix1,ix2)-Warray(it-1,ix1,ix2));
      end;
      Wtarray(it,ix1,ix2)=Wt;
      if ix1==1
        Wx1=idx1*(Warray(it,ix1+1,ix2)-Warray(it,ix1,ix2));
      elseif ix1==nx1p1
        Wx1=idx1*(Warray(it,ix1,ix2)-Warray(it,ix1-1,ix2));
      else
        Wx1=i2dx1*(Warray(it,ix1+1,ix2)-Warray(it,ix1-1,ix2));
      end
      if ix2==1
        Wx2=idx2*(Warray(it,ix1,ix2+1)-Warray(it,ix1,ix2));
      elseif ix2==nx2p1
        Wx2=idx2*(Warray(it,ix1,ix2)-Warray(it,ix1,ix2-1));
      else
        Wx2=i2dx2*(Warray(it,ix1,ix2+1)-Warray(it,ix1,ix2-1));
      end
      xvv=[x1;x2];
      gradv=[Wx1;Wx2];
      %
      % A check that W_x is what it should be.
      %
      %if it<ntfp1
      %  gradv
      %  wxchk=Pv(:,:,it)*xvv+qv(:,it)
      %  'pause in backsub...W'
      %  pause;
      %end;
      %
      % Do the backsubstitution check.
      %
      xqterm=0.5*xvv'*C*xvv;
      bqterm=gradv'*Amat*xvv;
      pqterm=0.5*gradv'*sigsigp*gradv;
      Ntildeatx=Ntildev(1,ix1); % Only depends on X_1
      backsubchkW(it,ix1,ix2)=Wt+xqterm+bqterm+Ntildeatx-pqterm;
      sumterms(it,ix1,ix2)=abs(Wt)+abs(xqterm)+abs(bqterm)+...
          abs(pqterm)+abs(Ntildeatx);
      %
      % Output some checks if desired.
      %
      printnumerical=0;
      if (it==1)&&(printnumerical==1)
        [Wt,xqterm,bqterm,Ntildeatx,pqterm,gradv']
        backsubchkW(it,ix1,ix2)
        %pause;
      end
      %
      % Check analytical if desired.
      %
      % NEEDS Wtchkarray as an input!
      %
      %checkanalytical=0;
      %if checkanalytical==1
      %  x1chk=xvv(1,1);
      %  qvanal=qvarray(:,it,ix1,ix2);
      %  gradWanal=Pv(:,:,it)*xvv+qvanal;
      %  bqanal=gradWanal'*Amat*xvv;
      %  pqanal=0.5*gradWanal'*sigsigp*gradWanal;
      %  Wtanal=Wtchkarray(it,ix1,ix2);
      %  alpchk=alpvarray(it,ix1,ix2);
      %  alpchk2=alphavfromx1v(1,ix1);
      %  [thetachk]=gettheta(alpchk,funcparams);
      %  thetachk2=thetafromx1v(1,ix1);
      %  [gradthetachk]=gradtheta(alpchk,funcparams);
      %  [gradthetachk2]=gradtheta(alpchk2,funcparams);
      %  chkN=invc1*gradthetachk-(alpchk-x1chk);
      %  chkN2=invc1*gradthetachk2-(alpchk2-x1chk);
      %  Ntildechk=thetachk-(c1/2)*(alpchk-xvv(1,1))^2;
      %  backanal(it,ix1,ix2)=Wtanal+xqterm+bqanal+Ntildeatx-pqanal;
      %  sumtermsanal(it,ix1,ix2)=abs(Wtanal)+abs(xqterm)+abs(bqanal)+...
      %    abs(pqanal)+abs(Ntildeatx);
      %  if it==1
      %    [Wtanal,xqterm,bqanal,Ntildeatx,pqanal,gradWanal']
      %    [qvanal',Ntildechk,backanal(it,ix1,ix2)]
      %    [thetachk,thetachk2,alpchk,alpchk2,chkN,chkN2]
      %    %pause;
      %  end;
      %end;
    end
  end
end
relbschkW=backsubchkW./sumterms;