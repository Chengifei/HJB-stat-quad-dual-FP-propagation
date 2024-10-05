function [backsubchkW,relbschkW,Wtarray,gradvarray]=backsubchkofWndim(Warray,...
    Abarv,Pv,Phibarst,DMz,nt1,...
    ycenter,icomp1,icomp2,x1v,x2v,inoncompv,ynoncompv,tfullv,...
    C,Amat,M0,c1,Gammatilde,funcparams)
%
%
% Does a backsubstitution check into the original HJ PDE.
%
% input variables:
%   Warray - W(t,x1,x2) on tfullv,x1v,x2v.
%   tfullv,x1v,x2v - associated domain arrays
%   funcparams - the parameters ([ktemp,eps]).
%   Others are parameters; see paper.
%
% Extract data
% (Get \sigma\sigma' from $\Gammatilde=\sigma\sigma'+c_2\Lzero\Lzerop$.)
%
M0p=M0';
dim=size(M0,2);
sigsigp=Gammatilde; % NOTE: Uses L^0=0!
Sc1=C-c1*M0p*M0;
invc1=1/c1;
ntfp1=size(tfullv,2);
nx1p1=size(x1v,2);
nx2p1=size(x2v,2);
%
nnoncomp=size(inoncompv,2);
%
dt=tfullv(1,2)-tfullv(1,1);
dx1=x1v(1,2)-x1v(1,1);
dx2=x2v(1,2)-x2v(1,1);
idt=1/dt;
idx1=1/dx1;
idx2=1/dx2;
i2dt=idt/2;
i2dx1=idx1/2;
i2dx2=idx2/2;
dx3=min(dx1,dx2);
idx3=1/dx3;
i2dx3=idx3/2;
%
% Get Ntilde.
%
% First we need to find the values of y_1.
%
if (icomp1~=1)&&(icomp2~=1)
  numi1=inoncompv==1;
  if size(numi1)>1
    'Error 1 in backsubchkofW... in index location. Pausing.'
    pause;
  elseif sum(numi1)==0
    'Error 2 in backsubchkofW... in index location. Pausing.'
    pause;
  else
    locof1=find(numi1);
  end;
  y1v=ynoncompv(1,locof1)*ones(size(x1v));
else
  %
  % The vector of y1 values (could be only one element, or entire vector of
  % such.
  %
  y1v=ycenter(1,1)+x1v;
end;
%['y1v= ',num2str(y1v)]
%'Pause after y1v in backsubchkofW...'
%pause;
%
% Get \tilde N(y_1) on y1v
%
[Ntildev,alphavfromy1v,thetafromy1v]=Ntildeviafxdpt(y1v,c1,funcparams);
%
%
% Loop over time and space variables.
%
allindices=cumsum(ones(1,dim));
sumterms=zeros(ntfp1,nx1p1,nx2p1);
Wtarray=zeros(ntfp1,nx1p1,nx2p1);
gradvarray=zeros(dim,ntfp1,nx1p1,nx2p1);
%Wx1a=zeros(ntfp1,nx1p1,nx2p1);
%Wx2a=zeros(ntfp1,nx1p1,nx2p1);
backsubchkW=zeros(ntfp1,nx1p1,nx2p1);
relbschkW=zeros(ntfp1,nx1p1,nx2p1);
%Wp1mat=zeros(dim,ntfp1);
%Wm1mat=zeros(dim,ntfp1);
Wx3v=zeros(dim,ntfp1);
'Would you like the order \delta method or the order \delta^2 method?'
'The first take M times as long as the original computation,'
'while the latter takes 2M times as long,'
'where M is the number of components of the state minus 2.'
prompt1="Would you the order delta^2 method? (y/n) ";
txtq1=input(prompt1,"s")
if txtq1=='y'
  orderdeltasquared=1;
else
  orderdeltasquared=0;
end
for ix1=1:nx1p1
  x1=x1v(1,ix1);
  for ix2=1:nx2p1
    x2=x2v(1,ix2);
    %
    % We need W off center in the other arguments for compuation of
    % partials.
    %
    for ii=1:dim
      %
      % The case where the 1st component one of the plotting axes.
      %
      if (ii~=icomp1)&&(ii~=icomp2)
        %
        % Get the components in the plotting grid.
        %
        y=ycenter;
        y(icomp1,1)=y(icomp1,1)+x1;
        y(icomp2,1)=y(icomp2,1)+x2;
        % 
        % Get the other components.
        %
        for iic=1:nnoncomp
           y(inoncompv(1,iic),1)=ynoncompv(1,iic);
        end;
        %
        % Compute the partials.
        %
        ycp1=y;
        ycp1(ii,1)=y(ii,1)+dx3;
        [Watycp1v,relbsalpproperr]=Watoneypoint(ycp1,tfullv,M0,Abarv,...
            Pv,Phibarst,DMz,...
            ycp1(1),c1,Gammatilde,nt1,funcparams);
        %Wp1mat(ii,:)=Watycp1v;
        if orderdeltasquared==0
          Wx3v(ii,:)=idx3*(Watycp1v-(Warray(:,ix1,ix2))');
        else
          ycm1=y;
          ycm1(ii,1)=y(ii,1)-dx3;
          ycm11=ycm1(1,1);
          [Watycm1v,relbsalpproperr]=Watoneypoint(ycm1,tfullv,M0,Abarv,...
            Pv,Phibarst,DMz,...
            ycm11,c1,Gammatilde,nt1,funcparams);
          Wx3v(ii,:)=i2dx3*((Watycp1v-Watycm1v)');
        end;
      end;
    end;
    for it=ntfp1:-1:1
      %
      % Get partials for those arguments over which there's a grid.
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
      yvv=y;
      gradv=Wx3v(:,it);
      gradv(icomp1,1)=Wx1;
      gradv(icomp2,1)=Wx2;
      gradvarray(:,it,ix1,ix2)=gradv(:,1);
      %
      % Do the backsubstitution check.
      %
      xqterm=0.5*yvv'*C*yvv;
      bqterm=gradv'*Amat*yvv;
      pqterm=0.5*gradv'*sigsigp*gradv;
      Ntildeatx=Ntildev(1,ix1); % Only depends on y_1
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
    end
  end
end
relbschkW=backsubchkW./sumterms;