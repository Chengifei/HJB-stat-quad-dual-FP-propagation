% lxlpalpnlpsiTbndim
%
% For n-dimensional case with one-dimensional nonlinearity in the cost.
% It si assumed that a rotation has been performed such that the
% nonlinearity is in the first component.
% Loops over two dimensions of y outside of
% a loop over time via propalphanl.
%
% functions called:
%   rk4revricandstmnstep[rk4revricandstm[revricandstmfunc]],
%   Chkpsi,
%   propalpnlpsiTbar[GetalpPhinlpsiTbar[
%                      getDs,Getalphafuncnl[
%                                           gradtheta,newtonforeta
%                                          ]
%                      extendfixedptalp[etainverse]
%                     ],
%                   gradeta,
%                   backsubsalpnlTbar[gradtheta],
%                   getdalpTdTnlpsiTbar,
%                   GetPhialp0psiTbar,
%                   getDs
%                  ],
%   qrintegrals[gettheta],
%   Checkqnadr[gettheta],
%   backsubchkofW[Ntildeviafxdpt[
%                                gettheta,gradtheta,gradeta]
%   
%
%
%
% Set up problem data.
%
Amat=[0 1; -1 0]; % defines the linear part of the nominal dynamics
Amat=[0 1.5; -1 0]; % defines the linear part of the nominal dynamics
Amat=[0 1.5 0; -1 0 0; 0.5 0 -0.5]
%Amat=[0 0 1.5; 0.5 0 -0.5; -1 0 0] %permutation of indexing of previous one
%Amat=[0 1.5 0; 0.5 0 -0.5; -1 0 0]  % Not a permuation of first one
%Amat=[0 1 0.05 0; -1 0 -0.1 -0.1; 0.0 0 0.1 0.25; 0 -0.1 0 -0.25]
%Amat=[0 1 0 0; -1 0 0.1 0; 0.5 0 0 -0.25; 0 0 1 0]
%Amat=[0 1 0 0; -0.25 0 0 0.5; 0 0.1 0 -1; 0 0 1 0];
Amatp=Amat';
dim=size(Amat,1);
M0=[1, zeros(1,dim-1)]; % Sets y_1 dim to be domain/range of  nonlinearity.
M0p=M0';
Gammatilde=eye(dim,dim);dim
Cmat=0.6*eye(dim,dim);
c1=2; % coefficient related to stat-quad duality
%c1=-1.5; % coefficient related to stat-quad duality
%c1=4; % coefficient related to stat-quad duality
c1=2.5; % coefficient related to stat-quad duality
c1=2.75; % coefficient related to stat-quad duality
c1=-3; % coefficient related to stat-quad duality
%c1=-2.1; % coefficient related to stat-quad duality
%c1=-2.75; % coefficient related to stat-quad duality
%
% parameters that define the nonllinear function.
%
%ktemp=2;
ktemp=-2;
eps=1;
%eps=1;
eps=0.75;
funcparams=[ktemp,eps];
%
% In this code, the original space is in "y" variables.
% We'll use "x_1,x_2" to denote the specific components over which we plot,
% regardless of which components of y these are.
%
% Set the grid to be used for plotting, and the y components that are
% constant.
% The grid will be over the components indicated by indices icomp1 and
% icomp2.
% The grid will be centered at ycenter, with the other components (not
% icomp1 and icomp2) being fixed at the "center" points.
%
%
icomp1=1;
icomp2=2;
inoncompv=[3]; % Vector indices that are not in the grid.
%
%icomp1=1;
%icomp2=3;
%inoncompv=[2]; % Vector indices that are not in the grid.
%
%icomp1=2;
%icomp2=3;
%inoncompv=[1]; % Vector indices that are not in the grid.
%
nnoncomp=size(inoncompv,2);
%
ycenter=zeros(dim,1);
%ycenter=[0.5;0;0;0];
%widthx1=0.0632;
%widthx1=0.2;
widthx1=0.4;
widthx1=3;
widthx1=5;
%widthx1=3;
x1min=-widthx1;
x1max=widthx1;
%x1max=-widthx1/2;
%x1max=0;
%
widthx2=0.2;
widthx2=0.4;
widthx2=2;
widthx2=4;
x2min=-widthx2;
x2max=widthx2;
%x2max=-widthx2/2;
%x2max=0;
%
nx1=20;
nx2=20;
%nx1=15;
%nx2=15;
nx1=30;
nx2=30;
nx1=50;
nx2=30;
x1v=linspace(x1min,x1max,nx1);
x2v=linspace(x2min,x2max,nx2);
%
ynoncompv=[-3]; % components of y that are not in the grid
if size(ynoncompv,2)~=nnoncomp
    'Error. Wrong ynoncompv size. Pausing.'
    pause;
end;
%
% We'll be looking at t0\le t \le s,T \le \Tbar
% where t is initial and T is (moving) terminal time.
% We expect that T will vary between t and \Tbar, as will s.
% We also expect that when T=\Tbar, we will have the solution
% we're actually looking for with \alpha^T_s for s\in [t,T=\Tbar].
%
t0=0;       % initial time (generally zero). This is "t" in the paper.
Tbar=0.5;
%
% We use the same time discretization over all of [t,\Tbar].
%
% Set the number of steps over [t0,\Tbar].
%
ntot=20;
ntot=40;
%ntot=60;
ntotp1=ntot+1;
%
% Set the time step at which we switch from fixed-point to propagation.
% At step nt1 (i.e., at index ntp1=nt1+1, the fixed-point method is used.
% This is at time t10=t0+dt*ntp1.)
% For all time steps, it, such that it>ntp1, the ODE starting from 
% initial time t10=t0+dt*ntp1 is used.
% That is, the value at any INDEX it>ntp1, the value has been obtained
% by the ODE propagation.
%
%
nt1=5;
nt1=10;
%nt1=15;
ntp1=nt1+1;
% 
% The number of propagations steps (post-fixed-pt sol region).
%
numstep=ntot-nt1;
%
% Get the full time period.
%
tfullv=linspace(t0,Tbar,ntotp1);
%
% Get time-step size.
%
dt=tfullv(1,2)-tfullv(1,1);
%
% Set terminal conditions.
%
Pend=zeros(dim,dim);
Phibarend=eye(dim,dim);
Pstmend=[Pend;Phibarend];
%
% Solve the DREs over the full time period.
%
Sc1=Cmat-c1*M0p*M0;
Amatp=Amat';
[t1v,Pstmv]=rk4revricandstmnstep(Tbar,Pstmend,t0,ntot,...
    Sc1,Amat,Amatp,Gammatilde);
Pv=Pstmv(1:dim,:,:);
Psiv=Pstmv(dim+1:2*dim,:,:);
%
% A check.
%
if size(tfullv)~=size(t1v)
    'Error in size of t1v vs tfullv in lxl... Pausing.'
    pause;
end;
%
% Get Abarv and check Psi's.
%
[Abarvchk,Psivchk,relPsivchk]=ChkPsi(Psiv,Pv,t1v,Gammatilde,Amatp);
Psivchk
relPsivchk
'That was Psivchk and relPsichkv in lxl.'
%pause;
%
% Get the Abar function and the transposes.
%
Abarv=zeros(dim,dim,ntotp1);
Abarpv=zeros(dim,dim,ntotp1);
for it=1:ntotp1
    Abarv(:,:,it)=Pv(:,:,it)*Gammatilde-Amatp;
    Abarpv(:,:,it)=Abarv(:,:,it)';
end;
%
% Get $\Phibar_\tau^s$ for all generic \tau, s.
% Also get $(\Phibar_\tau^s)'$.
% NOTE the ordering matches the subscript ordering of state transition
% matrices.
%
% Also, getting these for the entire duration, rather than only up to the
% initial terminal time, t10.
%
nfullp1=size(Psiv,3);
Phibarst=zeros(dim,dim,nfullp1,nfullp1);
Phibarstp=zeros(dim,dim,nfullp1,nfullp1);
Psiinvv=zeros(dim,dim,nfullp1);
for it=1:nfullp1
    Psiinvv(:,:,it)=inv(Psiv(:,:,it));
end
for it0=1:nfullp1
  for it1=1:nfullp1
      Phibarst(:,:,it1,it0)=Psiv(:,:,it1)*Psiinvv(:,:,it0);
      Phibarstp(:,:,it1,it0)=(Phibarst(:,:,it1,it0))';
  end
end
%
% Get Dssig and DPsigs
%
[Dssig,Dpsigs]=getDs(dim,t1v,Phibarst,Phibarstp);
%
% Get the over-all-time-part of the derivative.
% This is a "time-by-time-dimensioned" matrix.
%
DMz=zeros(dim,nfullp1,nfullp1);
for ik=1:nfullp1
 for ij=1:nfullp1
  if ij<ik
    DMz(:,ik,ij)=squeeze(Dpsigs(:,:,ij,ik)*M0p);
  else
    DMz(:,ik,ij)=squeeze(Dssig(:,:,ik,ij)*M0p);
  end;
 end;
end;
%
%chkAbar=Abarvchk-Abarv;
%for it=1:ntotp1
%    norm(chkAbar(:,:,it),'fro')
%end
%'Pausing after Abar chk in lxl.'
%pause;
%
% Loop over designated components.
%
Warray=zeros(ntotp1,nx1,nx2);
rarray=zeros(ntotp1,nx1,nx2);
qvarray=zeros(dim,ntotp1,nx1,nx2);
alparray=zeros(nx1,nx2);
alpvarray=zeros(ntotp1,nx1,nx2);
relbacksubserrpropmat=zeros(nx1,nx2);
maxrelbsalpproperr=0;
%
% If we want to check that the staticizing alphas we compute are truly
% stat, set checkingstat=1.
%
checkingstat=0;
if checkingstat==1
  statchksv=[];
end;
%
% If we are going to double-check W_t.
%
%chkWts=1;
%if chkWts==1
%  Wtchkarray=zeros(ntotp1,nx1,nx2);
%end;
%
% Now really beginning loop over designated components.
%
for ix1=1:nx1
  for ix2=1:nx2
    y=ycenter;
    y(icomp1,1)=ycenter(icomp1,1)+x1v(1,ix1);
    y(icomp2,1)=ycenter(icomp2,1)+x2v(1,ix2);
    for iic=1:nnoncomp
      y(inoncompv(1,iic),1)=ynoncompv(1,iic);
    end;;
    ['y= ',num2str(y')]
    y1=y(1,1);
    %
    % Get alpha0v (function of time, dependent on y).
    %
    [alpha0v]=Getalpha0psiTbar(M0,y,t1v,Phibarstp);
    %
    % Call propalpnlpsiTbar.
    %
    [allalpham,t1vfinal,Phibarst,relbsalpproperr]=...
        propalpnlpsiTbndim(M0,Abarv,Abarpv,Phibarst,Phibarstp,...
        Dssig,Dpsigs,DMz,alpha0v,...
        y1,c1,t0,Tbar,nt1,numstep,funcparams);
    numalpfuncs=size(allalpham,1);
    alphafinalv=allalpham(numalpfuncs,:);
    maxrelbsalpproperr=max(maxrelbsalpproperr,relbsalpproperr);
    relbacksubserrpropmat(ix1,ix2)=relbsalpproperr;
    %['The running maximum-over-y post-propagation rel alpha backsubs error = ',...
    %    num2str(maxrelbsalpproperr)]
    %if abs(maxrelbsalpproperr-relbsalpproperr)<1e-8
    %    ['Hit the running-maximum-over-y error at ',num2str(y')]
    %    %pause;
    %end;
    if maxrelbsalpproperr>0.45
        'Max rel backsub err from propalp...Tbar >0.45. Pausing in lxl...'
        pause;
    end;
    %
    % Check whether this is really stat or not, if desired.
    %
    if checkingstat==1
      alphav4chk=allalpham(numalpfuncs,:)
      [effectofhalving,shouldbeahalf]=checkstat(alphav4chk,y,t1v,Phibarst,...
      M0,c1,Gammatilde,funcparams);
      statchksv=[statchksv,effectofhalving];
      %['Checking if alpha really stat at y= ', num2str(y')]
      %['The effect of halving (should be near 0.25) is ',num2str(effectofhalving)]
      %['Check on halving (should be near 0.5) is ',num2str(shouldbeahalf)]
      'Pausing in lxl...'
      %pause;
    end;
    %
    %
    % Get q_s and r_s for the staticizing alpha
    %
    [qs,rs,qv,rv]=qrintegrals(t1vfinal,alphafinalv,Phibarst,...
        M0,c1,Gammatilde,funcparams);
    %
    % Compute W(t1v(is),y)=0.5*y'*Pv(:,:,is)*y+y'*qv(1,is)+0.5*rv(1,is)
    %
    Wv=zeros(1,ntotp1);
    for it=1:ntotp1
      Wv(1,it)=0.5*y'*Pv(:,:,it)*y+y'*qv(:,it)+0.5*rv(1,it);
      qvarray(:,it,ix1,ix2)=qv(:,it);
    end;
    %
    % Put result into big W array.
    %
    Warray(:,ix1,ix2)=Wv(1,:)';
    rarray(:,ix1,ix2)=rv(1,:)';
    alparray(ix1,ix2)=allalpham(numalpfuncs,1);
    alpvarray(:,ix1,ix2)=(allalpham(numalpfuncs,:))';
    %
    % Check qv,rv if desired.
    %
    chkqandr=0;
    if chkqandr==1
      [qchkv,rchkv,relqchkv,relrchkv]=Checkqandr(qv,rv,t1vfinal,Pv,alphafinalv,...
        Amat,Cmat,M0,Gammatilde,c1,Phibarst,Abarv,funcparams);
      'Checking q and r. Pausing.'
      relqchkv
      relrchkv      
      pause;
    end;
    %
    % Check that W_t is what we expect.
    %
    %if chkWts==1
    %  Amatp=Amat';
    %  alphav4chk=allalpham(numalpfuncs,:);
    %  [Thetav]=gettheta(alphav4chk,funcparams);
    %  for it=1:ntotp1
    %    Pmat=Pv(:,:,it);
    %    Pdotrhs=-(Sc1+Amatp*Pmat+Pmat*Amat-Pmat*Gammatilde*Pmat);
    %    qdotrhsalt=Abarv(:,:,it)*qv(:,it)-c1*M0p*alphav4chk(1,it);
    %    rdotrhs=qv(:,it)'*Gammatilde*qv(:,it)+c1*alphav4chk(1,it)^2-2*Thetav(1,it);
    %    Wtchkarray(it,ix1,ix2)=0.5*(y'*Pdotrhs*y+rdotrhs)+y'*qdotrhsalt;
    %  end;
    %end;
  end % loop over ix2
end % loop over ix1
%
% check of stats if it was desired.
%
if checkingstat==1
    'If the alphas were truly stat, the entries in the vector to follow'
    ' should be approximately 0.25'
    statchksv
    'Pausing in lxl...'
    pause;
end;
%
% PLot W
%
figure;
surf(x2v,x1v,squeeze(Warray(1,:,:)));
title('Value function');
ylabel(['y_',num2str(icomp1)])
xlabel(['y_',num2str(icomp2)])
%
% Plot some alpha values at some (which?!) time.
%
%figure;
%surf(x2v,x1v,squeeze(alparray(:,:)));
%title('Alphas')
%ylabel('x_1');
%xlabel('x_2');
%
% Plot backsubstitution errors from the alpha fxd-pt and ODE prop
% computations.
%
figure;
surf(x2v,x1v,squeeze(relbacksubserrpropmat(:,:)));
title('Relative back-substitution error post-prop')
ylabel(['y_',num2str(icomp1)])
xlabel(['y_',num2str(icomp2)])
%
% Backsubstitution check in original PDE if desired.
%
'You have the option of doing a backsubstitution check into the original PDE.'
'Note that this will take M (order \delta method)'
'or 2M (order \delta^2 method) times as long as the original computation,'
'where M is the number of components of the state minus 2.'
prompt1="Would you like check by backsubstitution into the PDE? (y/n) ";
txtq1=input(prompt1,"s")
if txtq1=='y'
%dobacksubspdechk=0;
%if dobacksubspdechk==1
  [backsubchkW,relbschkW,Wtarray,gradvibsubarray]=backsubchkofWndim(Warray,...
    Abarv,Abarpv,Pv,Phibarst,Phibarstp,Dssig,Dpsigs,DMz,nt1,numstep,...
    ycenter,icomp1,icomp2,x1v,x2v,inoncompv,ynoncompv,tfullv,...
    Cmat,Amat,M0,c1,Gammatilde,funcparams);
  figure;
  surf(x2v,x1v,squeeze(relbschkW(1,:,:)));
  title('Relative backsubstitution error in original PDE');
  ylabel(['y_',num2str(icomp1)])
  xlabel(['y_',num2str(icomp2)])
  %
  % Plot the backsubstitution check of the original PDE.
  %
  %
  % Wt check...
  %
  %if chkWts==1
  %  delWtchkm=Wtchkarray-Wtarray;
  %  reldelWtchk=delWtchkm./(abs(Wtchkarray)+abs(Wtarray));
  %  listall=0;
  %  if listall==1
  %    for ix1=1:nx1
  %      for ix2=1:nx2
  %        ['At ix1,ix2= ',num2str([ix1,ix2])]
  %        ['relative W_t errors are ',num2str((reldelWtchk(2:ntot,ix1,ix2))')]
  %        %pause;
  %      end;
  %    end;
  %  end;
  %end;
  %
  % Plot the backsubstitution check of the original PDE at the end of the
  % fxd-pt step region (i.e., before the ODE propagation).
  %
  %figure;
  %surf(x2v,x1v,squeeze(relbschkW(ntp1,:,:)));
  %title('Relative backsubstitution error in original PDE at start of propagation');
  %ylabel('x_1');
  %xlabel('x_2');
end;
