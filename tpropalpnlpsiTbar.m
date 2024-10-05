% tpropalpnlpsiTbar
%
% loop over time via propalphanl, using Abar obtained from solution of the
% DRE
%
% functions called:
%   rk4revricandstmnstep[rk4revricandstm[revricandstmfunc]],
%   Chkpsi,
%   propalpnlpsiTbar[GetalpPhinlpsiTbar[getDs,GetalpnlTbar[
%                                                          gradtheta,newtonforeta
%                                                          ]
%                                       ],
%                  gradeta,
%                  backsubsalpnlTbar[gradtheta],
%                  getdalpTdTnlpsiTbar,
%                  GetPhialp0psiTbar,
%                  getDs
%                  ],
%   qrintegrals[gettheta],
%   Checkqnadr[gettheta],
%   backsubchkofW[Ntildecomp]
%   
%
% Set up problem data.
%
%
% set up problem data
%
M0=[1 0]; % matrix defining lower dim'l set on which nonlinearities live
M0p=M0';
Amat=[0 1; -1 0]; % defines the linear part of the nominal dynamics
Amatp=Amat';
Gammatilde=eye(2,2);
Cmat=0.5*[1 0; 0 1];
Cmat=0.6*[1 0; 0 1];
dim=size(Cmat,1);
c1=2; % coefficient related to stat-quad duality
%c1=-1.5; % coefficient related to stat-quad duality
%c1=4; % coefficient related to stat-quad duality
c1=2.5; % coefficient related to stat-quad duality
c1=2.75; % coefficient related to stat-quad duality
c1=-3; % coefficient related to stat-quad duality
%c1=-2.1; % coefficient related to stat-quad duality
invc1=1/c1;
%
% parameters that define the nonllinear function.
%
%ktemp=2;
ktemp=-2;
eps=1;
%eps=0.5;
funcparams=[ktemp,eps];
%
% Set x 
%
x1=0.029937
x2=0.11579;
x1=0.0370;
x2=0.2;
x=[x1;x2]
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
ntot=80;
ntotp1=ntot+1;
%
% Set time discretization.
%
h=(Tbar-t0)/ntot;
%
% Set the time step at which we switch from fixed-point to propagation.
% At step nt1 (i.e., at index ntp1=nt1+1, the fixed-point method is used.
% This is at time t10=t0+h*ntp1.)
% For all time steps, it, such that it>ntp1, the ODE starting from 
% initial time t10=t0+h*ntp1 is used.
% That is, the value at any INDEX it>ntp1, the value has been obtained
% by the ODE propagation.
%
%
nt1=5;
%nt1=10;
ntp1=nt1+1;
%t10=t0+h*ntp1;
% 
% The number of propagations steps.
%
numstep=ntot-nt1;
%
% Get the full time period.
%
tfullv=linspace(t0,Tbar,ntotp1);
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
% GEt Abarv and check Psi's.
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
Phibarst=zeros(2,2,nfullp1,nfullp1);
Phibarstp=zeros(2,2,nfullp1,nfullp1);
Psiinvv=zeros(2,2,nfullp1);
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
%
%
%
% Begin time propagation.
%
x1=x(1,1);
%
% Get alpha0v (function of time, dependent on x).
%
[alpha0v]=Getalpha0psiTbar(M0,x,t1v,Phibarstp);
%
maxrelbsalpproperr=0;
%
% Call propalpnlpsiTbar.
%
[allalpham,t1vfinal,Phibarst,relbacksuberr]=...
    propalpnlpsiTbar(M0,Abarv,Abarpv,Phibarst,Phibarstp,...
    Dssig,Dpsigs,DMz,alpha0v,...
    x1,c1,t0,Tbar,nt1,numstep,funcparams);
%[allalpham,t1vfinal,Phibarst,relbsalpproperr]=...
%    propalpnlpsiTbar(M0,Abarv,Psiv,...
%    x,c1,t0,Tbar,nt1,numstep,funcparams);
%
    numalpfuncs=size(allalpham,1);
    alphafinalv=allalpham(numalpfuncs,:);
    maxrelbsalpproperr=max(maxrelbsalpproperr,relbsalpproperr);
    ['The post-propagation relative alpha backsubs error = ',...
        num2str(maxrelbsalpproperr)]
    %
    %
    % Check whether this is really stat or not, if desired.
    %
    checkingstat=1;
    if checkingstat==1
      %alphav4chk=allalpham(numalpfuncs,:)
      [effectofhalving,shouldbeahalf]=checkstat(alphafinalv,x,t1v,Phibarst,...
      M0,c1,Gammatilde,funcparams);
      'Checking if alpha really stat'
      ['The effect of halving (should be near 0.25) is ',num2str(effectofhalving)]
      ['Check on halving (should be near 0.5) is ',num2str(shouldbeahalf)]
      'Pausing in tprop...'
      pause;
    end;
    %
    % Get q_s and r_s for the staticizing alpha
    %
    [qs,rs,qv,rv]=qrintegrals(t1vfinal,alphafinalv,Phibarst,...
        M0,c1,Gammatilde,funcparams);
    %
    % Compute W(t1v(is),x)=0.5*x'*Pv(:,:,is)*x+x'*qv(1,is)+0.5*rv(1,is)
    %
    Wv=zeros(1,ntotp1);
    for it=1:ntotp1
      Wv(1,it)=0.5*x'*Pv(:,:,it)*x+x'*qv(:,it)+0.5*rv(1,it);
    end;
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
% Check whether the initial value of alpha is correct.
%
x1=x(1,1);
[alpha1,alphachk1]=Checkinitialalpha(alphafinalv,x1,c1,funcparams)
%it=1;
%x1=x(1,1);
%[Ntilde,alphavfromx,thetafromx]=Ntildeviafxdpt(x1,c1,funcparams);
%      alpchk=alphafinalv(1,1)
%      alpchk2=alphavfromx(1,1)
%      [thetachk]=gettheta(alpchk,funcparams)
%      thetachk2=thetafromx
%      [gradthetachk]=gradtheta(alpchk,funcparams);
%      [gradthetachk2]=gradtheta(alpchk2,funcparams);
%      chkN=invc1*gradthetachk-(alpchk-x1)
%      chkN2=invc1*gradthetachk2-(alpchk2-x1)


