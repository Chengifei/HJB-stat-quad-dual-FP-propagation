function [allalpham,t1vfinal]=propalphanl(M0,Abar,x,c1,t0,t10,...
    nt0,numstep,funcparams)
%
% function form of tpropagatealphanl
%
% First, for the "initial terminal time", it uses a fixed-point method
% to obtain alpha^T_\cdot where superscript T indicate the terminal time.
% Note that for T (initial terminal time) sufficiently small, the 
% fixed point method will work.
% Then, it
% propagates \alpha^T_\cdot as a function of terminal time T forward from
% the initial terminal value (t10) forward numstep steps.
%
% Currently assuming L^0=0, m=1 (and maybe n=2).
%
%
% input variables:
%   M0 - matrix that defines the state-components in the nonlinearity
%   Abar - nominal linear dynamics
%   x - a specific initial state
%   c1 - coefficient used in the stat-quad duality
%   t0 - initial time (constant)
%   t10 - initial terminal time
%   nt0 - number of steps over [t0,t10] for getting alpha with that
%           initial terminal time.
%   numstep - the number of steps the code will propagate t10 forward
%               to a later terminal time. Note that the step-size
%               is h=(t10-t0)/nt0, and hence it propagates forward to
%               terminal terminal time T_f=t10+numstep*h.
%   funcparams - the parameters ([ktemp,eps]).
%
% output variables:
%    allalpham - an (numstep+1)x(nto+numstep+1) array where the ith row
%                 will be alpha^T_\cdot. This row will be zeros after
%                   approximately the nt0+i+1 entry.
%   t1vfinal - 1x(nto+numstep+1) vector of timesteps corresponding to
%               the longest alpha, alpha^(t10+h*numstep)_\cdot
%
%
% Retrieve data.
%
%M0=[1 0]; % matrix defining lower dim'l set on which nonlinearities live
M0p=M0';
%Abar=[0 1; -1 0]; % defines the linear part of the nominal dynamics
Abarp=Abar';
%x=[3;4]; % initial state
%x=10*[3;4]; % initial state
dim=size(x,1);
%c1=2; % coefficient related to stat-quad duality
%c1=-1.5; % coefficient related to stat-quad duality
%c1=4; % coefficient related to stat-quad duality
%c1=2.5; % coefficient related to stat-quad duality
%c1=2.75; % coefficient related to stat-quad duality
%
% We'll be looking at t0\le t \le s \le T \le t1
% where t is initial and T is terminal.
%
%t0=0;       % initial time (generally zero)
%t10=0.5;    % base terminal time (terminal time will vary)
%
% Set number of time steps for approximation of the alpha function
% over the initial duration from t0 to t10.
%
%nt0=16;
%nt0=32;
%nt0=64;
%
% Number of steps in propagation of the terminal time.
% This is the number of steps forward from (terminal time) t10 that
% the code will propagate, thereby extending the terminal time from t10
% to a later time.
%
%numstep=16;
%numstep=32;
%numstep=64;
%
% Set up initial time arrays and step size.
%
nt=nt0;
ntp1=nt+1;
t1=t10; % terminal time.
t1v=linspace(t0,t1,ntp1); 
t0v=t1v;
del=t1v(1,2)-t1v(1,1);
invdel=1/del;
%
% Set up the array that will hold the output data.
%
allalpham=zeros(numstep+1,nt0+numstep+1);
%
% Get the alpha function at our initial time. Specifically,
% get Phi, Phi' matrices, the initial alpha values and alpha function
% (solution of (8.5) with L^0,\Theta=0).
%
% First select which fixed-point method to use.
%
whichmethod=2;
%
first=1;
alphastartv=zeros(1,ntp1);
[alphav,Phibarst,Phibarstp,Dssig,Dpsigs,alpha0v]=...
    GetalphaPhifuncnl(M0,c1,x,t1v,...
    first,alphastartv,funcparams,whichmethod);
alphaorigv=alphav;
%
% Check solution.
% Only works in linear case!!!!!!
%
% First, get gradient of \eta wrt \alpha.
%
%'gradetav in tprop...'
%[gradetav]=gradeta(alphav,c1);
[gradetav]=gradeta(alphav,c1,funcparams);
%
% Get diagonal matrix from gradient of \eta wrt \alpha.
%
matgradeta=diag(gradetav);
[alphacheckv]=chkalpha_lincase_only(M0,c1,x,t1v,Abar,...
    alphav,alpha0v,Dssig,Dpsigs,matgradeta);
%
% Do a back-subsitution check of the output of GetalphaPhifunc.
%
alphachkv=alphav;
[backsuberr,relbacksuberr]=backsubschkalphanl(alphachkv,...
    M0,c1,x,t1v,Dssig,Dpsigs,alpha0v,funcparams);
'Back substitution error and relative version at first step.'
backsuberr
relbacksuberr
%'Pausing after first backsuberror chk in tprop...'
%pause;
%
% Store output
%
allalpham(1,1:ntp1)=alphav;
%
% Begin propagation.
% Because of the odd aspect that the terminal time increases at
% each step, we are currently only trying an Euler's method.
%
for istep=1:numstep
  %
  % Include a one-step check if desired.
  % This assumes we can use the fixed-point method to get
  % alpha at these later times.
  % Appears in two places at each step.
  %
  havenumericalderapproxs=0;
  alphadirectchk=1;
  if (alphadirectchk==1)&&(istep==numstep)
  %if (alphadirectchk==1)
    %'Not skipping check. Pausing.'
    %pause;
    %if alphadirectchk==1
    %
    % Modify terminal time to look at rate of change of alpha function
    % wrt terminal time.
    %
    t1p=t1v(1,ntp1)+del; % new terminal time.
    ntp2=ntp1+1;
    t1vp=linspace(t0,t1p,ntp2); 
    t0vp=t1vp;
    %
    % Get Phi, Phi' matrices, the initial alpha values and alpha function
    % (solution of (8.5) with L^0,\Theta=0).
    %
    'Doing a check of (possibly propagated) alpha via a direct'
    'fixed-point computation.'
    if istep==1
      first=1;
      alphastartv=zeros(1,ntp2);
    else
      first=0;
      alphastartv=[alphav(1,1),alphav];
    end;
    [alphavp,Phibarstpnotp,Phibarstpp,Dssigp,Dpsigsp,alpha0vp]=...
      GetalphaPhifuncnl(M0,c1,x,t1vp,first,alphastartv,funcparams,whichmethod);
    %
    % Do a back-substitution check of the output of GetalphaPhifunc.
    %
    alphachkv=alphavp;
    [backsuberrp,relbacksuberrp]=backsubschkalphanl(alphachkv,...
      M0,c1,x,t1vp,Dssigp,Dpsigsp,alpha0vp,funcparams);
    'Back substitution error and relative version.'
    backsuberrp
    relbacksuberrp
    %'Pausing after backsuberror chk around lline 124 in tprop...'
    %pause;
    %
    % Get numerical derivative approximations.
    %
    havenumericalderapproxs=1;
    % %alpha0numdershouldbezero=invdel*(alpha0vp(1,2:ntp2)-alpha0v);
    alpha0numder=invdel*(alpha0vp(1,1:ntp1)-alpha0v);
    alphanumder=invdel*(alphavp(1,1:ntp1)-alphav);
    %(dalpdtv-alphanumder)./(abs(dalpdtv)+abs(alphanumder))
  end;
  %
  % Get gradient of \eta wrt \alpha.
  %
  %'gradetav in tprop...'
  %[gradetav]=gradeta(alphav,c1);
  [gradetav]=gradeta(alphav,c1,funcparams);
  %
  % Get diagonal matrix from gradient of \eta wrt \alpha.
  %
  matgradeta=diag(gradetav);
  %
  % Get the derivative of alpha wrt T (dalpdtv).
  % (Also derivative of alpah_0 wrt T (dadTxpartv).)
  %
  [dalpdtv,dadTxpartv]=getderofalphawrtTnl(M0,c1,x,t1v,Abar,...
      alphav,Phibarst,Phibarstp,Dssig,Dpsigs,matgradeta);
  %
  % Get $\frac{d\alpha^T_s}{ds}\vert_{s=T}.
  %
  [dalpTdsatT]=GetdalpTdsatTnl(c1,x,M0,Abar,t1v,alphav,Dpsigs,matgradeta);
  %
  % Propagate \alpha^T_\cdot forward.
  %
  alphanextv=alphav+del*dalpdtv;
  alphapreviousv=alphav;
  %
  % We need alphanextv at the newly added s=T+del point as well.
  %
  alphanextatT=alphanextv(1,ntp1)+del*dalpTdsatT;
  alphanextv=[alphanextv,alphanextatT];
  %
  % Continue the one-step check if it was desired.
  %
  if (alphadirectchk==1)&&(istep==numstep)
  %if (alphadirectchk==1)
    %
    % Compare with the one generated by direct computation.
    % Also include alpha at the original time for reference.
    %
    alphaorigv
    alphanextv
    alphavp
    reldiffalp=norm(alphanextv-alphavp,2)/...
      (norm(alphanextv,2)+norm(alphavp,2))
    %'Pausing in tprop...'
    %pause;
  end;
  %
  % Update the alpha function.
  %
  % First do the time vector.
  %
  t1endold=t1v(1,ntp1);
  t1vold=t1v;
  t1endnew=t1endold+del;
  t1v=[t1v,t1endnew];
  t0v=t1v;
  ntp1=size(t1v,2);
  nt=ntp1-1;
  %
  % Plot the alpha function if desired.
  %
  plotalpha=1;
  if plotalpha==1
    %'Plotting from within tpropagatealpha.'
    if istep==1
        figure;
        plot(t1vold,alphav);
        %axis([0 1 -3 3]);
        hold on;
    elseif istep==numstep
        plot(t1vold,alphav);
        hold off;
    else
        plot(t1vold,alphav);
    end;
  end;
  %
  % Now update the alpha vector.
  %
  alphav=alphanextv;
  %
  % Store output
  %
  allalpham(istep+1,1:ntp1)=alphav;
  if istep==numstep
    t1vfinal=t1v;
  end;
  %
  % We also need updated Phi, Phi', Dssig and Dpsigs matrices and
  % initial alpha (alpha_0) values.
  %
  %'The Phi and alpha0v update could be made more efficient.'
  %
  [Phibarst,Phibarstp,alpha0v]=GetPhialpha0(M0,c1,x,t1v);
  [Dssig,Dpsigs]=getDs(x,t1v,Phibarst,Phibarstp);
  %
  % Check solution if desired.
  % This checks how closely the solution is to the
  % one-step fixed-point solution in the linear case.
  % It is only valid in the linear case because it's obtained
  % from solution of the resulting linear system.
  % Again, this only works in linear case!!!!!!
  %
  linchkdesired=0;
  if (linchkdesired==1)&&(istep==numstep)
    %
    % First, get gradient of \eta wrt \alpha.
    %
    %'gradetav in tprop...'
    %[gradetav]=gradeta(alphav,c1);
    [gradetav]=gradeta(alphav,c1,funcparams);
    %
    % Get diagonal matrix from gradient of \eta wrt \alpha.
    %
    matgradeta=diag(gradetav);
    [alphacheckv]=chkalpha_lincase_only(M0,c1,x,t1v,Abar,...
      alphav,alpha0v,Dssig,Dpsigs,matgradeta);
  end; % linchk...
  %
  % A derivative check.
  %
  %'Comparison of num and anal der wrt s at s=t'
  %(alphav(1,ntp1)-alphav(1,ntp1-1))/del
  %dalpTdsatT
  %
  % Some output checks if available
  %
  % THESE DERIVATIVE CHECKS ARE NOT VALID. SEE TEXT-TO-SCREEN BELOW!
  %
  %if  havenumericalderapproxs==1
  %  'Relative difference in num and anal ders of alpha0v wrt T'
  %  norm(dadTxpartv-alpha0numder,2)/...
  %    (norm(dadTxpartv,2)+norm(alpha0numder,2))
  %  'Relative difference in numerical and analytical ders wrt T'
  %  %(dalpdtv-alphanumder)./(abs(dalpdtv)+abs(alphanumder))
  %  dalpdtv
  %  alphanumder
  %  norm(dalpdtv-alphanumder,2)/(norm(dalpdtv,2)+norm(alphanumder,2))
  %  istep
  %  'That is NOT a valid check after the first step!'
  %  'That is because the numerical derivative is obtained'
  %      'from solutions of the original f-p problem - NOT'
  %      'from the propagated solution!'
  %  pause;
  %end;
  ['At end of a loop step:',num2str(istep)]
  %
  % Do another fixed-point loop step.
  %
end;
%
% Do a back-subsitution check of the output of GetalphaPhifunc.
%
alphachkv=alphav;
[backsuberr,relbacksuberr]=backsubschkalphanl(alphachkv,...
    M0,c1,x,t1v,Dssig,Dpsigs,alpha0v,funcparams);
'Terminal back substitution error and relative version.'
backsuberr
relbacksuberr
