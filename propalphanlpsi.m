function [allalpham,t1vfinal,Phibarst,relbacksuberr]=...
    propalphanlpsi(M0,Abarv,Psiv,x,...
    c1,t0,t10,nt0,numstep,funcparams)
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
%   Abarv - array of time-dependent linear dynamics associated to q
%   Psiv - Fundamental matrix function (vector) associated Abarv.
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
% functions called:
%   GetalphaPhifuncnlpsi[
%                        getDs,Getalphafuncnl[
%                                             gradtheta,newtonforeta
%                                             ]
%                        ],
%   gradeta,
%   backsubschkalphanl[gradtheta],
%   getderofalpwrtTnlpsi,
%   GetdalpTdsatTnlpsi,
%   GetPhialpha0psi,
%   getDs
%
%
%
% Retrieve data.
%
M0p=M0';
dim=size(x,1);
%
% Set up initial time arrays and step size.
%
nt=nt0;
ntp1=nt+1;
nfull=nt0+numstep;
nfullp1=nfull+1;
t1=t10; % terminal time.
t1v=linspace(t0,t1,ntp1); 
t0v=t1v;
del=t1v(1,2)-t1v(1,1);
invdel=1/del;
%
% Get transposes of Abarv matrices.
%
Abarpv=zeros(dim,dim,nfullp1);
for it=1:nfullp1
    Abarpv(:,:,it)=Abarv(:,:,it)';
end;
%
% Set up the array that will hold the output data.
%
allalpham=zeros(numstep+1,nfullp1);
%
% Get the alpha function at our initial time. Specifically,
% get Phi, Phi' matrices, the initial alpha values and alpha function
% (solution of (8.5) with L^0,\Theta=0).
% This wil be done by a fixed-point method.
%
% First select which fixed-point method to use.
% The second one is the one that's actually guaranteed to work
% for sufficiently short time interval.
%
whichmethod=2;
%
% Now, actually do the fixed-point method for obtaining the alpha function
% at our initial time.
%
first=1;
alphastartv=zeros(1,ntp1);
[alphav,Phibarst,Phibarstp,Dssig,Dpsigs,alpha0v]=...
    GetalphaPhifuncnlpsi(M0,c1,x,t1v,Psiv,first,alphastartv,...
    funcparams,whichmethod);
alphaorigv=alphav;
%
% Check these state-transition matrices by derivatives if desired.
%
chkstmders=0;
if chkstmders==1
  for it0=1:nfullp1
    for it1=1:nfullp1
      'State-transition matrix check in GetalphaPhifuncnlpsi.'
      [it1,it0]
      der1Phibarst(:,:,it1,it0)=invdel*(Phibarst(:,:,it1+1,it0)-...
          Phibarst(:,:,it1,it0));
      der2Phibarst(:,:,it1,it0)=invdel*(Phibarst(:,:,it1,it0+1)-...
          Phibarst(:,:,it1,it0));
      der1Phibarst(:,:,it1,it0)
      Abarv(:,:,it1)*Phibarst(:,:,it1,it0)
      der2Phibarst(:,:,it1,it0)
      -Phibarst(:,:,it1,it0)*Abarv(:,:,it0)
      'Pausing.'
      pause;
    end
  end
end
%
% Check solution.
% Only works in linear case!!!!!!
%
% First, get gradient of \eta wrt \alpha.
%
[gradetav]=gradeta(alphav,c1,funcparams);
%
% Get diagonal matrix from gradient of \eta wrt \alpha.
%
matgradeta=diag(gradetav);
%
% Do a back-subsitution check of the output of GetalphaPhifuncnlpsi.
%
alphachkv=alphav;
[backsuberr,relbacksuberr]=backsubschkalphanl(alphachkv,...
    M0,c1,x,t1v,Dssig,Dpsigs,alpha0v,funcparams);
['Relative back substitution error at first step= ',num2str(relbacksuberr)]
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
  checkingthistime=0;
  alphadirectchk=1;
  if (alphadirectchk==1)&&(istep==numstep)
  %if (alphadirectchk==1)&&((istep==numstep)||(istep==1))
  %if (alphadirectchk==1)
    checkingthistime=1;
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
    if istep==1
      first=1;
      alphastartv=zeros(1,ntp2);
    else
      first=0;
      alphastartv=[alphav(1,1),alphav];
    end;
    [alphavp,Phibarstpnotp,Phibarstpp,Dssigp,Dpsigsp,alpha0vp]=...
      GetalphaPhifuncnlpsi(M0,c1,x,t1vp,Psiv,first,alphastartv,...
      funcparams,whichmethod);
    %
    % Do a back-substitution check of the output of GetalphaPhifunc.
    %
    alphachkv=alphavp;
    %
    % Might need to update backsubschkalphanl!
    %
    %[backsuberrp,relbacksuberrp]=backsubschkalphanl(alphachkv,...
    %  M0,c1,x,t1vp,Dssigp,Dpsigsp,alpha0vp,funcparams);
    %['(possibly) Post-prop relative back-subs error IF used fxd-pt = ',...
    %    num2str(relbacksuberrp)]
    %['(possibly) Post-prop relative back-subs error IF using fxd-pt = ',...
    %    num2str(relbacksuberrp),'. Pausing.']
    %pause;
    %
    % Get numerical derivative approximations.
    %
    havenumericalderapproxs=1;
    % %alpha0numdershouldbezero=invdel*(alpha0vp(1,2:ntp2)-alpha0v);
    alpha0numder=invdel*(alpha0vp(1,1:ntp1)-alpha0v);
    alphanumder=invdel*(alphavp(1,1:ntp1)-alphav);
    %(dalpdtv-alphanumder)./(abs(dalpdtv)+abs(alphanumder))
  end; % end of if for checking alpha.
  %
  % Check of our actual propagating solution.
  %
  alpha0vforchk=zeros(1,ntp1);
  for it0=1:ntp1
    alpha0vforchk(1,it0)=M0*Phibarstp(:,:,ntp1,it0)*x;
  end;
  %
  % Might need to update backsubschkalphanl!
  %
  %[backsuberrp,relbacksuberrp]=backsubschkalphanl(alphav,...
  %    M0,c1,x,t1v,Dssig,Dpsigs,alpha0vforchk,funcparams);
  %['(possibly) Post-prop relative back-subs error = ',...
  %      num2str(relbacksuberrp),' at istep= ',num2str(istep)]
  %'Pausing.'
  %pause;
  %
  %
  % Get gradient of \eta wrt \alpha.
  %
  [gradetav]=gradeta(alphav,c1,funcparams);
  %
  % Get diagonal matrix from gradient of \eta wrt \alpha.
  %
  matgradeta=diag(gradetav);
  %
  % Get the derivative of alpha wrt T (dalpdtv).
  % (Also derivative of alpah_0 wrt T (dadTxpartv).)
  %
  [dalpdtv,dadTxpartv]=getderofalpwrtTnlpsi(M0,c1,x,t1v,Abarpv,...
      alphav,Phibarst,Phibarstp,Dssig,Dpsigs,matgradeta);
  %
  % check alpha derivative at first step.
  %
  %if istep==1
  %    dalpdtv
  %    alphanumder
  %    'pausing after der chk.'
  %    pause;
  %end;
  %
  % Check solution of ODE.
  % Was fab. Currently not checking to save time.
  %[swrelnormsdererr,relnormdererr]=chkderwrtT(alphav,dalpdtv,...
  %  Phibarstp,Dssig,Dpsigs,t1v,x,c1,Abarv,M0,funcparams);
  %['relnormdererr = ',num2str(relnormdererr)]
  %' Pausing in propalphanlpsi.'
  %pause;
  %
  %
  % Get $\frac{d\alpha^T_s}{ds}\vert_{s=T}.
  %
  [dalpTdsatT]=GetdalpTdsatTnlpsi(c1,x,M0,Abarpv,t1v,alphav,Dpsigs,matgradeta);
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
  %if (alphadirectchk==1)&&(istep==numstep)
  if checkingthistime==1
    %
    % Compare with the one generated by direct computation.
    % Also include alpha at the original time for reference.
    %
    %alphaorigv
    %alphanextv
    %alphavp
    reldiffalp=norm(alphanextv-alphavp,2)/...
      (norm(alphanextv,2)+norm(alphavp,2));
    ['Comparison w/ fxd-pt comp (questionable if numstep high) is ',...
        num2str(reldiffalp)]
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
  plotalpha=0;
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
  %istep
  %
  % We also need updated Phi, Phi', Dssig and Dpsigs matrices and
  % initial alpha (alpha_0) values.
  %
  [alpha0v]=GetPhialpha0psi(M0,c1,x,t1v,Phibarst,Phibarstp);
  %
  [Dssig,Dpsigs]=getDs(x,t1v,Phibarst,Phibarstp);
  %
  % A derivative check.
  %
  %'Comparison of num and anal der wrt s at s=t'
  %(alphav(1,ntp1)-alphav(1,ntp1-1))/del
  %dalpTdsatT
  %pause;
  %
  %['At end of a loop step:',num2str(istep)]
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
%'Terminal back substitution error and relative version.'
%backsuberr
['Rel backsub err from backsubchkalphanl in propalphanlpsi= ',num2str(relbacksuberr)]
%relbacksuberr
%pause;
