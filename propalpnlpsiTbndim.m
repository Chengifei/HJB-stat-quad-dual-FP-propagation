function [allalpham,t1vfinal,relbacksuberr]=...
    propalpnlpsiTbndim(M0,Abarv,Phibarst,...
    DMz,alpha0v,c1,tfullv,nt0,funcparams)
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
%   Phibarst,Phibarstp,Dssig,Dpsigs,DMz - see paper.
%   alpha0v - initial guess for alpha
%   x1 - first component of a specific initial state
%        This is only needed as an input to Checkinitialalpha (a
%          function that is not needed for propagation, that is, it's only
%          a check that is used in debugging).
%   c1 - coefficient used in the stat-quad duality
%   tfullv - time discretization points
%   nt0 - number of steps over [0, tfullv(nt0)] for getting alpha
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
%   GetalpPhinlpsiTbar[
%                      getDs,Getalphafuncnl[
%                                           gradtheta,newtonforeta
%                                          ]
%                      extendfixedptalp[etainverse]
%                     ],
%   gradeta,
%   backsubsalpnlTbar[gradtheta],
%   gdalpTdTnlpsiTbar,
%   GetPhialp0psiTbar,
%   getDs
%
%
%
% Retrieve data.
%
Abarpv = permute(Abarv, [2 1 3]);
%dim=size(x,1);
dim=size(Abarpv,1);
%
% Set up initial time arrays.
%
nfullp1 = length(tfullv);
nfull = nfullp1 - 1; % full number of steps (-1 without 0).
ntp1=nt0+1;
numstep=nfull-nt0;
t0 = tfullv(1);
%
% Step size.
%
del=tfullv(2)-tfullv(1);
invdel=1/del;
%
% Get T (see note/paper).
%
T=tfullv(1,ntp1);
%['ntp1 (nt1+1) and T are ',num2str(ntp1),' and ',num2str(T)]
%pause;
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
% Now, actually do the fixed-point method for obtaining the alpha function
% at our initial time.
%
alphastartv=zeros(1,ntp1);
[alphav]=...
    GetalpPhinlpsiTbar(M0,c1,tfullv,Phibarst,DMz,...
    nt0,alpha0v,alphastartv,funcparams);
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
      pause
    end
  end
end
%
% Get gradient of \eta wrt \alpha.
%
[gradetav]=gradeta(alphav,c1,funcparams);
%
% Get diagonal matrix from gradient of \eta wrt \alpha.
%
matgradeta=diag(gradetav);
%
% Do a back-subsitution check of the output of GetalpPhinlpsiTbar if
% desired.
%
chkoutputat1ststep=0;
if chkoutputat1ststep==1
  [backsuberr,relbacksuberr]=backsubsalpnlTbar(alphav,...
       M0,c1,dim,tfullv,DMz,alpha0v,ntp1,funcparams);
  ['In propal...iTbar, rel back-sub err (from back...lTbar) at 1st step= ',num2str(relbacksuberr)]
  %'Pausing after first backsuberror chk in tprop...'
  %pause;
end;
%
% Store output
%
allalpham(1,1:nfullp1)=alphav;
%
% Begin propagation.
%
%'Start ODE propagation with indep var T.'
for istep=1:numstep
  %
  % Get T and its location.
  %
  locTp1=nt0+istep;
  T=t0+del*locTp1;
  %['in propalpnl..., locTp1 and T are ',num2str(locTp1),' and ',num2str(T)]
  %pause;
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
  [dalpdtv,dadtalpTpartv,matgradeta,derofintmat]=...
      getdalpTdTnlpsiTbar(M0,c1,tfullv,Abarpv,...
      alphav,DMz,matgradeta,locTp1);
  %
  % Include a one-step check if desired.
  % This assumes we can use the fixed-point method to get
  % alpha at these later times.
  % Appears in two places at each step.
  %
  checkingthistime=0;
  alphadirectchk=0;
  %if (alphadirectchk==1)&&(istep==numstep)
  if (alphadirectchk==1)&&((istep==numstep)||(istep==1))
    checkingthistime=1;
    %
    % Modify terminal time to look at rate of change of alpha function
    % wrt terminal time.
    %
    locT=locTp1-1;
    %
    % Get initial alpha values and alpha function
    % (solution of (8.5) with L^0,\Theta=0).
    % Get this at the NEXT time step!
    %
    alphastartv=zeros(1,locTp1+1);
    [alphavpchk]=...
      GetalpPhinlpsiTbar(M0,c1,tfullv,Phibarst,...
      DMz,locTp1,alpha0v,alphastartv,funcparams);
    %
    % Checks on derivative calculations.
    %
    % Get numerical derivative approximations.
    %
    alphanumder=invdel*(alphavpchk-alphav);
    'Num derivative by fixed-pt code. Will be poor if T too large.'
    alphanumder
    'Diff between numerical and analytical:'
    (dalpdtv-alphanumder)./(abs(dalpdtv)+abs(alphanumder))
    ['istep= ',num2str(istep)]
    %
    % Get the three terms individually at the next time step if desired.
    %
    wantmorechecks=0;
    if wantmorechecks==1
      invc1=1/c1;
      %
      [gradthetavpchk]=gradtheta(alphavpchk,funcparams);
      etavnextpchk=alphavpchk-invc1*gradthetavpchk;
      [gradthetavnow]=gradtheta(alphav,funcparams);
      etavnow=alphav-invc1*gradthetavnow;
      secnumderterm=invdel*(etavnextpchk-etavnow);
      'Eta/2nd terms (previously first terms) compared:'
      analderetachk=(matgradeta*(dalpdtv'))'
      numderetachk=(matgradeta*(alphanumder'))'
      secnumderterm
      norm(analderetachk-secnumderterm,2)/(norm(analderetachk,2)+...
          norm(secnumderterm,2))
      norm(numderetachk-secnumderterm,2)/(norm(numderetachk,2)+...
          norm(secnumderterm,2))
      'end of eta (2nd) terms compared.'
      %
      '3rd terms compared:'
      % Checking individual parts.
      [intonrhsvnow]=intonrhsofetaofalp(alphav,tfullv,locT,c1,M0,DMz);
      [intonrhsvmid]=intonrhsofetaofalp(alphavpchk,tfullv,locT,c1,M0,DMz);
      [intonrhsvpchk]=intonrhsofetaofalp(alphavpchk,tfullv,locTp1,c1,M0,DMz);
      numder3rdterm=invdel*(intonrhsvpchk-intonrhsvnow)
      %
      numder3part2=invdel*(intonrhsvpchk-intonrhsvmid)
      dadtalpTpartv
      %
      numder3part1=invdel*(intonrhsvmid-intonrhsvnow)
      (derofintmat*(alphanumder'))'
      %(derofintmat*(dalpdtv'))'
      fullanalder3rdterm=dadtalpTpartv+(derofintmat*(dalpdtv'))'
      norm(numder3rdterm-fullanalder3rdterm,2)/(norm(numder3rdterm,2)+...
          norm(fullanalder3rdterm,2))
      'end of 3rd terms compared:'
      pause;
    end; % end of wantmorechecks
    %
  end; % end of if for checking alpha.
  %'hello6 in prop...'
  %
  % Back subs error by backsubschkalphanl at current time step if desired.
  %
  wantbacksuberreachstep=0;
  if wantbacksuberreachstep==1
    [backsuberr,relbacksuberr]=backsubsalpnlTbar(alphav,...
      M0,c1,dim,tfullv,DMz,alpha0v,locTp1,funcparams);
    'Relative Back-substitution error.'
    relbacksuberr
    'alphav:'
    alphav
    istep
    pause;
  end;
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
  %
  % Propagate \alpha^T_\cdot forward as alphanextv.
  %
  alphanextv=alphav+del*dalpdtv;
  alphapreviousv=alphav;
  %
  % Continue the one-step check if it was desired.
  %
  %if (alphadirectchk==1)&&(istep==numstep)
  if checkingthistime==1
    %
    % Compare with the one generated by direct computation.
    % Also include alpha at the original time for reference.
    %
    reldiffalp=norm(alphanextv-alphavpchk,2)/...
      (norm(alphanextv,2)+norm(alphavpchk,2));
    ['Comparison w/ fxd-pt comp (questionable if numstep high) is ',...
        num2str(reldiffalp)]
    %'Pausing in tprop...'
    %pause;
  end;
  %
  % Plot the alpha function if desired.
  %
  plotalpha=0;
  if plotalpha==1
    if istep==1
        figure;
        plot(tfullv,alphav);
        %axis([0 1 -3 3]);
        hold on;
    elseif istep==numstep
        plot(tfullv,alphav);
        hold off;
    else
        plot(t1vold,alphav);
    end;
  end;
  %
  % Now update the alpha function/vector.
  %
  alphav=alphanextv;
  %
  % Store output
  %
  allalpham(istep+1,1:nfullp1)=alphav;
  if istep==numstep
    t1vfinal=tfullv;
  end;
  %
  % Also update locTp.
  %
  locTp2=locTp1+1;
  %
  % A derivative check.
  %
  %'Comparison of num and anal der wrt s at s=t'
  %(alphav(1,ntp1)-alphav(1,ntp1-1))/del
  %dalpTdsatT
  %pause;
  %['At end of a loop step:',num2str(istep)]
  %
end;
%
% Do a back-subsitution check of the output of GetalphaPhifunc.
%
[backsuberr,relbacksuberr]=backsubsalpnlTbar(alphav,...
    M0,c1,dim,tfullv,DMz,alpha0v,locTp1,funcparams);
%['Rel norm of backsub err from backsubchkalphanl in propalphanlpsi= ',num2str(relbacksuberr)]
%pause;
