function [alphav]=Getalphafuncnl(M0,c1,x,t1v,...
    DMz,alpha0v,alphastartv,funcparams,whichmethod)
%
% N.B.: HAS A HARDWIRED Theta FUNCTION (dependent on funcparams only)!
%
% Uses the fixed=point method for computation of alpha vector
% representing the alpha function of time.
% This is the nonlinear version, i.e., for $\grad\Theta\not=0$.
%
% input variables;
%   M0,c1,x,t1v - parameters
%   Dssig,Dpsigs - internal variables to the code
%   DMz - internal variable to the code
%   alpha0v - alpha function/vector at time zero (the function of x).
%   alphastartv - starting point for the fixed-point iteration.
%   funcparams - the parameters ([ktemp,eps]).
%   whichmethod - 1 for fixed-point method 1 and 2 for method 2.
%
% output variables:
%   alphav - 
%
% functions called:
%   gradtheta,newtonforeta
%
%
% Recover data.
%
M0p=M0';
invc1=1/c1;
t0v=t1v;
ntp1=size(t1v,2);
nt=ntp1-1;
t0=t1v(1,1);
t1=t1v(1,ntp1);
del=t1v(1,2)-t1v(1,1);
delo2=del/2;
%
% convergence criterion for fixed-point method.
%
convcrit=1e-6;
nfpmax=800;
normoversumov=[];
%
% First, get the over-all-time-part of the derivative.
% This is a "time-by-time-dimensioned" matrix.
%
%dim=size(x,1);
%DMz=zeros(dim,ntp1,ntp1);
%for ik=1:ntp1
% for ij=1:ntp1
%  if ij<ik
%    DMz(:,ik,ij)=squeeze(Dpsigs(:,:,ij,ik)*M0p);
%  else
%    DMz(:,ik,ij)=squeeze(Dssig(:,:,ik,ij)*M0p);
%  end
% end
%end
%
% If desired, check against analytical computation valid only
% for this Abar.
%
wantDscheck=0;
if wantDscheck~=0 
  Dssigchk=zeros(dim,dim,ntp1,ntp1);
  Dpsigschk=zeros(dim,dim,ntp1,ntp1);
  DMzchk=zeros(dim,ntp1,ntp1);
  cosinesmat=zeros(ntp1,ntp1);
  sinesmat=zeros(ntp1,ntp1);
  for isig=1:ntp1
    for is=1:ntp1
      cosinesmat(isig,is)=cos(t1v(1,isig)-t1v(1,is));
      sinesmat(isig,is)=sin(t1v(1,isig)-t1v(1,is));
    end
  end
  for isig=1:ntp1
    for is=1:ntp1
      Dpsigschk(1,1,isig,is)=t1v(1,isig)*cos(t1v(1,isig)-t1v(1,is));
      Dpsigschk(1,2,isig,is)=-t1v(1,isig)*sin(t1v(1,isig)-t1v(1,is));
      %Dpsigschk(1,1,isig,is)=sin(t1v(1,isig)-t1v(1,is));
      %Dpsigschk(1,2,isig,is)=cos(t1v(1,isig)-t1v(1,is));
      Dpsigschk(2,1,isig,is)=-Dpsigschk(1,2,isig,is);
      Dpsigschk(2,2,isig,is)=Dpsigschk(1,1,isig,is);
      Dssigchk(1,1,is,isig)=t1v(1,is)*cos(t1v(1,isig)-t1v(1,is));
      Dssigchk(1,2,is,isig)=-t1v(1,is)*sin(t1v(1,isig)-t1v(1,is));
      Dssigchk(2,1,is,isig)=-Dssigchk(1,2,is,isig);
      Dssigchk(2,2,is,isig)=Dssigchk(1,1,is,isig);
    end
  end
  for isig=1:ntp1
    for is=1:ntp1
      if isig<is
        DMzchk(:,is,isig)=Dpsigschk(:,:,isig,is)*M0p;
      else
        DMzchk(:,is,isig)=Dssigchk(:,:,is,isig)*M0p;
      end
    end
  end
end; % end of want DS check if
%'Dssig check in Getalphafunc...'
%sum(sum(sum(sum(abs(Dssigchk(:,:,2:ntp1,2:ntp1)-Dssig(:,:,2:ntp1,2:ntp1))))))
%'Dpsigs check in Getalphafunc...'
%sum(sum(sum(sum(abs(Dpsigschk-Dpsigs)))))
%'DMz check in Gcetalphafunc...'
%sum(sum(sum(abs(DMzchk-DMz))))
%pause;
%
% FIXED POINT LOOP TO SOLVE FOR ALPHA FUNCTION.
%
if whichmethod==1
%
% METHOD 1:
%   - Quicker computationally.
%   - Uses alpha=(1/c1)\grad\theta(\alpha)+[the other terms].
%
% Note(!!!): Method 1 is not currently being used, and the code within may
% be defunct.
%
%
for ifp=1:nfpmax
  if ifp==1
    alpha00v=alphastartv;
  else
    alpha00v=alphafixedptv;
  end;
  %
  % Get update of alphav in our fixed-point method.
  %
  % First, we need to do the integral that's on the RHS.
  %
  mc1delM0=-c1*del*M0;
  intforalp=zeros(1,ntp1);
  for ialp=2:ntp1
    if ialp<3
      for isig=1:nt
        intforalp(1,ialp)=intforalp(1,ialp)+mc1delM0*DMz(:,ialp,isig)*...
           alpha00v(1,isig);
      end;
    else
      intforalp(1,ialp)=0.5*mc1delM0*(DMz(:,ialp,1)*alpha00v(1,1)+...
        DMz(:,ialp,ntp1)*alpha00v(1,ntp1));
      for isig=2:nt
        intforalp(1,ialp)=intforalp(1,ialp)+mc1delM0*DMz(:,ialp,isig)*...
            alpha00v(1,isig);
      end;
     end;  % end of if ialp<3
  end; % end of for ialp=...
  %
  % Get $\grad\Theta(\alpha)$ for alpha in alpha00v.
  %
  [gradthetav]=gradtheta(alpha00v,funcparams);
  %
  % New iterate of alpha in the fixed-point method.
  %
  alphav=invc1*gradthetav+alpha0v+intforalp;
  %
  % Check convergence of the fixed point method.
  %
  if ifp>1
    normdiff=norm(alphafixedptv-alphav,2);
    sumonorms=norm(alphafixedptv,2)+norm(alphav,2);
    normoversumo=normdiff/sumonorms;
    normoversumov=[normoversumov,normoversumo];
    %pause;
  end;
  %
  alphafixedptv=alphav;
  %
  successoffp=0;
  if (ifp>1)&&(normoversumo<convcrit)
    successoffp=1;
    exitat=ifp;
    %
    % Exit F-P loop because we have convergence.
    %
    break;
    %
  end; % end of if ifp>1
end; % end of for ifp=
else % if for whichmethod.
%
% ********************************************************
%                                                        *
% Method 2:                                              *
%    - slower computationally                            *
%    - Uses Newton's method inside the fixed point loop  *
%      to solve \eta(\theta)=                            *
%                                                        *
% ********************************************************
%
for ifp=1:nfpmax
  if ifp==1
    alpha00v=alphastartv;
  else
    alpha00v=alphafixedptv;
  end;
  %
  % Get update of alphav in our fixed-point method.
  %
  % First, we need to do the integral that's on the RHS.
  %
  mc1delM0=-c1*del*M0;
  intforalp=zeros(1,ntp1);
  for ialp=2:ntp1
    if ialp<3
      for isig=1:nt
        intforalp(1,ialp)=intforalp(1,ialp)+mc1delM0*DMz(:,isig,ialp)*...
           alpha00v(1,isig);
      end;
    else
      intforalp(1,ialp)=0.5*mc1delM0*(DMz(:,1,ialp)*alpha00v(1,1)+...
        DMz(:,ntp1,ialp)*alpha00v(1,ntp1));
      for isig=2:nt
        intforalp(1,ialp)=intforalp(1,ialp)+mc1delM0*DMz(:,isig,ialp)*...
            alpha00v(1,isig);
      end;
     end;  % end of if ialp<3
  end; % end of for ialp=...
  %
  % Invoking Newton's method for the inside loop to solve \eta(alpha)=RHS.
  %
  rhsstuffv=alpha0v+intforalp;
  deltaf=1e-8;
  maxit=25;
  [alphasolv,falphasolv,successflagnewt]=newtonforeta(rhsstuffv,c1,...
    funcparams,deltaf,maxit);
  if successflagnewt==0
      'Newton method for the inside loop of fixed-point method failed'
      ' inside Getalphafuncnl. Pausing.'
      pause;
  end;
  %
  % If desired, check Newton sol relative to the other approach calcs.
  %
  checksol=0;
  if checksol==1
    [gradthetav]=gradtheta(alphasolv,funcparams);
     chkofnmsol=alphasolv-(invc1*gradthetav+alpha0v+intforalp);
     sumofmags=norm(alphasolv,2)+norm(invc1*gradthetav,2)+norm(alpha0v,2)+...
         norm(intforalp,2);
     'The check of the Newton sol rel to the other method equation'
     'Relative norm of error is'
     relchkofnmsol=norm(chkofnmsol,2)/sumofmags
     %pause;
  end;
  %
  % New iterate of alpha in the fixed-point method.
  %
  alphav=alphasolv;
  %
  % Check convergence of the fixed point method.
  %
  if ifp>1
    normdiff=norm(alphafixedptv-alphav,2);
    sumonorms=norm(alphafixedptv,2)+norm(alphav,2);
    normoversumo=normdiff/sumonorms;
    normoversumov=[normoversumov,normoversumo];
    %['In Getalphafuncnl, normoversumov= ',num2str(normoversumov)]
  end;
  %
  alphafixedptv=alphav;
  %
  successoffp=0;
  if (ifp>1)&&(normoversumo<convcrit)
    successoffp=1;
    exitat=ifp;
    %
    % Exit F-P loop because we have convergence.
    %
    break;
    %
  end; % end of if ifp>1
end; % end of for ifp=
end % if for whichmethod.
%
if successoffp==1
    ['In Getalphafuncnl exited successfully at it= ',num2str(exitat)...
       ' with rel err: ',num2str(normoversumo)]
    %exitat
    %normoversumo
    %pause;
else
    'Went to the max number of iterations in F-P loop. Relative diff:'
    'This implies that those checks that rquire direct f-p computation'
    'of alpha no longer work.'
    'One way this can happen is if the propagation has gone on so long'
    'that the fixed-point based check can no longer be used.'
    'The back-substitution check will still be valid.'
    'The sequence of relative errors, normoversumov is:'
    normoversumov
    'Pausing.'
    pause;
end;