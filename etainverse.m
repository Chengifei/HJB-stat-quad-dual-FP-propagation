function [alphav]=etainverse(etav,alphastartv,funcparams,c1)
%
% Performs $\eta^{-1}$ to obtain the alpha function from an eta function.
% These are stored as vectors, of course.
% Here, $\eta(\alpha) \doteq \grad\Theta(\alpha) - c_1\alpha $. 
% It performs the inverse via the fixed-point method.
%
%
% input variables:
%   etav - a vector of eta values.
%   alphastartv - starting guess of the alpha values for the fixed-point
%                 iteration.
%   funcparams - the parameters ([ktemp,eps]).
%   c1 - c_1.
%   % internal now. nmax - maximum allowable number of iterations.
%   % internal now. convcrit - convergence criterion for the iteration.
%
% output values:
%   alphav - a vector of alpha=eta^{-1}(eta) values.
%
% functions called:
%   none.
%
%
%
% N.B.: ASSUMES parituclar function form:
%
% output variables:
%   alphav - \eta^{-1}(\eta) for all the values in etav.
%
%
% Assign exit criteria.
%
nmax=100;
convcrit=1e-6;
normdifffailbnd=1e-10;
if convcrit<normdifffailbnd
  'Potential problem in chcecking stop condition in Ntildeviafxdpt.'
  'Pausing.'
  pause;
end;
%
% Starting value for fixed-point iteration.
%
alpha0v=alphastartv;
alphav=alphastartv;
%
% Begin iteration.
%
for ifp=1:nmax
  alpha1v=etav+c1 \ gradtheta(alphav,funcparams);
  %
  % Check convergence.
  %
  if ifp>1
    normdiff=norm(alpha1v-alphav);
    sumonorms=norm(alpha1v,2)+norm(alphav,2);
    if normdiff>normdifffailbnd
      normoversumo=normdiff/sumonorms;
    else
      %'Difference in etainverse so small that not using relative difference!'
      normoversumo=normdifffailbnd;
    end;
    %['in etainverse, normoversumo= ',num2str(normoversumo)]
  end;
  %
  % Update alphav
  %
  alphav=alpha1v;
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
  end;
end; % end of for ifp=1:nmax
if successoffp==1
    %['Exited f-p loop in etainverse with success at step ',num2str(exitat)]
    %exitat
else
    'Inside ETAINVERSE.M!'
    ['Went to the max number of iterations in F-P loop. Relative diff: ',...
        num2str(normoversumo)]
    'Pausing.'
    pause;
end;
%
% Assign result.
%
alphav=alpha1v;