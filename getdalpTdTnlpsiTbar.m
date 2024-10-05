function [dalpdtv,dadtalpTpartv,matgradeta,derofintmat]=...
    getdalpTdTnlpsiTbar(M0,c1,t1v,Abarpv,alphav,DMz,matgradeta,locTp1)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                              %
% !!! NEED TO FIX THIS CODE IN THE CASE OF C_2\NOT=0 !!!!!     %
% !!! NEED TO FIX THIS CODE IN THE CASE OF C_2\NOT=0 !!!!!     %
% !!! NEED TO FIX THIS CODE IN THE CASE OF C_2\NOT=0 !!!!!     %
% !!! NEED TO FIX THIS CODE IN THE CASE OF C_2\NOT=0 !!!!!     %
%                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
% N.B.: ASSUMES $\Theta(\alpha)=\atan(\alpha)$!!!!!!
%
% Note also that multiple items have only been obtained for the particular
% value of T (equiv locTp1) input - not for all T.+
%
% Obtains the derivative of $\alpha^T_\cdot$ wrt $T$.
% This is the nonlinear version, i.e., for $\grad\Theta\not=0$.
%
%
% input variables:
%
%   M0 - Matrix definining variable in nonlinearity
%   c1 - stat-quad coef.
%   dim - space dimension
%   t1v - (horizontal) vector of times
%   Abarpv - array of transposes of time-dependent linear dynamics
%            associated to q
%   alphav - time-indexed (hor) vector of alpha^T(\cdot) values.
%   matgradeta - diagonal matrix whose elements are the gradient of
%               eta wrt alpha.
%   locTp1 - the position of T in t1v, i.e., where T=t1v(1,locTp1).
%          Note that we have \alpha^T_s for t\le s,T \le \Tbar. Typically we
%          take t=0. Also, we're really only interested for s\in[t,T];
%          the use of \Tbar with s\in (T,\Tbar] is a mathematical
%          artifact that is useful for obtaing the function over [t,T].
%
% output variables:
%   dalpdtv - derivative of alpha wrt T.
%
% functions called:
%    none.
%
% First, recover data.
%
%dim=size(x,1);
M0p=M0';
nfullp1=size(t1v,2);
nfull=nfullp1-1;
locT=locTp1-1;
%
del=t1v(1,2)-t1v(1,1);
invdel=1/del;
delo2=del/2;
%
% Get the x part of the derivative.
% dadTxpartv zeroed due to correction of typo in the paper.
%dadTxpartv=zeros(1,nfullp1);
%
% Get the alpha_T part of the derivative.
% Note that this depends on T!
%
mdemic1M0=-0.5*c1*M0;
dadtalpTpartv=zeros(1,nfullp1);
for is=1:nfullp1
  dadtalpTpartv(1,is)=(mdemic1M0*DMz(:,locTp1,is))*alphav(1,locTp1)+...
      (mdemic1M0*DMz(:,locTp1+1,is))*alphav(1,locTp1+1);
end;
%
% Get matrix corresponding to RHS linear operators.
%
mc1delM0=-c1*del*M0;
derofintmat=zeros(nfullp1,nfullp1);
for is=1:nfullp1
  derofintmat(is,1:locTp1)=pagemtimes(mc1delM0, DMz(:,1:locTp1,is));
  derofintmat(is,[1 locTp1]) = derofintmat(is,[1 locTp1]) * 0.5;
end
%
% Obtain the linear system to be solved for the derivative.
% First get the RHS.
%
DeltaV=dadtalpTpartv';
%
% Next, get matrix in linear system we need to solve.
%
matmultdalpdT=matgradeta-derofintmat;
%
% Obtain the derivative from solution of (I-derofintmat)*Deltav,
% by solution of the linear system.
%
dalpdtv=(matmultdalpdT\DeltaV)';

