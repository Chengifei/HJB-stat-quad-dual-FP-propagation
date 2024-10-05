function [alpha1,alphachk1]=Checkinitialalpha(alphav,x1,c1,funcparams)
%
% Checks whether the value of alpha at the initial time and x_1
% matches what comes out of the stat-dual at that x_1.
% (This is for a cost nonlinearity that only depends on x_1 - the
% first component of x.)
%
% input variables:
%   alphav - the vector (time-indexed) of alpha values.
%   x1 - the first component of x (the one in which we have a nonlinearity).
%   c1 - c_1
%   funcparams - parameters in the form of the nonlinearity.
%
% output variables:
%   alpha1 - the first element of alphav.
%   alphachk1 - the alpha that it should be, i.e., the value that
%               would match the stat-dual value at the initial time/state.
%
% functions called: Ntildeviafxdpt [ gradtheta, gettheta ]
%
%
% Extract data.
%
it=1;
alpha1=alphav(1,1)
%
% Get the correct value of alpha at the initial time.
%
%
[Ntilde,alphavfromx,thetafromx]=Ntildeviafxdpt(x1,c1,funcparams);
alphachk1=alphavfromx(1,1);
%
% If desired, we could check the resulting difference in the nonlinearity.
% Currently commented out.
%
%[thetachk]=gettheta(alpchk,funcparams)
%thetachk2=thetafromx
%[gradthetachk]=gradtheta(alpchk,funcparams);
%[gradthetachk2]=gradtheta(alpchk2,funcparams);
%chkN=invc1*gradthetachk-(alpchk-x1)
%chkN2=invc1*gradthetachk2-(alpchk2-x1)

