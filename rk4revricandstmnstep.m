function [tv,Pstmv]=rk4revricandstmnstep(tend,Pstmend,t0,n,Sc1,A,Ap,Gammatilde)
%
% Objective: Perform 4th-order Runge-Kutta method for solution of
% % the REVERSE-TIME DRE
% \dot Prev = - [Prev*Gammatilde*Prev - (Sc1 + A'*Prev + Prev*A) ]
% along with the RHS of the state-transition matrix ODE
% for $\dot\Phibar = - \bar A \Phibar$ (this is reverse time!) i.e., 
% \dot Phibarrev = - [ Prev * Gammatilde - A'] Phibarrev.
% Both are over time interval [tend,t0]
% with terminal data Pstmend=[Pend;Phibarend].
%
% NOTE: The lower half of the output is \Psi_t matrices
% - NOT(!) \Phi_{s,t} state-transition matrices. You can obtain
% the latter by \Phi_{s,t}=\Psi_s[\Psi_t]^{-1}.
%
% input variables:
%   Pstmrev - terminal state (initial in reverse time) of
%               [ Prev ; Phibarrev]
%
% input variables:
%   tend - terminal time
%   Pstmend - terminal state. This is a (2*dim)\times(dim) matrix,
%             where dim is the state-space dimension.
%   t0 - initial time (that we propagate down to)
%   n - number of steps.
%   Sc1 - S_{c_1} see paper
%   A - A see paper
%   Ap - A'
%   Gammatilde - \tilde\Gamma see paper
%
% output variables:
%   tv - vector of times at which the solution is evaluated (including t0
%        and tend). This will be an 1 by n+1 array.
%        t0 is the first element and tend is the last element.
%   Pstmv - "vector" of the state values, [P(t);\bar\Phi(t)],
%           at the times in tv.
%           [P(t0);\bar\Phi(t0)] is the FIRST element,
%           and Pstmend=[P(tend);\bar\Phi(tend)] is the LAST element.
%
% functions called:
%    rk4revricandstm - one-step in RK4 [calls revricandstmfunc].
%
%
% Get array sizes and allocate memory.
%
np1=n+1;
dim=size(Pstmend,2);
twodim=2*dim;
Pstmv=zeros(twodim,dim,np1);
Pstmrevv=zeros(twodim,dim,np1);
%
% Initialize (reverse) time and (reverse-time) state.
%
t=0;
Pstmrev0=Pstmend;
Pstm=Pstmrev0;
%
% Put terminal state into the output state vector.
%
Pstmrevv(:,:,1)=Pstm;
%
% Compute the step-size, h.
%
h=(tend-t0)/n;
%
% Begin RK45 iteration (in reverse time).
%
for irk4=2:np1
    %
    % Update time and state.
    %
    %[P1,t1]=rk4revric(P,t,h,Sc1,A,Gammatilde);
    [Pstm1,t1]=rk4revricandstm(Pstm,t0,h,Sc1,A,Ap,Gammatilde);
    t=t1;
    Pstm=Pstm1;
    %
    % Load new reverse time and reverse-time state into the output vectors.
    %
    Pstmrevv(:,:,irk4)=Pstm;
end;
%
% Fix arrow of time.
% 
%
np2=np1+1;
for it=1:np1
    itrev=np1-it;
    Pstmv(:,:,it)=Pstmrevv(:,:,np2-it);
end;
tv=linspace(t0,tend,np1);
