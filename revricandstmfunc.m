function [Pstmrevdot]=revricandstmfunc(Pstmrev,Sc1,A,Ap,Gammatilde)
%
% Obtains the RHS of the DRE for reverse-time P (Prev),
% i.e., the RHS of
% \dot Prev = - [Prev*Gammatilde*Prev - (Sc1 + A'*Prev + Prev*A) ]
% along with the RHS of the state-transition matrix ODE
% for $\dot\Phibar = - \bar A \Phibar$ (this is reverse time!) i.e., 
% \dot Phibarrev = - [ Prev * Gammatilde - A'] Phibarrev.
%
% input variables:
%   Pstmrev - terminal state (initial in reverse time) of
%               [ Prev ; Phibarrev]
%   Sc1 - S_{c_1} see paper
%   A - A see paper
%   Ap - A'
%   Gammatilde - \tilde\Gamma see paper
%
% output variables:
%   Pstmrevdot - reverse-time derivative of [ Prev ; Phibarrev]
%
%
% Get dimensions
%
dim=size(Pstmrev,2);
%
% Get RHS of ODEs
%
Prev=Pstmrev(1:dim,:);
Phibarrev=Pstmrev(dim+1:2*dim,:);
Prevdot=Sc1+Ap*Prev+Prev*A-Prev*Gammatilde*Prev;
Phibarrevdot=(Ap-Prev*Gammatilde)*Phibarrev;
Pstmrevdot=[Prevdot;Phibarrevdot];
