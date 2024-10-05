function [Pstm1rev,t1]=rk4revricandstm(Pstm0rev,t0,h,Sc1,A,Ap,Gammatilde)
%
% Objective:
%    This function performs one step in the basic fourth-order
%    Runge-Kutta method for an ODE.
%
% input variables:
%   Pstm0rev - state at the current time step.
%   t0 - the current time step.
%   h - time-step size.
%   Sc1,A,Gammatilde - parameters on RHS of the ODE.
%
% output variables:
%   Pstm1rev - state at the next time step.
%   t1 - the next time step.
%
% functions called:
%   revricandstmfunc - the rhs of
%                      \dot Pstmrev = f(t,Pstmrev,Sc1,A,Ap,Gammatilde).
%
c1=1/6; 
%
% This one is for the REVERSE-TIME DRE!
%
Pstm0=Pstm0rev;
%
% Get f(t_0,Pstm_0), K_1 and the next (t,Pstm) point (tt,Pstmt).
%
%Pd0=revricfunc(P0,Sc1,A,Gammatilde);
Pstmd0=revricandstmfunc(Pstm0,Sc1,A,Ap,Gammatilde);
K1=h*Pstmd0;
tt=t0+0.5*h;
Pstmt=Pstm0+0.5*K1;
%
% Get f(tt,Pstmt), K_2 (labelled Pstmd1) and the
% next (t,Pstm) point (tt,Pstmt).
%
%Pstmd1=revricfunc(Pt,Sc1,A,Gammatilde);
Pstmd1=revricandstmfunc(Pstmt,Sc1,A,Ap,Gammatilde);
K2=h*Pstmd1;
Pstmt=Pstm0+0.5*K2;
%
% Get f(tt,Pstmt), K_3 (labelled Pstmd2) and the
% next (t,Pstm) point (tt,Pstmt).
%
%Pstm2=revricfunc(Pt,Sc1,A,Gammatilde);
Pstmd2=revricandstmfunc(Pstmt,Sc1,A,Ap,Gammatilde);
K3=h*Pstmd2;
tt=t0+h;
Pstmt=Pstm0+K3;
%
% Get f(tt,Pstmt) and K_4 (labelled Pstmd3).
%
%Pstmd3=revricfunc(Pt,Sc1,A,Gammatilde);
Pstmd3=revricandstmfunc(Pstmt,Sc1,A,Ap,Gammatilde);
K4=h*Pstmd3;
%
% Get (t1,Pstm1), i.e., the next time and the next state value.
t1=tt;
Pstm1=Pstm0+c1*(K1+K4+2*(K2+K3));
%
% This one is for the REVERSE-TIME DRE!
%
Pstm1rev=Pstm1;

