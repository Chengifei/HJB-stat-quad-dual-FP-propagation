% trk4revricandstmnstep
%
% Code for testing rk4revricnstep
%
% Set time range.
%
t0=0;
tend=1;
n=20;
%
% Set problem data.
%
dim=2;
Sc1=2*eye(dim,dim);
A=[0 1; -1 0];
Ap=A';
Gammatilde=eye(dim,dim);
%
% Set terminal conditions.
%
Pend=eye(dim,dim);
Phibarend=eye(dim,dim);
Pstmend=[Pend;Phibarend];
%
% Propagate
%
[tv,Pstmv]=rk4revricandstmnstep(tend,Pstmend,t0,n,Sc1,A,Ap,Gammatilde);
%
% Backsubstitution check.
%
h=(tend-t0)/n;
invh=1/h;
chkPdiffv=zeros(1,n);
chkstmdiffv=zeros(1,n);
for it=1:n
    % RHS of all
    %[RHSPstm]=-revricfunc(Pstmv(:,:,it),Sc1,A,Gammatilde);
    [RHSPstm]=revricandstmfunc(Pstmv(:,:,it),Sc1,A,Ap,Gammatilde)
    % change in P adn comparison
    delP=invh*(Pstmv(1:dim,:,it+1)-Pstmv(1:dim,:,it));
    chkPdiffv(1,it)=norm(delP+RHSPstm(1:dim,:),'fro');
    % change in \Phibar (stm) and comparison
    delstm=invh*(Pstmv(dim+1:2*dim,:,it+1)-Pstmv(dim+1:2*dim,:,it));
    chkstmdiffv(1,it)=norm(delstm+RHSPstm(dim+1:2*dim,:),'fro');
end;
chkPdiffv
Pstmv(1:dim,:,n+1)
chkstmdiffv
Pstmv(dim+1:2*dim,:,n+1)