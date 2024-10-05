function [Watptv,relbsalpproperr]=Watoneypoint(y,tfullv,M0,Abarv,...
    Pv,Phibarst,DMz,...
    y1,c1,Gammatilde,nt1,funcparams, chkqandr)
if ~exist("chkqandr", "var")
    chkqandr = false;
end
%
% Obtains W at a single point in space, for a vector of times.
%
%
% Extract data.
%
ntotp1=length(tfullv);
%
% Get the correct alpha0v vector.
%
%y1=y(1,1);
[alpha0v]=Getalpha0psiTbar(M0,y,tfullv,Phibarst);
%
% Call propalpnlpsiTbar.
%
[allalpham,t1vfinal,relbsalpproperr]=...
    propalpnlpsiTbndim(M0,Abarv,Phibarst,...
    DMz,alpha0v,c1,tfullv,nt1,funcparams);
numalpfuncs=size(allalpham,1);
alphafinalv=allalpham(numalpfuncs,:);
%
% Check whether this is really stat or not, if desired.
%
checkingstat=0;
if checkingstat==1
  [effectofhalving,shouldbeahalf]=checkstat(alphafinalv,y,tfullv,Phibarst,...
     M0,c1,Gammatilde,funcparams);
  ['Checking if alpha really stat at y= ', num2str(y')]
  ['The effect of halving (should be near 0.25) is ',num2str(effectofhalving)]
  ['Check on halving (should be near 0.5) is ',num2str(shouldbeahalf)]
  %'Pausing in Watoneypoint'
  %pause;
end;
%
%
% Get q_s and r_s for the staticizing alpha
%
[qs,rs,qv,rv]=qrintegrals(t1vfinal,alphafinalv,Phibarst,...
    M0,c1,Gammatilde,funcparams);
%
% Compute W(t1v(is),y)=0.5*y'*Pv(:,:,is)*y+y'*qv(1,is)+0.5*rv(1,is)
%
Watptv=zeros(1,ntotp1);
for it=1:ntotp1
  Watptv(1,it)=0.5*y'*Pv(:,:,it)*y+y'*qv(:,it)+0.5*rv(1,it);
end;


if chkqandr==1
    [qchkv,rchkv,relqchkv,relrchkv]=Checkqandr(qv,rv,t1vfinal,Pv,alphafinalv,...
        A,C,M0,Gammatilde,c1,Phibarst,Abarv,funcparams);
    'Checking q and r. Pausing.'
    relqchkv
    relrchkv
    pause;
end
%size(Watptv)
