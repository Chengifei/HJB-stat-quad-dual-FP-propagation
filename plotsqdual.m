% plotsqdual
%
%
xlo=-1.5;
xhi=1.5;
nx=80;
nxp1=nx+1;
xv=linspace(xlo,xhi,nxp1);
c1=-2.1;
twoc1=2*c1;
c1o2=c1/2;
alpl=-4;
alph=4;
nalp=2000;
nalpp1=nalp+1;
alpv=linspace(alpl,alph,nalpp1);
del=1e-6;
atdv=zeros(1,nalpp1);
foxv=zeros(1,nxp1);
for ix=1:nxp1
  x=xv(1,ix);
  for ia=1:nalpp1
    alp=alpv(1,ia);
    fm=funcforstatchk(alp-del);
    fp=funcforstatchk(alp+del);
    fderapp=(fp-fm)/(2*del);
    abstotder=abs(fderapp-c1*(alp-x));
    atdv(1,ia)=abstotder;
  end;
  [mintd,imin]=min(atdv)
  if (imin==1)||(imin==nalpp1)
      'argmin at an endpoint.'
      ['x= ',num2str(x)]
      'Pausing.'
      pause;
  end;
  alpmin=alpv(imin);
  fofx=funcforstatchk(alpmin)+c1o2*(alpmin-x)^2;
  foxv(1,ix)=fofx;
end;
figure;
plot(xv,foxv);

