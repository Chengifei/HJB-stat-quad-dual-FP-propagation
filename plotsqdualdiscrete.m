% plotsqdualdiscrete
%
%
xlo=-1.5;
xhi=1.5;
nx=200;
nxp1=nx+1;
xv=linspace(xlo,xhi,nxp1);
c1=-2.1;
twoc1=2*c1;
c1o2=c1/2;
del=1e-6;
alpl=-4;
alph=4;
nalpv=[25 100 400 1600]
nnalp=size(nalpv,2);
colorsv='brgkcobrgkco';
for inum=1:nnalp
  nalp=nalpv(1,inum);
  nalpp1=nalp+1;
  alpv=linspace(alpl,alph,nalpp1);
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
  %if inum==1
    figure;
    plot(xv,foxv,colorsv(1,inum));
    titletextv=['Discrete stat-quad dual with ',num2str(nalpv(1,inum)),' functions.']
    title(titletextv);
  %  hold on;
  %elseif in==nnalp
  %  plot(xv,foxv);
  %  hold off;
  %else
  %  plot(xv,foxv);
  %end
end;