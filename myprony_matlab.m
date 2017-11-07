%%% Algorithme de Prony historique pour la d\'etermination des param\`etres
%%%
%%% entr\'ees :
%%%   - xx : signal
%%%   - pp : ordre du mod\`ele
%%%       choisi de mani\`ere ind\'ependante ou \'egal \`a la demi-longueur du signal
%%%   - fe : fr\'equence d'\'echantillonnage
%%%   - rec: 0 : DSP
%%%          1 : spectre de raies
%%%          2 : spectre de raies + interpolation
%%%
%%% sorties :
%%%   - aa     : les param\`etres de Prony
%%%   - ff     : fr\'equences auxquelles la dsp a \'et\'e calcul\'ee
%%%   - mydsp  : la dsp elle-m\^eme
%%%
%%% exemple : fe=32000;f0=440;tsig=1280;xx=cos(2*pi*f0/fe*[1:tsig]+2*pi*rand(1,1))+sqrt(1e-2)*randn(1,tsig);myprony(xx,40,10,fe,1);
%%%
%%% S. Rossignol -- 2012


function [aa, ff, mydsp] = myprony_matlab(xx, pp, nco, fe, rec)

%clear;rand('seed',100*sum(clock));fe=32000;f0=440;tsig=1280;xx=sqrt(1e-2)*randn(1,tsig)+cos(2*pi*f0/fe*[1:tsig]+2*pi*rand(1,1))+cos(2*pi*(f0+26)/fe*[1:tsig]+2*pi*rand(1,1));pp=400;
%clear;rand('seed',100*sum(clock));fe=32000;f0=440;tsig=1280;xx=sqrt(1e-2)*randn(1,tsig);for(ii=1:10)  xx=xx+cos(2*pi*ii*f0/fe*[1:tsig]+2*pi*rand(1,1));end;pp=10;


tsig = length(xx);

df=0.9765625; %%% la dsp est calcul\'ee tous les df Hz
ff=-fe/2:df:fe/2;


strt=0;
count=0;
done=0;
while (done<1)
  ok1=0;
  ok2=0;

  %%% alpha et f

  yy=[];
  for (ii=1:pp)
    yy = [yy xx(pp+ii-1+strt:-1:ii+strt)'];
  end;
  yy=yy';

  y1 = xx(pp+1+strt:2*pp+strt)';

  condit=rcond(yy);
  if (condit>1.52721e-17)
    invaa=inv(yy);
    ok1=1;
    aa = -invaa*y1;
    aa = [1 aa']';

    racines=roots(aa)';

    alpha = log10(abs(racines));
    frequ = 1/2/pi*atan2(imag(racines),real(racines));


    %%% rho et theta

    zz=[];
    for (ii=0:pp-1)
      zz = [zz racines.^ii'];
    end;
    zz=zz';

    y2 = xx(1:pp)';   %%% 'strt' is not used here; should it be changed?
    
    condit=rcond(zz);
    if (condit>1.52721e-17)
      invzz=inv(zz);
      ok2=1;
    end;
  end;

  if (ok1==1 & ok2==1)
    hh = invzz*y2;

    rho = abs(hh);
    theta = atan2(imag(hh),real(hh));

    done=1;
  else
    fprintf(1,'pc ');
    count=count+1;
    if (count==nco)
      done=2;
    end;
    strt=round(rand(1,1)*(tsig-2*pp));
  end;
end;


%%% densit\'e spectrale de puissance

ff=ff/fe;
if (done==1)
  if (rec==0)
    %%% DSP Prony
    sff=zeros(1,length(ff));
    for ii=1:pp
      sff = sff + hh(ii) *( 1./(1 - racines(ii)*exp(-j*2*pi*ff) ) - 1./(1 - ( conj(racines(ii))*exp(j*2*pi*ff) ).^(-1) ) );
    end;
    mydsp=abs(sff).^2/fe;
  elseif (rec==1)
    %%% spectre de raies
    mydsp=zeros(1,length(ff))+1e-10;
    for ii=1:pp
      [mini, posi] = min(abs(ff-frequ(ii)));
      mydsp(posi) = abs(hh(ii));

      [mini, posi] = min(abs(ff+frequ(ii)));
      mydsp(posi) = abs(hh(ii));
    end;
  else
    %%% spectre de raies interpol\'e
    [freqs, pos] = sort([-frequ' frequ']);
    ampl1 = abs([hh' hh']);
    ampl = ampl1(pos);

    freqs = [min(ff) freqs max(ff)];
    ampl = [0 ampl 0];

    mydsp = interp1 (freqs, ampl, ff, 'linear');
  end;
else
  fprintf(1,'unable to comply\n');
  mydsp=zeros(1,length(ff))+1e-15;
  mydsp(round(length(ff)/2))=1e-5;
  aa=[];
end;
ff=ff*fe;

figure(1);
clf;
grid on;
hold on;
%semilogy(ff,mydsp,'linewidth',2);
plot(ff,mydsp,'linewidth',2);
xlabel('frequency (in Hz)','fontsize',20);
ylabel('magnitude','fontsize',20);
hold off;

figure(2);
clf;
grid on;
hold on;
%semilogy(ff,mydsp,'linewidth',2);
plot(ff,mydsp,'linewidth',2);
xlabel('frequency (in Hz)','fontsize',20);
ylabel('magnitude','fontsize',20);
xlim([400 506]);
hold off;
drawnow;

