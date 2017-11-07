%%% Algorithme de Pisarenko pour la d\'etermination des param\`etres de Pisarenko
%%%
%%% entr\'ees :
%%%   - xx : signal
%%%   - pp : ordre du mod\`ele (choisi de mani\`ere ind\'ependante)
%%%   - fe : fr\'equence d'\'echantillonnage
%%%   - rec: 0 : toutes les valeurs propres                           + DSP
%%%          1 : recursion pour avoir la valeur propre la plus petite + spectre de raies
%%%          2 : recursion pour avoir la valeur propre la plus petite + enveloppe interpol\'ee
%%%
%%% sorties :
%%%   - ff     : fr\'equences auxquelles la dsp a \'et\'e calcul\'ee
%%%   - mydsp  : la dsp elle-m\^eme
%%%
%%% exemple : fe=32000;f0=440;xx=cos(2*pi*f0/fe*[1:1280]+2*pi*rand(1,1));mypisarenko(xx,1,fe,1);
%%%
%%% S. Rossignol -- 2012



function [ff, mydsp] = mypisarenko(xx, pp, fe, rec)
%clear;rand('seed',100*sum(clock));fe=32000;f0=440;tsig=1280;xx=cos(2*pi*f0/fe*[1:tsig]+2*pi*rand(1,1));pp=2;


MM1 = 2*pp;
MM=MM1+1;
xx = xx-mean(xx);


%%% corr\'elation
NN=length(xx);%round(length(xx)/2);

acf = xcorr(xx(1:NN), MM1, 'biased');
lMM=length(acf);
rrr1 = acf(MM1+1:lMM)';
for ii=1:MM1
  rrr1 = [rrr1 acf(MM1+1-ii:lMM-ii)'];
end;
rrr1 = rrr1';

rrr = rrr1;
irrr = inv(rrr);


if (rec==1 || rec==2)
  %%% recursion pour trouver la plus petite valeur propre

  mu=1;
  muold=mu;
  aaa=rand(1,2*pp+1);

  kk=1;
  while kk==1
    mu = aaa*rrr*aaa'/(aaa*aaa');
    aaa = irrr*mu*aaa';
    aaa = aaa';

    if (abs(muold-mu)<1e-10)
      kk=0;
    end;
    muold=mu;
  end;
  if(mu<0)fprintf(1,'[mu=%e]',mu);fflush(1);end;
else
  %%% m\'ethode directe pour trouver toutes les valeurs propres

  [vec, lambda] = eig(rrr);
  lamb = diag(lambda);
  [vl,pl] = min(lamb);
  aaa = vec(:,pl);
  if(vl<0)fprintf(1,'[vl=%e]',vl);fflush(1);end;
end;


%%%

racines = roots(aaa);
frequ = log(racines)/j/2/pi*fe;
frequ = abs(frequ(1:2:end));


%%% densit\'e spectrale de puissance

df=0.9765625; %%% la dsp est calcul\'ee tous les df Hz
ff=-fe/2:df:fe/2;

musicway=1;
if (rec==1)
  musicway=0;
elseif (rec==2)
  musicway=2;
end;
if (musicway==1)
  %%% comme dans MUSIC

  mydsp=zeros(length(ff),1);
  deni = vec(:,pl)*conj(vec(:,pl))';
  for ii=1:length(ff)
    ee = cos(2*pi*ff(ii)*[-1:MM-2]/fe);   %%% attention : ici, [-1:MM-2] au lieu de [0:MM-1]
                                          %%%   chez MUSIC
    den = conj(ee)*deni*ee';
    mydsp(ii) = abs(1/den);
  end;
elseif (musicway==0)
  %%% standard Pisarenko

  for ii=1:pp
    for jj=1:pp
      coscos(ii,jj) = cos(2.*pi*ii*frequ(jj)/fe);
    end;
  end;

  alpha = inv(coscos)*acf(2*pp+1:3*pp)';
  amp = sqrt(2*abs(alpha));

  mydsp=zeros(length(ff),1)+1e-10;
  for ii=1:pp
    [mini, posi] = min(abs(ff-frequ(ii)));
    mydsp(posi) = amp(ii);

    [mini, posi] = min(abs(ff+frequ(ii)));
    mydsp(posi) = amp(ii);
  end;
elseif (musicway==2)
  %%% standard Pisarenko + interpolation

  for ii=1:pp
    for jj=1:pp
      coscos(ii,jj) = cos(2.*pi*ii*frequ(jj)/fe);
    end;
  end;

  alpha = inv(coscos)*acf(2*pp+1:3*pp)';
  amp = sqrt(2*abs(alpha));

  [freqs, pos] = sort([-frequ' frequ']);
  ampl1 = [amp' amp'];
  ampl = ampl1(pos);

  freqs = [min(ff) freqs max(ff)];
  ampl = [0 ampl 0];

  mydsp = interp1 (freqs, ampl, ff, 'linear');
end;

% figure(1);
% clf;
% grid on;
% hold on;
% %semilogy(ff,mydsp,'linewidth',2);
% plot(ff,mydsp,'linewidth',2);
% xlabel('frequency (in Hz)','fontsize',20);
% ylabel('magnitude','fontsize',20);
% %xlim([-1000 1000]);
% hold off;
% drawnow;

