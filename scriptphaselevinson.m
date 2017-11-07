%%% Algorithme de Levinson-Durbin pour la d\'etermination des param\`etres AR
%%% Effet de la phase d'une sinuso\"{\i}de sur la DSP obtenue
%%%
%%% S. Rossignol -- 2012

clear;
close all;

fe=32000;   %% fr\'equence d'\'echantillonnage  (en Hz)
f0=440;     %% fr\'equence de la sinuso\"{\i}de (en Hz)
tsig=1280;  %% longueur du signal               (en nombre d'\'echantillons)
ordre=150;  %% ordre du mod\`ele AR

nphi=100;   %% nombre de phases consid\'er\'ees
phase=2*pi*[0:nphi-1]/nphi;

%%% boucle sur les phases
for ii=1:length(phase)
  
  xx=cos(2*pi*f0/fe*[1:tsig]+phase(ii));                        %%% construction du signal

  [aa, sigma2, ref, ff, mydsp] = mylevinsondurbin(xx,ordre,fe); %%% calcul du mod\`ele AR

  [mmax, pmax] = max(mydsp);
  rapport(ii) = mydsp(round(length(ff)-1)/2+1)/mmax;            %%% rapport entre l'amplitude \`a la
                                                                %%% fr\'equence 0 et l'amplitude du maximum
  position(ii)= abs(ff(pmax));                                  %%% position fr\'equentielle du max
                                                                %%% note : elle peut \^etre en -f0 ou en +f0,
                                                                %%%   d'o\`u l'abs

  figure(2);
  clf;
  semilogy(phase(1:ii),rapport,'linewidth',2);
  hold on;
  grid on;
  xlabel('phase (in rad)','fontsize',20);
  title('rapport dsp(f=0)/max(dsp)','fontsize',20);
  hold off;

  figure(3);
  clf;
  plot(phase(1:ii),position,'linewidth',2);
  hold on;
  grid on;
  xlabel('phase (in rad)','fontsize',20);
  title('frequence pour laquelle on a le max','fontsize',20);
  hold off;

  drawnow;
end;

fprintf(1, 'nombre de succes %d (sur %d)\n', sum(position>=400 & position<=480), nphi);

