%%% Algorithme de Burg pour la d\'etermination des param\`etres AR
%%%
%%% entr\'ees :
%%%   - xx : signal
%%%   - pp : ordre du mod\`ele AR (choisi de mani\`ere ind\'ependante)
%%%   - fe : fr\'equence d'\'echantillonnage
%%%
%%% sorties :
%%%   - aa     : les param\`etres AR
%%%   - sigma2 : variance du bruit
%%%   - ref    : coefficients de r\'eflexion
%%%   - ff     : fr\'equences auxquelles la dsp a \'et\'e calcul\'ee
%%%   - mydsp  : la dsp elle-m\^eme
%%%
%%% exemple : fe=32000;f0=440;xx=cos(2*pi*f0/fe*[1:1280]+2*pi*rand(1,1));myburg(xx,4,fe);
%%%
%%% S. Rossignol -- 2012

function [aa, sigma2, kk, ff, mydsp] = myburg(xx, pp, fe)

%%% Initialisations
NN  = length(xx);
xx  = xx(:);
ef = xx;
eb = xx;
aa = 1;

%%% Erreur initiale
EE = xx'*xx./NN;

%%% Pr\'e-allocation de 'kk' (pour am\'eliorer la vitesse d'ex\'ecution)
kk = zeros(1, pp);

for tt=1:pp
   %%% Calcule le coefficient de r\'eflexion suivant
   efp = ef(2:end);
   ebp = eb(1:end-1);
   num = -2.*ebp'*efp;       %%% c'est -Nk (F5)
   den = efp'*efp+ebp'*ebp;  %%% c'est Dk  (F6)

   kk(1,tt) = num./den;      %%% c'est -gamma dans le cours (F7)

   %%% Mise \`a jour des erreurs de pr\'ediction 'avant' et 'arri\`ere'
   ef = efp + kk(tt)*ebp;    %%% F3
   eb = ebp + kk(tt)'*efp;   %%% F4

   %%% Mise \`a jour des coefficients
   aa=[aa;0] + kk(1,tt)*[0;conj(flipud(aa))];   %%% F1 et F2

   %%% Mise \`a jour de l'erreur de pr\'ediction globale ; E5
   EE(tt+1) = (1 - kk(1,tt)'*kk(1,tt))*EE(tt);
end
aa = aa(:).';

%%% pas forc\'ement bien estim\'e ; note : c'est EE
%sigma2=1;
sigma2=EE(end);

%%% densit\'e spectrale de puissance
interm2=-j*2*pi/fe*[1:pp];
df=0.9765625; %%% ls dsp est calcul\'ee tous les df Hz
ff=-fe/2:df:fe/2;

interm3=interm2'*ff;
interm=1.+aa(2:pp+1)*exp(interm3);
mydsp = sigma2./(interm.*conj(interm));

% figure(1);
% clf;
% grid on;
% hold on;
% plot(ff,mydsp,'linewidth',2);
% xlabel('frequency (in Hz)','fontsize',20);
% ylabel('magnitude','fontsize',20);
% xlim([-1000 1000]);
% hold off;
% drawnow;

