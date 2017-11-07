%%% Algorithme de Marple pour la d\'etermination des param\`etres AR
%%%
%%% entr\'ees :
%%%   - xx : signal
%%%   - pp : ordre du mod\`ele AR (choisi de mani\`ere ind\'ependante)
%%%   - fe : fr\'equence d'\'echantillonnage
%%%
%%% sorties :
%%%   - aa     : les param\`etres AR
%%%   - sigma2 : variance du bruit
%%%   (- ref    : coefficients de r\'eflexion -- ce n'est le cas pour l'instant)
%%%   - ff     : fr\'equences auxquelles la dsp a \'et\'e calcul\'ee
%%%   - mydsp  : la dsp elle-m\^eme
%%%
%%% exemple : fe=32000;f0=440;tsig=1280;xx=cos(2*pi*f0/fe*[1:tsig]+2*pi*rand(1,1))+sqrt(1e-4)*randn(1,tsig);pp=4;mymarple(xx,pp,fe);
%%%
%%% S. Rossignol -- 2012


function [aa, sigma2, kk, ff, mydsp] = mymarple_matlab(xx, pp, fe)


%%% Initialisations
NN = length(xx);
xx  = xx(:);
tol1=0.0;    %%% tol\'erance 1 (0.0 : pas de tol\'erance)
tol2=0.0;    %%% tol\'erance 2 (0.0 : pas de tol\'erance)
etat = 1;    %%% l'algorithme ne converge pas toujours
OrdreMax = pp;
aa(1) = 1.0;
f1=0.0;
for kk=1:NN
  f1 = f1 + real(xx(kk)*conj(xx(kk)));
end;
e0 = 2.0*f1;
Q1 = 1.0/e0;
Q2 = conj(Q1*xx(1));
ga = Q1*real(xx(1)*conj(xx(1)));
wa = Q1*real(xx(NN)*conj(xx(NN)));
den = 1.0 - ga - wa;
fa = xx(1);
ba = xx(NN);
ha = Q2*conj(xx(NN));
sa = Q2*xx(NN);
va = Q2*conj(xx(1));
ua = xx(NN)*xx(NN);
ua = Q1*ua;
Q4 = 1.0/den;
Q5 = 1.0 - ga;
Q6 = 1.0 - wa;
ea = e0 * den;
Q1 = 1.0/ea;
ca(1) = conj(Q1*xx(1));
da(1) = Q1*xx(NN);
MM=1;
c1 = 0.0;
for kk=2:NN
  c1 = c1 + xx(kk)*conj(xx(kk-1));
end;
c1 = 2.0*c1;
ra(1) = c1;
aa(2) = -Q1*ra(1);
f1 = real(aa(2)*conj(aa(2)));
ea = ea*(1.0-f1);

if ( f1 > 1.0 )
  etat=3;
end;

%%% d\'ebut de boucle principale
while ( (MM<OrdreMax) & (etat==1) )
   %%% mise \`a jour des erreurs de pr\'ediction
   evieux = ea;

   c1 = xx(MM+1);

   for kk=2:MM+1
     c1 = c1 + xx(MM-kk+2)*aa(kk);
   end;

   fa = c1;

   c1 = xx(NN-MM);
   for kk=1:MM
     c1 = c1 + xx(NN-MM+kk)*conj(aa(kk+1));
   end;

   ba = c1;

   %%% mise \`a jour des vecteurs auxiliaires
   Q1 = 1.0/ea;
   Q2 = conj(Q1*fa);
   Q3 =  Q1*ba;

   for kk=MM+1:-1:2
     ca(kk) = ca(kk-1) + Q2*aa(kk);
     da(kk) = da(kk-1) + Q3*aa(kk);
   end;

   ca(1) = Q2;
   da(1) = Q3;


   %%% mise \`a jour des scalaires auxiliaires
   Q7 = real(sa*conj(sa));

   f1 = real(fa*conj(fa));
   f2 = real(va*conj(va));
   f3 = real(ba*conj(ba));
   f4 = real(ua*conj(ua));

   c1 = conj(va)*ha*sa;
   ga = ga + f1*Q1 + Q4*( f2*Q6 + Q7*Q5 + 2.0*real(c1) );

   c1 = conj(sa)*ha*ua;
   wa = wa + f3*Q1 + Q4*( f4*Q5 + Q7*Q6 + 2.0*real(c1) );

   ha = 0.0;
   for kk=1:MM+1
     ha = ha + conj(xx(NN-MM+kk-1))*ca(kk);
   end;

   sa = 0.0;
   for kk=1:MM+1
     sa = sa + xx(NN-kk-1)*ca(kk);
   end;

   ua = 0.0;
   for kk=1:MM+1
     ua = ua + xx(NN-kk-1)*da(kk);
   end;

   va = 0.0;
   for kk=1:MM+1
     va = va + conj(xx(kk))*ca(kk);
   end;


   %%% mise \`a jour du d\'enominateur
   Q5  = 1.0 - ga;
   Q6  = 1.0 - wa;
   den = Q5*Q6 - real(ha*conj(ha));

   if ( den<=0.0 )
     etat=2;
     MM=MM+1;
     %fprintf(1,'\netat=2 : donnees numeriques mal conditionnees\n\n');
     break;
   end;


   %%% mise \`a jour des param\`etres de d\'ecalage temporel
   Q4 = 1.0/den;
   Q1 = Q1*Q4;
   c1 = ha*fa*ba;
   alpha = 1.0 / (1.0 + Q1*( f1*Q6 + f3*Q5 + 2.0*real(c1)));
   ea = alpha * ea;

   c1 = Q6*fa;
   c1 = c1 + conj(ba*ha);
   c1 = Q4*c1;

   c2 = conj(Q5*ba);
   c2 = c2 + ha*fa;
   c2 = Q4*c2;

   c3 = Q6*va;
   c3 = c3 + ha*sa;
   c3 = Q4*c3;

   c4 = Q5*sa;
   c4 = c4 + va*conj(ha);
   c4 = Q4*c4;

   c5 = Q6*sa;
   c5 = c5 + ha*ua;
   c5 = Q4*c5;

   c6 = Q5*ua;
   c6 = c6 + sa*conj(ha);
   c6 = Q4*c6;

   for kk=2:MM+1
      aa(kk) = aa(kk) + c1*ca(kk) + c2*da(kk);
      aa(kk) = alpha*aa(kk);
   end;

   M2 = MM/2 + 1;
   for kk=1:M2
      MK = MM + 2 - kk;

      save1 = conj(ca(kk));
      save2 = conj(da(kk));
      save3 = conj(ca(MK));
      save4 = conj(da(MK));

      ca(kk) = ca(kk) + c3*save3 + c4*save4;
      da(kk) = da(kk) + c5*save3 + c6*save4;

      if ( MK~=kk )
	 ca(MK) = ca(MK) + c3*save1 + c4*save2;
         da(MK) = da(MK) + c5*save1 + c6*save2;
      end;
   end;


   %%% mise \`a jour de l'ordre
   MM=MM+1;

   delta = 0.0;
   for kk=MM:-1:2
      c1     = ra(kk-1) - xx(NN-kk+2)*conj(xx(NN-MM+1));
      ra(kk) = c1       - xx(MM-1)*conj(xx(kk-1));

      delta = delta + ra(kk)*aa(kk);
   end;

   c1 = 0.0;
   for kk=2:(NN-MM+1)
      c1 = c1 + xx(kk+MM-1)*conj(xx(kk-1));
   end;

   ra(1) = 2.0*c1;

   delta = delta + ra(1);

   Q2 = -delta/ea;

   aa(MM+1) = Q2;

   M2 = floor(MM/2);
   for kk=1:M2
      MK = MM - kk;
      save1 = conj(aa(kk+1));
      aa(kk+1) = aa(kk+1) + Q2*conj(aa(MK+1));

      if ( MK~=kk )
	    aa(MK+1) = aa(MK+1) + Q2*save1;
      end;
   end;

   f1 = real(Q2*conj(Q2));
   ea = ea*( 1.0 - f1);    

   if ( f1>1.0 )
     etat=3;
   elseif ( (ea/e0)<tol1 )
     etat=4;
   elseif ( ((evieux-ea)/evieux)<tol2 )
     etat=5;
   end;
end; %%% fin de boucle

npp=MM;
if ( etat~=1 )
  npp=MM-1;
end;
%fprintf(1, '%d ',npp);
aa = aa(:).';

sigma2 = ea;


%%% densit\'e spectrale de puissance
interm2=-j*2*pi/fe*[1:npp];
df=0.9765625; %%% la dsp est calcul\'ee tous les df Hz
ff=-fe/2:df:fe/2;

interm3=interm2'*ff;
interm=1.+aa(2:npp+1)*exp(interm3);
mydsp = sigma2./(interm.*conj(interm));

% figure(1);
% clf;
% grid on;
% hold on;
% plot(ff,mydsp,'linewidth',2);
% xlabel('frequency (in Hz)','fontsize',20);
% ylabel('magnitude','fontsize',20);
% hold off;
% drawnow;

