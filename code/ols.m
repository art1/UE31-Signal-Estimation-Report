% function x = ols(H,z,mu,test);
% 
% Minimisation du crit�re J(x) = ||z-Hx||_2^2
% par l'algorithme OLS
%
% - Prend en compte une matrice H sparse
% - Argument optionnel qinit pour choisir une initialisation

% H. Carfantan d'après J. Idier, juin 2008
% Modif Herve par rapport a SMLR : 
% abs(.)^2 pour donnees z et matrice H complexes

function [x, vecm1]= ols(H,z,Nb_it,test);

VERB = 1;        

HH = H'*H;
Hz = H'*z;
h2 = diag(HH);
z2 = z'*z;
M = size(H,2);

% Solution initiale nulle
q = spalloc(M,1,M);
vecm1 = [];
lvecm0 = M;
lvecm1 = 0;
F = [];
crit = z2;

fprintf('0:\t\t%f\n',crit)

it=0;
while ((it<=Nb_it) && (crit>test))
  it=it+1;
  dcritopt = 0;
  mopt = [];

% Balayage: on envisage les ajouts
  vecm0 = find(q==0);
%  FGH = F*HH(vecm1,vecm0);
  for ind = 1:lvecm0,
    m = vecm0(ind);
    if(lvecm1)
      Gh = HH(vecm1,m);
      hGF = (F*Gh)';
      F22 = 1/(h2(m)-hGF*Gh);
      dcritnew = - F22*abs((hGF*Hz(vecm1)) - Hz(m))^2;
    else
      F22 = 1/h2(m);
      dcritnew = - F22*abs(Hz(m))^2;
    end
    if(dcritnew<dcritopt)
      mouv = 1;
      mopt = m;
      dcritopt = dcritnew;
    end
  end

      crit = crit + dcritopt;
      % Remise a jour de F
      Gh = HH(vecm1,mopt);
      FGh = F*Gh;
      F22 = 1/(h2(mopt)-Gh'*FGh);
      F12 = -FGh*F22;
      F11 = F-FGh*F12';
      F = [F11 F12;F12' F22];

      vecm1 = [vecm1; mopt];
      lvecm0 = lvecm0 - 1;
      lvecm1 = lvecm1 + 1;
      if(VERB)
        fprintf('%d:\t+%d\t%f\n',it,mopt,crit)
      end
    q(mopt) = ~q(mopt);

end
x = spalloc(M,1,M);
x(vecm1) = F*Hz(vecm1);


