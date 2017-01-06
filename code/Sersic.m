function Im = Sersic(nu,Lig,Col)

% Im = Sersic(nu,Lig,Col)
%
% Calcul de l'image d'un modèle de Sersic sur une grille (définie par les
% matrices Lig et Col) de paramètres contenus dans le vecteur nu:
%  - nu(1) = l_0 : position en ligne
%  - nu(2) = c_0 : position en colonne
%  - nu(3) = sigma_l : écart-type en ligne
%  - nu(4) = sigma_c : écart-type en colonnne 
%  - nu(5) = angle : angle en radian   
%  - nu(6) =  n : indice de Sersic

% H. Carfantan, IRAP, novembre 2014

l_0 = nu(1);
c_0 = nu(2);
sigma_l=nu(3);
sigma_c=nu(4);
n = nu(6);
alpha = nu(5);
Im = exp(-  (  (((Col-c_0)*cos(alpha)-(Lig-l_0)*sin(alpha))/sigma_l).^2 ...
             + (((Col-c_0)*sin(alpha)+(Lig-l_0)*cos(alpha))/sigma_c).^2).^(1/(2*n)) );

