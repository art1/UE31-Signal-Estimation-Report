clear all, clc, close all;

% Modeling the problem
L = 30;
C = 30;
l0 = floor(L/2);
c0 = floor(C/2);
sigma_l = 10;
sigma_c = 5;
alpha = 0.3;
n = 0.5;
a = 10;
s = 3;

% noise parameter
sigma = 0.8;

% creating the galaxy
lin = 1:L; 
col = 1:C;
[Col,Lin] = meshgrid(col,lin);
nu = [l0; c0; sigma_l; sigma_c; alpha; n];
Gal = Sersic(nu,Lin,Col);

% adding random noise to the initial Galaxy Image
D = s + a*Gal + (sigma*randn(L,C));
d = D(:);
D = reshape(d,L,C);

% plot the initial data sets
figure
subplot(2,2,1), imagesc(Gal), title('initial data set')
subplot(2,2,2), imagesc(D), title('initial set with added noise')
colorbar

% now we write matrix H
H = [ones(L*C,1), Gal(:)];
ML = inv(transpose(H)*H)*transpose(H)*d

% question 3, has not to be done here - theoretical question
real = [s a].'
bias = ML - real

% we create our p-vector that contains theta and nu parameters
% p_init = [s; a; nu];
p_init = [1;20;14;17;9;4;1;0.5];

% start optimisation using the crit_J cost function
options = optimset('fminsearch');
p_opt = fminsearch(@(p) crit_J(p,D), p_init,options)
Gal = Sersic(p_opt(3:end),Lin,Col);
D_opt = p_opt(1) + p_opt(2)*Gal ;

d = D_opt(:);
D_opt = reshape(d,L,C);

subplot(2,2,3), imagesc(D_opt), title('estimated data set')
colorbar

% calculate image residuals from estimated parameters and show histogram
Res = D - D_opt;
stdRes = std2(Res)
subplot(2,2,4)
imhist(Res)
ylim([0 max(Res(:))+10])
title(sprintf('Residual histogram, stdDev: %0.3f',stdRes))