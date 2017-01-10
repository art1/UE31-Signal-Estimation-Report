% Compilation du programme C/Matlab
%  mex min_L2_L1_0.c
%%

N = 100;
t = sort(rand(N,1));
y = cos(2*pi*10*t + pi/3) + .1*randn(N,1);

%%
K = 200;
fmax=20;
freq = (-K:K)/K*fmax;
W = exp(2*j*pi*t*freq);
lambda_max = max(abs(W'*y));
lambda = .3* lambda_max;

itmax = 50000;
% Appel du programme avec tous les arguments de sortie
  [X0, J0, NZ0] = min_L2_L1_0(y,W,lambda,itmax); 

% Appel du programme avec un seul argument de sortie
  X0 = min_L2_L1_0(y,W,lambda,itmax); 

