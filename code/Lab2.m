close all; clc; clear all;
load('data.mat')

x = 0;
SNR = 10;

for k=1 : 5
   x = x + (A_th(k) * sin(2*pi*f_th(k)*t + phi_th(k)))
end

% calculate noise amplitude
pms = sumsqr(x)/length(x);
% standard Deviation
sigma = sqrt(pms/10);

noise =  sigma * randn(1,length(x));
x_n = x + transpose(noise)

figure
plot(t,x,'g+')
hold on;
plot(t,x_n,'r+')


% irregular sampling case
fmax = 100;
M = 1024;
N = length(x_n);
freq =(-M:M)/M*fmax;
W=exp(2*j*pi*t*freq);
periodogram =  abs(W'*x_n)/N;

figure
plot(freq,periodogram)
% plot real frequencies
hold on
plot(f_th,A_th/2,'+')
hold on
plot((-1.*f_th),A_th/2,'+')


Win=W'*ones(N,1)/N;

figure
plot(freq,abs(Win))

% 3.1 Matching Pursuit Algorithm
r_n = x_n;
Gamma0 = [];
a = zeros(2049,1);
tau = chisqq(0.95,N)
T = tau +1;
k = 1;

while T > tau
    W_current = W(:,k);
    [val, k] = max(abs(W'*r_n));
    Gamma0 = [Gamma0 k];
    a(k) = a(k) + (1/((W_current')*W_current))*W_current'*r_n;
    r_n = r_n - (((1/(W_current'*W_current)).*W_current'*r_n).*W_current);
    T = (norm(r_n)^2)/(sigma^2);
end

% figure, plot(freq,abs(a));
% figure
% plot(t,W*a,'g+')
% hold on;
% plot(t,x_n,'r+')


MethodOneIterations = size(Gamma0);
%ylim([-4 4])


% 3.2 Orthogonal Matching Pursuit Algorithm

r_n = x_n;
Gamma0 = [];
W_g = [];
a = [];
tau = chisqq(0.95,N)
T = tau +1;
k = 1;

while T > tau
    [val, k] = max(abs(W'*r_n))
    Gamma0 = [Gamma0 k];
    
    W_g = [];
    for l=1:length(Gamma0)
        W_g = [W_g W(:,Gamma0(l))];
    end
    
    a = ((W_g'*W_g)^(-1))*W_g'*x_n;
    %size(a)
    %a_vec = [a_vec a];
    
    r_n = x_n - W_g*a;
    T = (norm(r_n)^2)/(sigma^2)
end
a_plot = zeros(2049,1);
for ind=1:length(Gamma0)
   a_plot(Gamma0(ind)) = a(ind);
end

%figure, plot(freq,abs(a_plot));
MethodTwoIterations = size(Gamma0);

% 3.3 Orthogonal Least Square
r_n = x_n;
Gamma0 = [];
W_g = [];
a = 0;
tau = chisqq(0.95,N)
T = tau +1;
k = 1;
test = (sigma^2)*tau;
a_vec = [];
while T > tau
    [val, k] = ols(W,x_n,Inf,test);
    Gamma0 = [Gamma0 k];
    
    W_g = [];
    for l=1:length(Gamma0)
        W_g = [W_g W(:,Gamma0(l))];
    end
    
    a = ((W_g'*W_g)^(-1))*W_g'*x_n;
    a_vec = [a_vec a];
    
    r_n = x_n - W_g*a;
    T = (norm(r_n)^2)/(sigma^2)
end

a_plot = zeros(2049,1);
for ind=1:length(Gamma0)
   a_plot(Gamma0(ind)) = a_vec(ind);
end

% figure, plot(freq,abs(a_plot))
% figure
% plot(t,W_g*a,'g+')
% hold on;
% plot(t,x_n,'r+')
MethodThrIterations = size(Gamma0);


% 4 Sparse representation with convex relations

lambda_max = max(abs(W'*x_n));
lambda = 0.05 * lambda_max;
n_it_max = 100000;

a1 = min_L2_L1_0(x_n,W,lambda,n_it_max);
figure, plot(t,W*a1,'+')
figure, plot (freq,abs(a1));
MethodFourIterations = size(Gamma0);