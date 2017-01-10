clear; clc; clear all;
load('data.mat')

x = 0;
x_n = 0;
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
plot(t,x,'+')
hold on;
plot(t,x_n,'+')


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

% 3
% Gamma is a set, Gamma = [Gamma, index]
k = 0;
r_n = x;
Gamma0 = [];
a = 0;

%[val, index] = max(abs(W'*r_n))
tau = chisqq(0.95,N)
T = 0;

while T < tau
    [val, index] = max(abs(W'*r_n));
    Gamma0 = [Gamma0 index];
    a = a + (1./(W'*W))*W'*r_n;
    r_n = r_n - ((1./(W'*W))*W'*r_n).'*W.'
end

