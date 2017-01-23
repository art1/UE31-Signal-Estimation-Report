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
legend('real data',sprintf('data with gaussian noise, stdDev: %0.3f',sigma))
title('Herbig star HD 104237 data samples')
ylabel('signal amplitude')
xlabel('time')
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPosition', [0 0 900 450]);
saveas(gcf,'../images/data.png')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.2  irregular sampling case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fmax = 100;
M = 1024;
N = length(x_n);
freq =(-M:M)/M*fmax;
W=exp(2*j*pi*t*freq);
periodogram =  abs(W'*x_n)/N;

figure
subplot(2,2,1) 
plot(freq,periodogram)
% plot real frequencies
hold on
plot(f_th,A_th/2,'r+')
hold on
plot((-1.*f_th),A_th/2,'r+')
xlabel('frequency')
ylabel('amplitude')
title('frequency representation')

subplot(2,2,3)
plot(freq,periodogram)
hold on;
plot(f_th,A_th/2,'r+')
xlim([28 40])
xlabel('frequency')
ylabel('amplitude')
legend('frequency distribution','original frequencies')

% plot spectral window
subplot(2,2,[2 4])
Win=W'*ones(N,1)/N;
plot(freq,abs(Win))
xlim([-20 20])
title('spectral window')
xlabel('frequency')
ylabel('amplitude')

set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPosition', [0 0 900 450]);
saveas(gcf,'../images/data_freq.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.1 Matching Pursuit Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
MethodOneIterations = size(Gamma0);
[MaxMP,MaxIdxMP] = findpeaks(abs(a),'threshold',0.13);

figure
subplot(2,1,1)
plot(freq,abs(a));
hold on;
plot(f_th,A_th/2,'r+')
hold on;
plot((freq(MaxIdxMP)),MaxMP,'bo')
for num=1:length(f_th)
    hold on;
   line([f_th(num) f_th(num)], [0 A_th(num)/2], 'Color','r','LineStyle','--') 
end
xlim([28 40])
legend('reconstructed frequency','original frequencies')
xlabel('frequency')
ylabel('amplitude')
subplot(2,1,2)
plot(t,W*a,'+')
suptitle('Matching Pursuit (pre-whitening) algorithm')
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPosition', [0 0 900 450]);
saveas(gcf,'../images/mp.png')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.2 Orthogonal Matching Pursuit Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
MethodTwoIterations = size(Gamma0);
[MaxOMP,MaxIdxOMP] = findpeaks(abs(a_plot),'threshold',0.1);

figure
subplot(2,1,1)
plot(freq,abs(a_plot));
hold on;
plot(f_th,A_th/2,'r+')
hold on;
plot((freq(MaxIdxOMP)),MaxOMP,'bo')
for num=1:length(f_th)
    hold on;
   line([f_th(num) f_th(num)], [0 A_th(num)/2], 'Color','r','LineStyle','--') 
end
xlim([28 40])
legend('reconstructed frequency','original frequencies')
xlabel('frequency')
ylabel('amplitude')
subplot(2,1,2)
plot(t,W*a_plot,'+')
xlabel('time')
ylabel('amplitude')
suptitle('Orthogonal Matching Pursuit algorithm');
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPosition', [0 0 900 450]);
saveas(gcf,'../images/omp.png')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.3 Orthogonal Least Square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
MethodThrIterations = size(Gamma0);
[MaxOLS,MaxIdxOLS] = findpeaks(abs(a_plot),'MinPeakHeight',0.1);

figure
subplot(2,1,1)
plot(freq,abs(a_plot));
hold on;
plot(f_th,A_th/2,'r+')
hold on;
plot((freq(MaxIdxOLS)),MaxOLS,'bo')
for num=1:length(f_th)
    hold on;
   line([f_th(num) f_th(num)], [0 A_th(num)/2], 'Color','r','LineStyle','--') 
end
xlim([28 40])
xlabel('frequency')
ylabel('amplitude')
legend('reconstructed frequency','original frequencies')
subplot(2,1,2)
plot(t,W*a_plot,'+')
xlabel('time')
ylabel('amplitude')
suptitle('Orthogonal Least Square');
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPosition', [0 0 900 450]);
saveas(gcf,'../images/ols.png')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 Sparse representation with convex relation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda_max = max(abs(W'*x_n));
lambda = 0.06 * lambda_max;
n_it_max = 100000;

a1 = min_L2_L1_0(x_n,W,lambda,n_it_max);

[MaxSparse,MaxIdxSparse] = findpeaks(abs(a1),'MinPeakHeight',0.03);

figure
subplot(2,1,1)
plot(freq,abs(a1));
hold on;
plot(f_th,A_th/2,'r+')
hold on;
plot((freq(MaxIdxSparse)),MaxSparse,'bo')
for num=1:length(f_th)
   hold on;
   line([f_th(num) f_th(num)], [0 A_th(num)/2], 'Color','r','LineStyle','--') 
end
xlim([28 40])
xlabel('frequency')
ylabel('amplitude')
legend('reconstructed frequency','original frequencies')
subplot(2,1,2)
plot(t,W*a1,'+')
xlabel('time')
ylabel('amplitude')
suptitle('Sparse representation with convex relaxation');
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPosition', [0 0 900 450]);
saveas(gcf,'../images/convex.png')

% Save detected data to file:
fileID = fopen('../images/img_data.txt','w');
fprintf(fileID, 'frequency ; amplitude\n ');
fprintf(fileID, '# Matching Pursuit, %d iterations\n', length(MethodOneIterations(1)));
fprintf(fileID, [repmat(' %0.3f ; %0.3f \n', 1, length(MaxIdxMP))] , [ freq(MaxIdxMP).' MaxMP].');
fprintf(fileID, '# Orthogonal Matching Pursuit, %d iterations\n', length(MethodTwoIterations(2)));
fprintf(fileID, [repmat(' %0.3f ; %0.3f \n', 1, length(MaxIdxOMP))] , [ freq(MaxIdxOMP).' MaxOMP].');
fprintf(fileID, '# Orthogonal Least Square, %d iterations\n', length(MethodThrIterations(1)));
fprintf(fileID, [repmat(' %0.3f ; %0.3f \n', 1, length(MaxIdxOLS))] , [ freq(MaxIdxOLS).' MaxOLS].');
fprintf(fileID, '# Convex Relaxation\n');
fprintf(fileID, [repmat(' %0.3f ; %0.3f \n', 1, length(MaxIdxSparse))] , [ freq(MaxIdxSparse).' MaxSparse].');

fclose(fileID);

