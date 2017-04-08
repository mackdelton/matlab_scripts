%%
%%% As taken from MULTI-ARMED BANDIT ALLOCATION INDICES GITTINS ET AL.
%%% Initialize %%%
clear all;close all;clc;
% Set parameters for search routine
beta = 0.95;    %Discount rate (=a in text)
N = 750;        %Starting value for a+b *alpha + beta in text)
step = 0.0001;   %step size

% Initialize matrices
G = zeros(N-1,N-1);         %Array with intermediate values
Gittins = zeros(N-1,N-1);   %array with final Gittins indices

% Initialize endpoints (starting points for backward induction)

for a = 1:N-1
    G(a,N-a) = a/N;     %Initialize endpoints at E[beta(a,b)]
end

%%%%%%%%%%%%%%%%%
%%% Main Loop %%%
%%%%%%%%%%%%%%%%%
h = waitbar(0, 'Calculating Gittins'' Indices');
for p = step/2:step:1
    waitbar(p/1,h);
% Continuous value for 'safe' arm

    safe = p/(1-beta);
    
    for nn = N-1:-1:2
        for a = 1:nn-1
            
            % Continuation value for 'risky' arm
            risky = (a/nn) * (1 + beta*G(a+1,nn-a)) + (nn-a)/nn *(beta*G(a,nn-a+1));
            
            %'safe' increases faster than 'risky' in p, so
            % the first time safe = risky, we get the index value
            if (Gittins(a,nn-a) == 0) && (safe > risky)
                Gittins(a,nn-a) = p - step/2;
            end
            
            % Update recursion
            G(a,nn-a) = max(safe,risky);
            
        end
    end
end
close(h)
%%%%%%%%%%%%%%%%%%%%%%
%%% Output Results %%%
%%%%%%%%%%%%%%%%%%%%%%

% Set output ranges
a_out = 1:40;
b_out = 1:40;

Results = double([0 a_out ; b_out' Gittins(a_out,b_out)']);

%%

%Identifying the trend used in the Gittins' Index table of values
%(reference [25] from M. Cheung's MIT MS thesis.

%n = [1:9,10:10:90,100:100:1000];
n = [10:100];
nLen = length(n);
%Use the following values for a = 0.5:0.1:0.9,0.95,0.99,0.995
%a = [0.5:0.1:0.9,0.95,0.99,0.995];
a = [0.5 0.6 0.8 0.9 0.995]; %trial of different beta values
% From Table 8.1 of Gittins et al. (1989)
% Use values of beta == 0.95 for n = [100 1000]
v = [0.22263 0.28366 0.32072 0.34687 0.36678 0.38267 0.39577 0.40682 0.41631... %1-9
     0.42458 0.47295 0.49583 0.50953 0.51876 0.52543 0.53050 0.53449 0.53771... %10-90
    0.54037 0.55344 0.55829 0.56084 0.56242 0.56351 0.56431 0.56493 0.56543 0.56583]; %100 - 1000
vi = [1:10,20:10:100,200:100:1000]; %Generate vector of iterations associated with pre-calculated Gittins' Indices

vi_interp = linspace(1,1000,1000);
%vi_interp = vi_interp(vi_interp ~= vi);

vNew = interp1(vi,v,vi_interp,'cubic');
% %Plot of Gittins' Index values
%figure;
%plot(vi,v,'r');hold on
%plot(vi_interp,vNew,'b');

%Plot Gittins Index values (scaled version)
% Presuming a = 0.95, 0.2236 = (1-0.95)^0.5 and vi_interp = n
vNewScaled = vNew./(vi_interp*0.2236);
figure;
semilogy(vi_interp,vNewScaled,'g','linewidth',3);grid on;
xlabel('Number of Trials','fontsize',12,'fontweight','bold')
ylabel('Gittins'' Index','fontsize',12,'fontweight','bold')
set(gca,'fontsize',12,'fontweight','bold')
%Generate Gittins' Index reference for use in MAB multi-agent sim.
save('gittins.mat','vi_interp','vNew','vNewScaled');
% %plotchr = ['d';'*';'*';'*';'*'];
% plotchr = ['g';'r';'k';'c';'b'];
% aLen = length(a);
% v = zeros(aLen,nLen);
% figure;ski = 50;
% for p = 1:aLen
%     %v(p,:) = (1./n)*(1-a(p))^(-0.5);
%     v(p,:) = n*(1-a(p))^(0.5);
%     semilogy(n,v(p,:),plotchr(p));hold on;
%     %semilogy(n(1:ski:end),v(p,1:ski:end),plotchr(p),'linewidth',2);
% end
% grid on
% legend('\beta = 0.5','\beta = 0.6','\beta = 0.8','\beta = 0.9','\beta = 0.995');
% xlabel('Number of Trials','fontsize',12,'fontweight','bold')
% ylabel('Gittins'' Index','fontsize',12,'fontweight','bold')
% set(gca,'fontsize',12,'fontweight','bold')