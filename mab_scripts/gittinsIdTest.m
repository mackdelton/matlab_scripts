%Identifying the trend used in the Gittins' Index table of values
%(reference [25] from M. Cheung's MIT MS thesis.

%n = [1:9,10:10:90,100:100:1000];
n = [1:500];
nLen = length(n);
%Use the following values for a = 0.5:0.1:0.9,0.95,0.99,0.995
%a = [0.5:0.1:0.9,0.95,0.99,0.995];
a = [0 0.5 0.8 0.9 0.995]; %trial of different beta values
plotchr = ['d';'<';'*';'^';'>'];
aLen = length(a);
v = zeros(aLen,nLen);
figure;ski = 50;
for p = 1:aLen
    v(p,:) = (1./n)*(1-a(p))^(-0.5);
    semilogy(n,v(p,:));hold on;
    semilogy(n(1:ski:end),v(p,1:ski:end),plotchr(p),'linewidth',2);
end
grid on
legend('\beta = 0','\beta = 0.5','\beta = 0.8','\beta = 0.9','\beta = 0.995');
xlabel('Number of Trials','fontsize',12,'fontweight','bold')
ylabel('Gittins'' Index','fontsize',12,'fontweight','bold')
set(gca,'fontsize',12,'fontweight','bold')