clear all; close all;clc

betaVal = 0;
%**********OPTIONS FOR DATA COLLECTION*****************
N = 20
iter = [500 1000 2000 5000]; %Number of iterations ("t")
v = [0 1 2]; %Solution version: 0 - Varaiya, 1 - Baseline(random),
      %2 - Semi-intelligent, 3 - Parker test
kk = 1; %Minimum number of successes required to be considered successful
        %comms
noN = 1; %Total number of "trials" or attempts at communicating between 
         %agents
%**********SETTINGS FOR DATA COLLECTION*****************
hh = waitbar(0,'Please wait...running through iterations');
for iter_i = iter
    iter_i
    hhh = waitbar(0,'Please wait...running through types of methods');
    for ii = v
        ii
        rng(iter_i)
        maxSpace = 20; %Max space of navigation area (nmi)
        maxR = sqrt((2*maxSpace)^2+maxSpace^2); %Defines the max distance between two agents

        locsA = maxSpace*randn(N,2); %Define N random locations, set at 20nmi x 20nmi max
        locsA = [(locsA(:,1) - maxSpace) (locsA(:,2) - maxSpace/2)];
        %This is the new branch
        locsB = maxSpace*randn(N,2); %Define N random locations, set at 20nmi x 20nmi max
        locsB = [locsB(:,1) (locsB(:,2) - maxSpace/2)];

        %figure;plot(locsA(:,1),locsA(:,2),'r*');hold on;plot(locsB(:,1),locsB(:,2),'b*')

        %Use for TDMA version
        %[histA histB] = scheduleCalc_tdma(betaVal,locsA,locsB,[v iter],maxR);

        %Use for single agent version
        %[histA aId] = scheduleCalc(betaVal,locsA,locsB,[v iter],maxR);

        %Use for single agent Bernoulli version
        [histA aId aB gRef distMax] = scheduleCalc_bern(betaVal,locsA,locsB,[ii iter_i kk noN],maxR);

        %Use for single agent Bernoulli version and binomial-defined rewards
        %[histA aId aB gRef] = scheduleCalc_bernbino(betaVal,locsA,locsB,[v iter],maxR);

        %Plot results for visual
        %eval(['save(''bernoulliGittens_' num2str(N) '_' ...
        %    num2str(iter) '.mat'',''histA'',''aId'',''aB'',''locsA'',''locsB'',''gRef'');']);
        eval(['save(''./tempdata/bernoulliGittins_' num2str(N) '_' num2str(iter_i) '_' num2str(ii) '.mat'');']);
        waitbar(ii/length(v));
    end
    close(hhh)
    waitbar(iter_i/length(iter));
end
close(hh)
