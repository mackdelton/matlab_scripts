clear all; close all;clc

betaVal = 0;
%**********SETTINGS FOR DATA COLLECTION*****************
N = 10
iter = [50 100 150 200]; %Number of iterations ("t")
v = [0 1 2]; %Solution version: 0 - Varaiya, 1 - Baseline(random),
      %2 - Semi-intelligent, 3 - Parker test
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

        locsA = maxSpace*rand(N,2); %Define N random locations, set at 20nmi x 20nmi max
        locsA = [(locsA(:,1) - maxSpace) (locsA(:,2) - maxSpace/2)];

        locsB = maxSpace*rand(N,2); %Define N random locations, set at 20nmi x 20nmi max
        locsB = [locsB(:,1) (locsB(:,2) - maxSpace/2)];

        %figure;plot(locsA(:,1),locsA(:,2),'r*');hold on;plot(locsB(:,1),locsB(:,2),'b*')

        %Use for TDMA version
        %[histA histB] = scheduleCalc_tdma(betaVal,locsA,locsB,[v iter],maxR);

        %Use for single agent version
        %[histA aId] = scheduleCalc(betaVal,locsA,locsB,[v iter],maxR);

        %Use for single agent Bernoulli version
        [histA aId aB gRef distMax] = scheduleCalc_bern(betaVal,locsA,locsB,[ii iter_i],maxR);

        %Use for single agent Bernoulli version and binomial-defined rewards
        %[histA aId aB gRef] = scheduleCalc_bernbino(betaVal,locsA,locsB,[v iter],maxR);

        %Plot results for visual
        %eval(['save(''bernoulliGittens_' num2str(N) '_' ...
        %    num2str(iter) '.mat'',''histA'',''aId'',''aB'',''locsA'',''locsB'',''gRef'');']);
        eval(['save(''bernoulliGittins_' num2str(N) '_' num2str(iter_i) '_' num2str(ii) '_test.mat'');']);
        waitbar(ii/length(v));
    end
    close(hhh)
    waitbar(iter_i/length(iter));
    close(hh)
end
%**************************************************************************
%Use for simple-case example to explore estimation of success of
%communication curve
runIter = 20;
histAVec = zeros(runIter,N)
for u = 1:runIter
    [histA aId aB gRef distMax] = scheduleCalc_bern(betaVal,locsA,locsB,[v iter],maxR);
    histAVec(u,:) = histA(1,:);
end
yA = mean(histAVec,1);
eA = std(histAVec,1,1);
compA = [yA;eA;histA(1,:);histA(3:4,:)];

compA = sortrows(compA',5)';
errorbar(compA(1,:)/N,compA(2,:),'rx');hold on
plot([1:N],compA(4,:),'k')


%%
%Load the desired dataset
%Format --> bernoulliGittens_AA_BB_C.mat
%AA: Number of arms considered
%BB: Number of trials run
%C: Arm selection method: 0 - Varaiya, 1 - Baseline(random),
%2 - Semi-intelligent

% figure;
% subplot(1,2,1);stem(aId);
% xlabel('Iteration No.','FontSize',textS,'FontWeight','bold')
% ylabel('Arm Id','FontSize',textS,'FontWeight','bold')
% subplot(1,2,2);hist(aId)
% xlabel('Arm Id','FontSize',textS,'FontWeight','bold')
% ylabel('Arm Frequency','FontSize',textS,'FontWeight','bold')
saveFlag = 1;
textS = 14;
N=10;
iter = [50 100 150 200]; %Number of iterations ("t")
%iter = [50]; %Number of iterations ("t")

v = [0 1 2]; %Solution version: 0 - Varaiya, 1 - Baseline(random),
      %2 - Semi-intelligent, 3 - Parker test
%v = [0]; %Solution version: 0 - Varaiya, 1 - Baseline(random),

for iter_i = iter
    for ii = v
        
         eval(['load bernoulliGittins_' num2str(N) '_' num2str(iter_i) '_' num2str(ii) '.mat;']); 
%         %Redundant plots
%         stemFig = figure;
%         stem(aId);
%         xlabel('Iteration No.','FontSize',textS,'FontWeight','bold')
%         ylabel('Arm Id','FontSize',textS,'FontWeight','bold')
%         if(saveFlag == 1)
%             eval(['print(stemFig,''stemFig_' num2str(N) '_' num2str(iter_i) '_' num2str(ii) ''',''-dpng'');']);
%         end
% 
%         histFig = figure;
%         hist(aId)
%         xlabel('Arm Id','FontSize',textS,'FontWeight','bold')
%         ylabel('Arm Frequency','FontSize',textS,'FontWeight','bold')
%         if(saveFlag == 1)
%             eval(['print(histFig,''histFig_' num2str(N) '_' num2str(iter_i) '_' num2str(ii) ''',''-dpng'');']);
%         end

        %Plot resulting best arm selection on C(r), used primarily for
        %SMC2016 draft
        linew = 2;
        matchFig = figure;
        plot(gRef(1,:),gRef(2,:),'b','linewidth',linew);hold on;
        plot(gRef(1,:),gRef(3,:),'r--','linewidth',linew)
        mkRef = linspace(0,1);
        maxArm = mode(aId);
        distPt = pdist([locsA(maxArm,:);aB]);
        mkRefVal = distPt*ones(1,length(mkRef));hold on;
        plot(mkRefVal,mkRef,'k:','linewidth',linew);
        %Identify the location of the best arm on C(r)
        intersectPt = interp1(gRef(1,:),gRef(2,:),distPt,'linear');
        h = plot(distPt,intersectPt,'*');
        set(h, 'MarkerSize', 15, 'MarkerEdgeColor','k', 'LineWidth',3);
%         ptRef = [(distPt+5),(intersectPt+0.05)]
%         if (intersectPt < gRef(3,1))
%             annotation('textarrow',[ptRef(1) distPt+2]/max(gRef(1,:)),[ptRef(2) ptRef(2)],'String','Below')
%         else
%             annotation('textarrow',[ptRef(1) distPt+2]/max(gRef(1,:)),[ptRef(2) ptRef(2)],'String','Above')
%         end
        xlabel('Separation, r[nmi]','FontSize',textS,'FontWeight','bold');
        ylabel('Success of Communication, C(r)','FontSize',textS,'FontWeight','bold');
        FE = legend('C(r)','\gamma','r_{Best Arm}');
        LEG = findobj(FE,'type','text');
        set(LEG,'FontSize',14,'FontWeight','bold');

        if(saveFlag == 1)
            eval(['print(matchFig,''matchFig_' num2str(N) '_' num2str(iter_i) '_' num2str(ii) ''',''-dpng'');']);
        end
    end
end

%% Print script for comparing methods (primarily used for data for SMC2016)
close all;clc;
N = 20; iter_i = 150;
ii=0;
eval(['load bernoulliGittins_' num2str(N) '_' num2str(iter_i) '_' num2str(ii) '.mat;']);
v0Dist = distMax

ii=1;
eval(['load bernoulliGittins_' num2str(N) '_' num2str(iter_i) '_' num2str(ii) '.mat;']);
v1Dist = distMax

ii=2;
eval(['load bernoulliGittins_' num2str(N) '_' num2str(iter_i) '_' num2str(ii) '.mat;']);
v2Dist = distMax

methName = {'Gittens Index','Uniformed Random','Educated Guess'};
methCmp = figure;
ax = gca;
bar([0 1 2],[v0Dist v1Dist v2Dist]);
%ax.XTickLabel = {methName};
set(ax,'XTickLabel',methName);

eval(['print(methCmp,''methCmpFig_' num2str(N) '_' num2str(iter) ''',''-dpng'');']);