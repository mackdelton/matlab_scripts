%%
%Load the desired dataset
%Format --> bernoulliGittens_AA_BB_C.mat
%AA: Number of arms considered
%BB: Number of trials run
%C: Arm selection method: 0 - Varaiya, 1 - Baseline(random),
%2 - Semi-intelligent
%Run first to create file
%fid = fopen('statistics_012_rndgrid_exp_11_15.txt','w')
%fprintf(fid,'Iteration 1 \t Iteration 2 \t Iteration 3 \n');fclose all

clear all
clc
close all
%Determines whether or not the figures generated will be saved to the
%working directory
saveFlag = 0;

textS = 14; %Define plot text font size
N=20; %Number of agents
%iter = [50 100 150 200]; %Number of iterations ("t")
iter = [100 500 1000]; %Number of iterations ("t")
v = [0 1 2]; %Solution version: 0 - Varaiya, 1 - Baseline(random),
      %2 - Semi-intelligent
analysisOpt = 1;
%Collection of statistics
statVecMedian = zeros(length(v),length(iter)); %Use to collect median
statVecMean = zeros(length(v),length(iter)); %Use to collect mean

%Setup writing to file
%fid = fopen('statistics_012_rndgrid_exp_11_15.txt','a+');
fid = fopen('bogus.txt','a+');
for ii = v
    fprintf(fid,'%d\t',ii);
    for iterate_i = 1:3
         eval(['load ./rndgrid_gamma/tempdata11/bernoulliGittins_' num2str(N) '_' num2str(iter(iterate_i)) '_' num2str(ii) '.mat;']); 
         iter = [100 500 1000]; %Number of iterations ("t")
         switch(analysisOpt)
            case {0} %Produce graphics highlighting how the arm most often selected
                     %exists at a range corresponding to a value on the
                     %heuristic that is greater than the pre-defined
                     %threshold, gamma
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
                    eval(['print(matchFig,''matchFig_' num2str(N) '_' num2str(iter(iterate_i)) '_' num2str(ii) ''',''-dpng'');']);
                end
             case{1} %Produce graphics highlighting the arms selected and
                     %how often messages were successfully transmitted.
                    ratioTrial = [];
                    for n=1:N
                        %For a given location, collect output of all trials
                        outTruth = histA((aId==n),4);
                        
                        %Record [no. of successes; no. of failures]
                        ratioTrial = [ratioTrial;[sum(outTruth) (length(outTruth)-sum(outTruth))]];
                    end
                    %Record statistics
                    %Calculate the percentage of successes out of total
                    %trials
                    summaryVec = ratioTrial(:,1)./(ratioTrial(:,1)+ratioTrial(:,2));
                    %Remove cases where there is neither a success or
                    %failure, i.e., the location was visited.
                    summaryVec(isnan(summaryVec)) = 0;
                    %Rows represent method (e.g., MAB, Random, etc.),
                    %columns represent the number of iterations attempted
                    statVecMedian(ii+1,iterate_i) = median(summaryVec);
                    statVecMean(ii+1,iterate_i) = mean(summaryVec);
                    fprintf(fid,'%3.3f\t',statVecMean(ii+1,iterate_i));
                    %Output statistics/plots
                    figure;subplot(1,2,1);bar(ratioTrial,'stacked');xlabel('Candidate Location');
                    legend('Successes','Failures');
                    subplot(1,2,2);hist((ratioTrial(:,1)+ratioTrial(:,2))/iter(iterate_i),20);
                    xlabel('% of Trials');ylabel('No. of locations');
             otherwise
                disp{'Unknown selection'}
         end
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
%statVecMedian
statVecMean

statDiff = [(statVecMean(2,:)-statVecMean(1,:))./statVecMean(1,:);(statVecMean(3,:)-statVecMean(1,:))./statVecMean(1,:)]

%%
% Commands for printing output for all results for all conditions as .csv output for analysis
% 1) Load data
% 2) Customize file name
% 3) Run
% Example:
N=[20 100]; %Number of agents
%iter = [50 100 150 200]; %Number of iterations ("t")
iter = [100 500 1000]; %Number of iterations ("t")
v = [0 1 2]; %Solution version: 0 - Varaiya, 1 - Baseline(random),
      %2 - Semi-intelligent
bern = [11 12 13 14 15 22 23 24 25 33 34 35 44 45 55];
h = waitbar(0, 'loading and writing data, WAIT!');
for ii = v
    for iterate_i = 1:3
        for b = 1:length(bern) %For bernoulli condition of success
            for n = N
                load(strcat({'/home/lparker/matlab_scripts/mab_scripts/meshgrid_gamma/tempdata',...
                     num2str(bern(b)),...
                    '/bernoulliGittins_'
                    num2str(n), '_', num2str(iterate_i(iter)), '_', num2str(ii), '.mat'}));
                    ob = ones(length(histA),1); %Generic 1's vector
                    csvwrite(strcat({'testdata_',...
                        num2str(n), '_', num2str(iterate_i(iter)), '_',...
                        num2str(ii), '_', 'meshgamma.csv'}),...
                        [aId' histA n*ob, ii*ob, iterate_i(iter)*ob, b*ob]);
                    csvwrite(strcat({'testloc_',...
                        num2str(n), '_', num2str(iterate_i(iter)), '_',...
                        num2str(ii), '_', 'meshgamma.csv'}),[locsA locsB])
            end
        end
        clear all %Clear previous data
    end
    waitbar(ii/length(v),h);
end
close(h)

%%
% %% Print script for comparing methods (primarily used for data for SMC2016)
% % Used to compare total distance accumulated by each method
% close all;clc;
% N = 20; iter_i = 200;
% ii=0;
% eval(['load bernoulliGittins_' num2str(N) '_' num2str(iter_i) '_' num2str(ii) '_101816.mat;']);
% v0Dist = distMax
% 
% ii=1;
% eval(['load bernoulliGittins_' num2str(N) '_' num2str(iter_i) '_' num2str(ii) '_101816.mat;']);
% v1Dist = distMax
% 
% ii=2;
% eval(['load bernoulliGittins_' num2str(N) '_' num2str(iter_i) '_' num2str(ii) '_101816.mat;']);
% v2Dist = distMax
%  
% methName = {'Gittins Index','Uniformed Random','Educated Guess'};
% methCmp = figure;
% ax = gca;
% bar([0 1 2],[v0Dist v1Dist v2Dist]);
% %ax.XTickLabel = {methName};
% set(ax,'XTickLabel',methName);
% 
% eval(['print(methCmp,''methCmpFig_' num2str(N) '_' num2str(iter) ''',''-dpng'');']);
% 
% %% Print script used to validate the usefulness of location selection by an MAB-inspired approach