function [vAvec agentId] = scheduleCalc(beta,armLocA,armLocB,opts,mR)
%Test script for schedule algorithm based on Weber Tutorial and L. Gregorio
%This script is modified to account/test for the impact of alternating
%arm/location selection to mimic TDMA communication protocol.
%Created by: Lonnie Parker
%Created on: 02/19/16
%Last modified on: 03/16/16
%Naval Undersea Warfare Center DIVNPT

%INPUTS
% beta: A n x 1 vector of discount values to be tested for each solution
% armLocA: Specifies the location of each arm within Agent A's scope
% armLocB: Specifies the location of each arm within Agent B's scope
% opts: b x 1 vector containing the options with which a solution should
% be tested [Solution type, length of run], where Solution type can be:
% (0 - Varaiya, 1 - Sonin, 2 - Katehakis, 3 - Parker)
% mR: Maximum possible distance between agents
% OUTPUTS
% NONE
rng(1); %Set seed value for reproduceability
load socs.mat %Load SoC options
%Calculate derivative of SoC
c = ss_est_gamma;
xx = linspace(0,mR,length(c));
dc = diff(c)./diff(xx);dc = dc/max(dc);
d2c = diff(dc)./diff(xx(1:end-1));

figure;
subplot(3,1,1);plot(xx,c,'r')
subplot(3,1,2);plot(xx(1:end-1),dc/max(dc),'r')
subplot(3,1,3);plot(xx(1:end-2),d2c/(2*max(d2c)),'r')

dc = dc/max(dc); %Normalize

armNumA = length(armLocA); %Total number of arms Agent A will consider
armNumB = length(armLocB); %Total number of arms Agent B will consider
vA = zeros(1,armNumA); % Gittins indices for AgentA
vB = zeros(1,armNumB); % Gittins indices for AgentB

%History of Gittins indices
vAvec = zeros(opts(2),armNumA);
vBvec = zeros(opts(2),armNumB);
agentB = armLocB(1,:); %Freeze location of AgentB

%Define max reward as a function of max distance between AgentA arm
%locations
localMr = max(pdist(armLocA));
agentId = zeros(1,opts(2)); %Record agentA arm selection

%History of all arms selected
allArmHistory = zeros(1,armNumA);
allArmReward = allArmHistory;
%Calculate distances between potential AgentA locations and expected
%location of AgentB
r = pdist([mean(armLocB,1);armLocA],'euclidean');
r = r(1:armNumA); %Remove excess values measuring inter-A arm distances
%Define reward for state s1 as the inverse distance to potential location
%rA = rand; %Test reward
%rA = pdist([agentA+eps;armLocA(a,:)])/localMr;
rA = interp1(xx,c,r,'linear');

%Interested to see progression of rand impact
allArmRewardVec = zeros(opts(2),armNumA);
switch(opts(1))
    case {0} %Solution by Varaiya et al.
        h = waitbar(0,'Please wait...running through iterations');
        initAgentId = ceil(10*rand);
        agentA = armLocA(initAgentId,:); %Initialize current location of Agent A
%         allArmReward(initAgentId) = rA(initAgentId);
        allArmReward = rA;
        for t = 1:opts(2)
                %Determine AgentA's next location for reception
                for a = 1:armNumA
                    %dcdr = abs(interp1(xx(1:end-1),dc,r,'linear')); %Use this value to define the state transition matrix

                    %Define 2-state MDP for Varaiya solution
                    % s0-Poor/no communication (below threshold)
                    % s1-Good/yes communication (above threshold)
                    % NOTE: Threshold does not have to be defined.
                    %pA = [(1-dcdr) dcdr;dcdr (1-dcdr)];
                    pA = stmCalc(agentA, armLocB, [xx;c], 0.5*max(c));
                    %rBias = (pdist([agentB;armLocA(a,:)],'euclidean')/mR);
                    vA(a) = max(gittins_index_by_varaiya(beta,[0;allArmReward(a)],pA));
                end
                maxV = (vA==max(vA)); %Identify max Gittin's index
                agentA = armLocA(maxV,:); %Initialize current location of Agent A
                vAvec(t,:) = vA; %Store history
                allArmHistory = allArmHistory + maxV; %Keep running sum of which arm is selected
                agentId(t) = find(maxV == 1);
                %Calculate running sum of arm reward history based on
                %prior C(r) values
                %allArmReward(agentId(t)) = allArmReward(agentId(t)) + ...
                %                (rand^(t)*rA(agentId(t)))^(allArmHistory(agentId(t)));
                randBias = rand(1,armNumA).^(0.1*t);
                allArmReward = allArmReward + randBias.*rA;
                allArmRewardVec(t,:) = allArmReward;
                %Update waitbar
                waitbar(t/opts(2));
        end
        
%     case {1}
%         
%     case {2}
%         
%     case {3}
%         
%     case {4}
%         
    otherwise
        disp{'Unknown selection'}
end
figure;plot(allArmRewardVec(:,1),'r');hold on;
plot(allArmRewardVec(:,3),'b');plot(allArmRewardVec(:,9),'g');
plot(allArmRewardVec(:,10),'k');

pdist([agentB;agentA],'euclidean')
dList = pdist([agentB;armLocA],'euclidean');
dlist = dList(1:armNumA)
close(h)
        



