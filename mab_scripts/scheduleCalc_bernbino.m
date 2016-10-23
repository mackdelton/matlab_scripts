function [bestArmHistory agentId agentB cofr] = scheduleCalc_bernbino(beta,armLocA,armLocB,opts,mR)
%Test script for schedule algorithm based on Weber Tutorial and L. Gregorio
%This script is modified to account/test for the impact of alternating
%arm/location selection to mimic TDMA communication protocol.
%Created by: Lonnie Parker
%Created on: 02/19/16
%Last modified on: 03/27/16
%Modification History: Created Bernoulli version of scheduleCalc.m
%(scheduleCalc_bern.m) to test unity reward and select arm based on
%proximity.
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
% bestArmHistory: Physical locations of the best arm selection at each time
% step, the separation distance at which that occurs.
% agentId: The id of the best arm selection at each time step.
% agentB: Send back the location of agentB (for now, since all
% agentA locations will be determined relative to a stationary AgentB).

rng(opts(2)); %Set seed value for reproduceability
load socs.mat %Load SoC options
%Calculate derivative of SoC
c = ss_est_gamma;
xx = linspace(0,mR,length(c));
dc = diff(c)./diff(xx);dc = dc/max(dc);
d2c = diff(dc)./diff(xx(1:end-1));
gThresh = 0.5*max(c); %Threshold for successful acomms, "Gamma"
gPlot = gThresh*ones(1,length(xx));

%Return C(r) data along with threshold used for decisions
cofr = [xx;c;gPlot];
%subplot(3,1,2);plot(xx(1:end-1),dc/max(dc),'r'); %First derivative (no
%longer relevant)
%subplot(3,1,3);plot(xx(1:end-2),d2c/(2*max(d2c)),'r') %Second derivative
%no longer relevant)

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
bestArmHistory = zeros(opts(2),3); %[X-pos Y-pos (Separation Dist)]
%allArmReward = allArmHistory;
%Calculate distances between potential AgentA locations and expected
%location of AgentB
r = pdist([mean(armLocB,1);armLocA],'euclidean');
r = r(1:armNumA); %Remove excess values measuring inter-A arm distances
rA = interp1(xx,c,r,'linear');

switch(opts(1))
    case {0} %Solution by Varaiya et al.
        h = waitbar(0,'Please wait...running through iterations');
        initAgentId = ceil(10*rand);
        agentA = armLocA(initAgentId,:); %Initialize current location of Agent A
        for t = 1:opts(2)
                %Determine AgentA's next location for reception
                for a = 1:armNumA
                    %dcdr = abs(interp1(xx(1:end-1),dc,r,'linear')); %Use this value to define the state transition matrix

                    %Define 2-state MDP for Varaiya solution
                    % s0-Poor/no communication (below threshold)
                    % s1-Good/yes communication (above threshold)
                    % NOTE: Threshold does not have to be defined.
                    %pA = [(1-dcdr) dcdr;dcdr (1-dcdr)];
                    %pA = stmCalc(agentA, armLocB, [xx;c], gThresh);
                    %rBias = (pdist([agentB;armLocA(a,:)],'euclidean')/mR);
                    
                    %LTP, 04/2/16: Rewards based on binomial distribution
                    %and C(r). Presume 5 transmissions at each arm location
                    %with a likelihood of C(r)
                    ll = 5; %Number of potential transmissions received at
                    %arm location
                    rew = [0;(1-binopdf(0,ll,rA(a)))];
                    
                    vA(a) = max(gittins_index_by_varaiya(beta,rew,pA));
                    %Scale by distance from current AgentA position
                    
                    %vA(a) = vA(a) * pdist([agentA;armLocA(a,:)],'euclidean');
                end
                % Look for closest of the best arms
                bestV = (vA == max(vA))
                %if ~isempty(min(vA(vA>0)))
                %    bestV = vA == (min(vA(vA>0))); %Identify closest location
                %    %Update current location of Agent A, otherwise, leave
                %    %as is.
                agentA = armLocA(bestV,:);
                %end
                agentId(t) = find(bestV==1);
                %Save [x-best y-best separation-best]
                bestArmHistory(t,:) = [agentA r(bestV)];
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
close(h)