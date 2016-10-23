function [bestArmHistory, agentId, agentB, c_r, distTot] = scheduleCalc_bern(beta,armLocA,armLocB,opts,mR)
%Test script for schedule algorithm based on Weber Tutorial and L. Gregorio
%This script is modified to account/test for the impact of alternating
%arm/location selection to mimic TDMA communication protocol.
%Created by: Lonnie Parker
%Created on: 02/19/16
%Last modified on: 03/27/16
%Modification History:
%03/27/16: Created Bernoulli version of scheduleCalc.m
%(scheduleCalc_bern.m) to test unity reward and select arm based on
%proximity.
%10/10/16: Cleaned up for checking into git repo.
%Naval Undersea Warfare Center DIVNPT

%INPUTS
% beta: A n x 1 vector of discount values to be tested for each solution
% armLocA: Specifies the location of each arm within Agent A's scope
% armLocB: Specifies the location of each arm within Agent B's scope
% opts: 4 x 1 vector containing the options with which a solution should
% be tested [Solution type, length of run, minimum number of successes,
% total number of communication "trials" considered].
% Solution type can be:
% (0 - Varaiya, 1 - Baseline(random), 2 - Parker, 3 - Sonin, 3 - Katehakis)
% mR: Maximum possible distance between agents
% OUTPUTS
% bestArmHistory: Physical locations of the best arm selection at each time
% step, the separation distance at which that occurs.
% agentId: The id of the best arm selection at each time step.
% agentB: Send back the location of agentB (for now, since all
% agentA locations will be determined relative to a stationary AgentB).
% cofr: A 3 x n vector representing r, C(r), and the threshold used
% distTot: a 1x1 value representing the total distance traveled throughout
% all iterations.

%rng(opts(2)); %Set seed value for reproduceability
load socs.mat %Load SoC options
%Calculate derivative of SoC
c = ss_est_gamma;
xx = linspace(0,mR,length(c));
%dc = diff(c)./diff(xx);dc = dc/max(dc); %Unused (artifact from previous idea for implementation)
%d2c = diff(dc)./diff(xx(1:end-1)); %Unused (artifact from previous idea for implementation)
gThresh = 0.5*max(c); %Threshold for successful acomms, "Gamma"
gPlot = gThresh*ones(1,length(xx));

%Return C(r) data along with threshold used for decisions
c_r = [xx;c;gPlot];
%subplot(3,1,2);plot(xx(1:end-1),dc/max(dc),'r'); %First derivative (no
%longer relevant)
%subplot(3,1,3);plot(xx(1:end-2),d2c/(2*max(d2c)),'r') %Second derivative
%no longer relevant)

%dc = dc/max(dc); %Normalize %Unused (artifact from previous idea for implementation)

armNumA = length(armLocA); %Total number of arms Agent A will consider
armNumB = length(armLocB); %Total number of arms Agent B will consider
vA = zeros(1,armNumA); % Gittins indices for AgentA
vB = zeros(1,armNumB); % Gittins indices for AgentB %Unused since Agent B is considered stationary

%History of Gittins indices
%vAvec = zeros(opts(2),armNumA); %Unused
%vBvec = zeros(opts(2),armNumB); %Unused
agentB = armLocB(1,:); %Freeze location of AgentB
%Define max reward as a function of max distance between AgentA arm
%locations
%localMr = max(pdist(armLocA)); %Unused
agentId = zeros(1,opts(2)); %Record agentA arm selection

%History of all arms selected
bestArmHistory = zeros(opts(2),4); %[X-pos Y-pos (Separation Dist) success]

%Calculate distances between potential AgentA locations and expected
%location of AgentB
r = pdist([mean(armLocB,1);armLocA],'euclidean');
r = r(1:armNumA); %Remove excess values measuring inter-A arm distances
rA = interp1(xx,c,r,'linear');
distTot = 0; %Initialize total distance traversed between arm selections
switch(opts(1))
    case {0} %Solution by Varaiya et al.
        h = waitbar(0,'Please wait...running through intelligent iterations');
        initAgentId = ceil(armNumA*rand);
        agentA = armLocA(initAgentId,:); %Initialize current location of Agent A
        for t = 1:opts(2)
                %Determine AgentA's next location for reception
                for a = 1:armNumA
                    %dcdr = abs(interp1(xx(1:end-1),dc,r,'linear')); %Use this value to define the state transition matrix

                    %Define 2-state State Transition Matrix (STM) for Varaiya solution
                    % s0-Poor/no communication (below threshold)
                    % s1-Good/yes communication (above threshold)
                    % NOTE: Threshold does not have to be defined.
                    % (10/07/16, LTP: I have no idea why I said the
                    % threshold doesn't have to be defined)
                    %Unused, (artifact from previous idea for implementation)
                    %Thought that I would use the derivative of the C(r) to
                    %define the STM, but did not take that route.
                    %pA = [(1-dcdr) dcdr;dcdr (1-dcdr)];
                    pA = stmCalc(agentA, armLocB, [xx;c], gThresh);
                    %rBias = (pdist([agentB;armLocA(a,:)],'euclidean')/mR);
                    
                    %LTP, 03/27/16: Attempt at implementing Bernoulli-like
                    %reward system
                    %10/7/16, LTP: Potentially error prone way of
                    %determining rewards. May have just been a way to
                    %generate variety in definition of reward system.
                    if ((rand(1) < rA(a)) && (rand(1) > gThresh))
                        rew = [0;1];
                    else
                        rew = [0;0];
                    end
                    %Identify the location(s) with the max Gittins' Index
                    %(there may be more than one)
                    vA(a) = max(gittins_index_by_varaiya(beta,rew,pA));
                    
                    %Scale by distance from current AgentA position
                    %101816, LTP: 
                    vA(a) = vA(a) * pdist([agentA;armLocA(a,:)],'euclidean')/max(pdist(armLocA));
                end
                agentAold = agentA;
                % Look for closest point
                if ~isempty(min(vA(vA>0)))
                    bestV = vA == (min(vA(vA>0))); %Identify closest location
                    %Update current location of Agent A, otherwise, leave
                    %as is.
                    agentA = armLocA(bestV,:);
                end
                distTot = distTot + pdist([agentAold;agentA]);
                agentId(t) = find(bestV==1);
                                
                %Evaluate success or failure of decision
                mark = envReward(rA(agentId(t)),opts(3),opts(4));
                
                %Save [x-best y-best separation-best success/failure]
                bestArmHistory(t,:) = [agentA r(bestV) mark];
                
                %Update waitbar
                waitbar(t/opts(2));
        end
        
     case {1} %Uniformed Random (UR)
        h = waitbar(0,'Please wait...running through baseline iterations');
        initAgentId = ceil(armNumA*rand);
        agentA = armLocA(initAgentId,:); %Initialize current location of Agent A
        for t = 1:opts(2)
                %Determine AgentA's next location for reception
                agentAold = agentA;
                ranSel = ceil(armNumA*rand); %Choose random location
                agentA = armLocA(ranSel,:); %Move agentA to new arm location
                distTot = distTot + pdist([agentAold;agentA]);
                agentId(t) = ranSel;
                
                %Evaluate success or failure of decision
                mark = envReward(rA(agentId(t)),opts(3),opts(4));
                
                %Save [x-best y-best separation-best]
                bestArmHistory(t,:) = [agentA r(ranSel) mark];
                %Update waitbar
                waitbar(t/opts(2));
        end
        
     case {2} %Educated Guess (EG)
        h = waitbar(0,'Please wait...running through educated iterations');
        initAgentId = ceil(armNumA*rand);
        candidateNum = 100; %Candidate number of arms considered
        agentA = armLocA(initAgentId,:); %Initialize current location of Agent A
        for t = 1:opts(2)
                %Determine AgentA's next location for reception
                agentAold = agentA;
                %Determine AgentA's next location for reception
                for a = 1:armNumA
                    %Get min/max range of AgentB locations
                    minBx = min(armLocB(:,1));
                    maxBx = max(armLocB(:,1));
                    minBy = min(armLocB(:,2));
                    maxBy = max(armLocB(:,2));
                    
                    %Generate candidate locations of AgentB from which agentA
                    %will consider.
                    bCandidate = [(minBx + (abs(minBx-maxBx)*rand(1,candidateNum)));
                    (minBy + (abs(minBy-maxBy)*rand(1,candidateNum)))]';
                    %Calculate the proximity between AgentA and all
                    %candidate locations of AgentB
                    bCandidater = pdist([agentA;bCandidate]);
                    bCandidater = bCandidater(1:candidateNum); %Parse
                    
                    cofrCandidate = interp1(c_r(1,:),c_r(2,:),bCandidater,'linear');
                    %Identify the number of candidate locations greater
                    %than gamma 
                    vA(a) = length(cofrCandidate((cofrCandidate>gThresh)));
                end
                %Parse all arms under consideration for Agent A that
                %project maximum number of arms which AgentB could consider
                %that allow AgentA to remain within state S1
                vA = vA==max(vA);
                %Scale by distance between current AgentA location and
                %candidate arm location
                distA = pdist([agentAold;armLocA]);
                distA = distA(1:armNumA);
                
                %%%If one of the candidate arms for AgentA is the current
                %%%location, then remain
                %%%if sum(vA & (distA==0)) 
                %%%    agentId(t) = find(distA==0);
                    %Save [x-best y-best separation-best(last)]
                %%%    bestArmHistory(t,:) = [agentA r(distA==0)];
                %%%else
                    %Only re-assign AgentA a new location when a better
                    %nearest location exists.
                    vA = vA.*distA; %Scale
                    if ((~isempty(min(vA(vA>0)))) || ~(sum(vA)==0))
                        bestV = vA == (min(vA(vA>0))); %Identify closest location
                        %Update current location of Agent A, otherwise, leave
                        %as is.
                        agentA = armLocA(bestV,:);
                        agentId(t) = find(bestV==1);                        
                        %Save [x-best y-best separation-best mark(initialized with placeholder as 0)]
                        bestArmHistory(t,:) = [agentA r(bestV) 0];
                    else if(t==1)
                        agentId(t) = initAgentId;
                        else
                            agentId(t) = agentId(t-1);
                        end
                        %Save [x-best y-best separation-best mark(initialized with placeholder as 0)]
                        bestArmHistory(t,:) = [agentA r(agentId(t)) 0];
                    end
                    distTot = distTot + pdist([agentAold;agentA]);
                %%%end
                
                %Evaluate success or failure of decision
                mark = envReward(rA(agentId(t)),opts(3),opts(4));
                
                %Update environment output
                bestArmHistory(t,4) = mark;
                
                %Update waitbar
                waitbar(t/opts(2));
        end 
    otherwise
        disp{'Unknown selection'}
end
close(h)