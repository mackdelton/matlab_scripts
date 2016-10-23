function [vAvec, vBvec] = scheduleCalc_tdma(beta,armLocA,armLocB,opts,mR)
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

load socs.mat %Load SoC options
%Calculate derivative of SoC
xx = linspace(0,mR,length(ss_est_gamma));
dc = diff(ss_est_norm)./diff(xx);
d2c = diff(dc)./diff(xx(1:end-1));

figure;
subplot(3,1,1);plot(xx,ss_est_gamma,'r')
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
toggle = 1;
switch(opts(1))
    case {0} %Solution by Varaiya et al.
        agentA = armLocA(1,:); %Initialize current location of Agent A
        agentB = armLocB(1,:); %Initialize current location of Agent B
        for t = 1:opts(2)
            toggle = (-1)^t;
            if toggle < 0
                %Determine AgentA's next location for reception
                for a = 1:armNumA
                    %Calculate distance between potential AgentA locations
                    %and AgentB
                    r = pdist([agentB;armLocA(a,:)],'euclidean');
                    dcdr = abs(interp1(xx(1:end-1),dc,r,'linear')); %Use this value to define the state transition matrix

                    %Define 2-state MDP for Varaiya solution
                    % s0-Poor/no communication (below threshold)
                    % s1-Good/yes communication (above threshold)
                    % NOTE: Threshold does not have to be defined.
                    pA = [(1-dcdr) dcdr;dcdr (1-dcdr)];
                    
                    %Define reward for state s1 as the inverse distance to potential location
                    rA = 1/pdist([agentA+eps;armLocA(a,:)]);
                    vA(a) = max(gittins_index_by_varaiya(beta,[0;rA],pA));
                end
                maxV = (vA==max(vA)); %Identify max Gittin's index
                agentA = armLocA(maxV,:); %Initialize current location of Agent A
                vAvec(t,:) = vA; %Store history
            else %if (toggle > 0)
                %Determine AgentB's next location for reception
                for b = 1:armNumB
                    %Calculate distance between potential AgentB locations
                    %and AgentA
                    r = pdist([agentA;armLocB(b,:)],'euclidean');
                    dcdr = abs(interp1(xx(1:end-1),dc,r,'linear')); %Use this value to define the state transition matrix

                    %Define 2-state MDP for Varaiya solution
                    % s0-Poor/no communication (below threshold)
                    % s1-Good/yes communication (above threshold)
                    % NOTE: Threshold does not have to be defined.
                    pB = [(1-dcdr) dcdr;dcdr (1-dcdr)];
                    
                    %Define reward for state s1 as the inverse distance to potential location
                    rB = 1/pdist([agentB+eps;armLocB(b,:)]);
                    vB(b) = max(gittins_index_by_varaiya(beta,[0;rB],pB));
                end
                maxV = (vB==max(vB)); %Identify max Gittin's index
                agentB = armLocB(maxV,:); %Initialize current location of Agent B
                vBvec(t,:) = vB; %Store history
            end
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
        



