function [stmOut] = stmCalc(locStart,locEnd,cofr,gamma)
%This script is use to estimate a state transition matrix of arm locations
%defined in locX relative to locY.
%Created by: Lonnie Parker
%Created on: 03/20/16
%Last modified on: 03/20/16
%Naval Undersea Warfare Center DIVNPT

%INPUTS
% locStart: Receive arm locations (n x 2 vector) to test relative to changing state impacted
% by distance to locEnd.
% locEnd: Transmit arm locations (n x 2 vector)
% cofr: A n x 2 vector of values measuring the success of communication for
% a given range of separation values. First row are the SoC values, second
% row are the separation values.
% gamma: The threshold of communication dictating whether the system is in
% state s0 or s1
% OUTPUTS
% stmOut: A 2 x 2 state transition matrix

N = 100; %Size of sample space
stmOut = zeros(2,2); %Initialize state-transition matrix
a = [0 0]; %Counter of state transitions

% 10FEB17, lparker: Removed to eliminate confusion. State transition matrix
% will take all of SoC into account.
% 23FEB17, LTP: Additional modification will calculate STM as a function of
% two agents 
%Identify the range of separation distances for locStart
rRangeVec = [];
for ee = 1:length(locStart)
    rRange = pdist([locStart(ee,:);locEnd]);
    rRange = rRange(1:length(locEnd));
    rRangeVec = [rRangeVec rRange];
end

%rRange = pdist([locStart;locEnd]);
%rRange = rRange(1:length(locEnd));
subRange = [min(rRangeVec) max(rRangeVec)];
subIndCofr = (cofr(1,:) >= subRange(1)) & (cofr(1,:) <= subRange(2));
cofr = cofr(:,subIndCofr); %Redefine the range of C(r) values that should be considered for the given arm

for w = 1:2
    %For N iterations, select a deltaR and record whether or not the state has
    %changed.
    if w == 1
        %Randomly select an r value that determines a state of s0
        sInd = cofr(2,:) < gamma;
    else %Randomly select an r value that determines a state of s1
        sInd = cofr(2,:) >= gamma;
    end
    if (sum(sInd) ~= 0)
        stateCofr = cofr(:,sInd);
        rRef = datasample(stateCofr(1,:),1);
        for n = 1:N
            %Select delta r, according to uniform distribution and scaled to
            %the range of r values provided in cofr
            deltaR = (max(cofr(1,:))*rand)-rRef;
            %THIS IS A PLACEHOLDER FOR TESTING THE SELECTION OF DELTAR ACCORDING TO
            %A NORMAL DISTRIBUTION WHERE THE MEAN IS THE SEPARATION THAT ACHIEVES
            %THE MAXIMUM PROBABILITY.
            % deltaR = (max(rVec)*randn)-rRef;
            %
            %    
            %Define C(r+deltaR)
            newC = interp1(cofr(1,:),cofr(2,:),(rRef+deltaR),'linear');

            if w == 1
                %Evaluate new C(r) value for which state it is in
                if newC < gamma %NO state change from s0
                    a(w) = a(w) + 1;
                end %else, a change DID occur to s1
            else %w == 2
                %Evaluate new C(r) value for which state it is in
                if newC > gamma %NO state change from s1
                    a(w) = a(w) + 1;
                end %else, a change DID occur to s0
            end
        end
    else
        if w == 1    % If searching for candidate s0 states, but none exist
            a(2) = N;% then s1 is absolute
        else
            a(1) = N;% If searching for candidate s1 states, but none exist
        end          % then s0 is absolute
    end
end
stmOut = [a(1), (N-a(1)); ...
        (N-a(2)), a(2)]/N;        
    
