function [outcome] = envReward(sofr, kSuccess, numTrials)
%This function generates an output (success or failure) of communication
%based on the conditions of the input.
%Created by: Lonnie Parker
%Created on: 10/18/16
%Last modified on: N/A
%Modification History:
%

%INPUTS
% sofr: A heuristic value representing the number of successes as a
% function of the communication range between two agents.
% kSuccess: The minimum number of successes that will determine the
% probability of whether or not a message was successful between two agents
% numTrials: Specifies the number of trials used to determine the
% probability of at least kSuccess successes.
% 
%OUTPUTS
% outcome: Scalar value (0 - failure or 1 - success)

outcome = 0; %Initialize as 0.
% Determine the probability of success based on sofr and the goal of at
% least kSuccess successes.

%Calculate the cummulative sum of probabilities of exactly 1 through
%kSuccess successes.
cummulativeSum = 0;
for f=1:kSuccess
    cummulativeSum = cummulativeSum + (nchoosek(numTrials,(f-1))*(sofr^((f-1)))*(1-sofr)^(numTrials-(f-1)));
end

p = 1 - cummulativeSum;

%Generate outcome as a function of p based on a uniform randomly generated
%number.

%When p is high (i.e. close to "1"), then the randomly generated number is
%more likely to fall under p and is considered a success (i.e., a "1")
if rand <= p 
    outcome = 1;
end