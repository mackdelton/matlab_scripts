function [g] = calculateGidx(be, infoSet, tStep)
%Test script for schedule algorithm based on Weber Tutorial
%Created by: Lonnie Parker
%Created on: 02/22/16
%Naval Undersea Warfare Center DIVNPT

%INPUTS
% be: Designated beta value to be used for calculating index option
% infoSet: Set of all locations with reward info used within scenario
% [xloc yloc reward]
% tStep: Time epoch of index evaluation
% OUTPUTS
% g: Greatest index.

len = length(infoSet); %Record size of index vector
gSet = zeros(1,len); %Define index vector for record keeping

for ii = 1:len
    gSet(ii) = [( infoSet(ii,3) * (be^(tStep)) * (1-be) ) ...
                / (1-(be^(tStep))) ii];
end

sortrows(gSet)
g = gSet(1,3); %Return the location with the highest index. 