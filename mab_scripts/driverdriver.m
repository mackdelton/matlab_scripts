clear all
clc
close all
%Quick script to driver long-term data collection
maxNumTrials = 5;
for j = 1:maxNumTrials
%For each maxNumTrials, evaluate mabdriver for 1 to j total successes
    for i = 1:j
        %eval(['mkdir(''./meshgrid_exp'',''tempdata' num2str(i) '' num2str(j) ''');']);
        mabdriver(i,j)
    end
end

% mabdriver(1,4)
% mabdriver(1,5)
% mabdriver(3,3)
% mabdriver(2,4)
% mabdriver(2,5)
% mabdriver(3,4)
% mabdriver(3,5)
% mabdriver(4,4)
% mabdriver(4,5)
% mabdriver(5,5)
% 
