clear all
clc
close all   
%Load data
load('bernoulliGittins_20_150_2.mat')
%Generate visualization
% Prepare the new file.
eval(['vidObj = VideoWriter(''bernoulliGittens_' num2str(N) '_' num2str(iter_i) '_' num2str(ii) '.avi'');']);
vidObj.FrameRate = 1;
open(vidObj);
%set(gca,'nextplot','replacechildren');
%h = waitbar(0,'Running CPA calculation...');
agentA = 'AUV';

width = 500;
height = 500;
hfig = figure;
set(hfig, 'Position', [100 100 width height])
% plot(locsA(:,1),locsA(:,2),'o','MarkerSize', 10,'LineWidth',2,...
%                  'MarkerFaceColor','b');hold on;
% plot(aB(1),aB(2),'d','MarkerSize', 10,'LineWidth',2,...
%                  'MarkerFaceColor','m');
for t=1:length(histA)
    plot(locsA(:,1),locsA(:,2),'o','MarkerSize', 10,'LineWidth',2,...
                'MarkerFaceColor','b');hold on;
    plot(histA(t,1), histA(t,2),'r*',...
                'MarkerSize', 10,'LineWidth',2,'MarkerFaceColor','g');hold on;
    text(histA(t,1)+1, histA(t,2),num2str(aId(t)),'FontSize',15,'FontWeight','Bold');
    plot(aB(1),aB(2),'d','MarkerSize', 10,'LineWidth',2,...
                'MarkerFaceColor','g');
    xlabel('Easting [m]','FontSize',12,'FontWeight','Bold');
    ylabel('Northing [m]','FontSize',12,'FontWeight','Bold');
    legend('Arm Locations','Agent A','Agent B')
    rngOfData = [locsA;aB];
    minX = min(rngOfData(:,1));
    maxX = max(rngOfData(:,1));
    minY = min(rngOfData(:,2));
    maxY = max(rngOfData(:,2));
    axis([minX maxX minY maxY]);    
    %set(gca,'nextplot','replacechildren'); %This line removes previous plotting history
    % Write each frame to the file.
    %currFrame = getframe(gcf,[74 47 800 600]);
    currFrame = getframe(hfig);%(fig);
    writeVideo(vidObj,currFrame);
    clf(hfig);
    %waitbar(t/length(tsPos_rx));
end
% Close the file.
close(vidObj);