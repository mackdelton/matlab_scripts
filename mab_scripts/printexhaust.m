bern = {'1of1','1of5','1of10'};
socType = [0 1]; %Type of Success of Comms curve (0-gamma, 1-exp)
socTypeS = {'gamma','exp'};
NN=[20 100]; %Number of agents
distType = [0 1]; % Spatial distribution type (0-Random, 1-Mesh)
iters = [100]; %Number of iterations ("t")
vs = [0 1]; %Solution version: 0 - Varaiya, 1 - Baseline(random),
      %2 - Semi-intelligent
h = waitbar(0, 'loading and writing data, WAIT!');

data_name = strcat('testdata03APR17_stationaryB_0.txt');
sB = 0;

%fid_data = fopen(data_name,'a+');
%fid_loc = fopen(loc_name,'a+');

for bn = 1:length(bern) %For bernoulli condition of success
    for dT = distType % Spatial distribution of candidate locations
        for sT = socType %Type of SoC curve
            for nn = NN % Number of candidate locations
                for iterate_i = iters % Number of iterations
                    for ii = vs % Solution type
                        [ii iterate_i bn nn]
                        load(char(strcat('/home/lparker/matlab_scripts/mab_scripts/data_stationaryB_0/cond_',...
                             bern(bn),'/',...
                            'dataout_',...
                            num2str(dT),'_',...
                            num2str(sT), '_', num2str(nn), '_',...
                            num2str(iterate_i), '_', num2str(ii), '.mat')));
                        ob = ones(length(histA),1); %Generic 1's vector
                        dlmwrite(data_name,...
                            [aId' histA bn*ob, dT*ob, sT*ob, nn*ob, iterate_i*ob, ii*ob distMax*ob],'-append');
                        %dlmwrite(loc_name,[locsA locsB],'-append');
                    end
                end
            end
        end
        clearvars -except bn bern dT distType sT socType nn NN iterate_i iters ii vs data_name loc_name h %Clear previous data
    end
    waitbar(ii/length(vs),h);
end
loc_name_rnd20 = strcat('testloc_rnd20_cmp.txt');
loc_name_rnd100 = strcat('testloc_rnd100_cmp.txt');
loc_name_mesh20 = strcat('testloc_mesh20_cmp.txt');
loc_name_mesh100 = strcat('testloc_mesh100_cmp.txt');
load('/home/lparker/matlab_scripts/mab_scripts/data_stationaryB_0/cond_1of1/dataout_0_0_20_100_0.mat')
dlmwrite(loc_name_rnd20,[locsA locsB]);
load('/home/lparker/matlab_scripts/mab_scripts/data_stationaryB_0/cond_1of1/dataout_1_0_20_100_0.mat')
dlmwrite(loc_name_mesh20,[locsA locsB]);
load('/home/lparker/matlab_scripts/mab_scripts/data_stationaryB_0/cond_1of1/dataout_0_0_100_100_0.mat')
dlmwrite(loc_name_rnd100,[locsA locsB]);
load('/home/lparker/matlab_scripts/mab_scripts/data_stationaryB_0/cond_1of1/dataout_1_0_100_100_0.mat')
dlmwrite(loc_name_mesh100,[locsA locsB]);


close(h)

%% Generate quick visual tool for analyzing evolution of Gittins Indices
% Load data of interest

%Select statistic to plot
% q = 1 -- Mean
% q = 2 -- StdDev
% q = 3 -- Number of arm pulls
% q = 4 -- Gittins Index

q = 3;

%Select the arms of interest (up to 4)

a1 = 2;
a2 = 3;
a3 = 4;
a4 = 5;

y1 = squeeze(gittinsHist(a1,a,:));
y2 = squeeze(gittinsHist(a2,a,:));
y3 = squeeze(gittinsHist(a3,a,:));
y4 = squeeze(gittinsHist(a4,a,:));
figure;
plot(y1,'r');hold on;
plot(y2,'g');
plot(y3,'c');
plot(y4,'c');