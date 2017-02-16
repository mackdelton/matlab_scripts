bern = {'1of1','1of5','1of10'};
socType = [0 1]; %Type of Success of Comms curve (0-gamma, 1-exp)
socTypeS = {'gamma','exp'};
NN=[20 100]; %Number of agents
distType = [0 1]; % Spatial distribution type (0-Random, 1-Mesh)
iters = [100 500 1000]; %Number of iterations ("t")
vs = [0 1 2]; %Solution version: 0 - Varaiya, 1 - Baseline(random),
      %2 - Semi-intelligent
h = waitbar(0, 'loading and writing data, WAIT!');

data_name = strcat('testdata12FEB17.txt');

%fid_data = fopen(data_name,'a+');
%fid_loc = fopen(loc_name,'a+');

for bn = 1:length(bern) %For bernoulli condition of success
    for dT = distType % Spatial distribution of candidate locations
        for sT = socType %Type of SoC curve
            for nn = NN % Number of candidate locations
                for iterate_i = iters % Number of iterations
                    for ii = vs % Solution type
                        [ii iterate_i bn nn]
                        load(char(strcat('/home/lparker/matlab_scripts/mab_scripts/data/cond_',...
                             bern(bn),'/',...
                            'dataout_',...
                            num2str(dT),'_',...
                            num2str(sT), '_', num2str(nn), '_',...
                            num2str(iterate_i), '_', num2str(ii), '.mat')));
                        ob = ones(length(histA),1); %Generic 1's vector
                        dlmwrite(data_name,...
                            [aId' histA bn*ob, dT*ob, sT*ob, nn*ob, iterate_i*ob, ii*ob],'-append');
                        %dlmwrite(loc_name,[locsA locsB],'-append');
                    end
                end
            end
        end
        clearvars -except bn bern dT distType sT socType nn NN iterate_i iters ii vs data_name loc_name h %Clear previous data
    end
    waitbar(ii/length(vs),h);
end
loc_name_rnd20 = strcat('testloc_rnd20.txt');
loc_name_rnd100 = strcat('testloc_rnd100.txt');
loc_name_mesh20 = strcat('testloc_mesh20.txt');
loc_name_mesh100 = strcat('testloc_mesh100.txt');
load('/home/lparker/matlab_scripts/mab_scripts/data/cond_1of1/dataout_0_0_20_100_0.mat')
dlmwrite(loc_name_rnd20,[locsA locsB]);
load('/home/lparker/matlab_scripts/mab_scripts/data/cond_1of1/dataout_1_0_20_100_0.mat')
dlmwrite(loc_name_mesh20,[locsA locsB]);
load('/home/lparker/matlab_scripts/mab_scripts/data/cond_1of1/dataout_0_0_100_100_0.mat')
dlmwrite(loc_name_rnd100,[locsA locsB]);
load('/home/lparker/matlab_scripts/mab_scripts/data/cond_1of1/dataout_1_0_100_100_0.mat')
dlmwrite(loc_name_mesh100,[locsA locsB]);


close(h)