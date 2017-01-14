NN=[20 100]; %Number of agents
%iter = [50 100 150 200]; %Number of iterations ("t")
iters = [100 500 1000]; %Number of iterations ("t")
vs = [0 1 2]; %Solution version: 0 - Varaiya, 1 - Baseline(random),
      %2 - Semi-intelligent
bern = [11 12 13 14 15 22 23 24 25 33 34 35 44 45 55];
h = waitbar(0, 'loading and writing data, WAIT!');

data_name = strcat('testdata_','rndexp.txt')
loc_name = strcat('testloc_','rndexp.txt')
%fid_data = fopen(data_name,'a+');
%fid_loc = fopen(loc_name,'a+');
for ii = vs
    for iterate_i = iters
        for bn = bern %For bernoulli condition of success
            for nn = NN
                [ii iterate_i bn nn]
                load(char(strcat('/home/lparker/matlab_scripts/mab_scripts/rndgrid_exp/tempdata',...
                     num2str(bn) ,...
                    '/bernoulliGittins_',...
                    num2str(nn), '_', num2str(iterate_i), '_', num2str(ii), '.mat')));
                    ob = ones(length(histA),1); %Generic 1's vector
                    dlmwrite(data_name,...
                        [aId' histA nn*ob, ii*ob, iterate_i*ob, bn*ob],'-append');
                    dlmwrite(loc_name,[locsA locsB],'-append');
            end
        end
        clearvars -except ii vs iterate_i iters bn bern nn NN data_name loc_name h %Clear previous data
    end
    waitbar(ii/length(vs),h);
end
close(h)