101816_data
This dataset contains initial evaluations of updates to the MAB script where a the environment yields an output. This allows us to determinine whether or not the agent's chosen location resulted in a success or failure.
****OPTIONS DEFINED IN mabdriver.m*****
N = 20
iter = [50 100 150 200]; %Number of iterations ("t")
v = [0 1 2]; %Solution version: 0 - Varaiya, 1 - Baseline(random),
      %2 - Semi-intelligent, 3 - Parker test
kk = 1; %Minimum number of successes required to be considered successful
        %comms
noN = 5; %Total number of "trials" or attempts at communicating between 
         %agents

101916_data
This dataset contains the same evaluations of updates found in 101816_data, but using 500, 1000, 2000, and 5000 iterations.
****OPTIONS DEFINED IN mabdriver.m*****
N = 20
iter = [500 1000 2000 5000]; %Number of iterations ("t")
v = [0 1 2]; %Solution version: 0 - Varaiya, 1 - Baseline(random),
      %2 - Semi-intelligent, 3 - Parker test
kk = 1; %Minimum number of successes required to be considered successful
        %comms
noN = 5; %Total number of "trials" or attempts at communicating between 
         %agents


101916b_data
This dataset contains the same evaluations of updates found in 101916_data, but the total number of trials (for the environmental response) is lowered to 1, i.e., noN=1 and only a single trial is allowed (instead of noN=5 for the previous two datasets).
****OPTIONS DEFINED IN mabdriver.m*****
N = 20
iter = [500 1000 2000 5000]; %Number of iterations ("t")
v = [0 1 2]; %Solution version: 0 - Varaiya, 1 - Baseline(random),
      %2 - Semi-intelligent, 3 - Parker test
kk = 1; %Minimum number of successes required to be considered successful
        %comms
noN = 1; %Total number of "trials" or attempts at communicating between 
         %agents

