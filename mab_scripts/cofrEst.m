%%
%John Baylog test
%Collecting samples at different locations. Use C(r) to dictate
%whether or not those samples yield a 1 or 0. Visit each location
%an equal number of times (consider variable number of visits in
%case {4}. Use the average location of agentB (mean(agentB)) as
%reference.
%Use for simple-case example to explore estimation of success of
%communication curve. Case 3
runIter = 20;
N = 10;
betaVal = 0;
v = 3;
iter = 200;

maxSpace = 20; %Max space of navigation area (nmi)
maxR = sqrt((2*maxSpace)^2+maxSpace^2); %Defines the max distance between two agents

locsA = maxSpace*rand(N,2); %Define N random locations, set at 20nmi x 20nmi max
locsA = [(locsA(:,1) - maxSpace) (locsA(:,2) - maxSpace/2)];

locsB = maxSpace*rand(N,2); %Define N random locations, set at 20nmi x 20nmi max
locsB = [locsB(:,1) (locsB(:,2) - maxSpace/2)];

armNumA = length(locsA); %Total number of arms Agent A will consider
armNumB = length(locsB); %Total number of arms Agent B will consider

load socs.mat %Load SoC options
%Calculate distances between potential AgentA locations and expected
%location of AgentB
r = pdist([mean(locsB,1);locsA],'euclidean');
r = r(1:armNumA); %Remove excess values measuring inter-A arm distances
c = ss_est_gamma;
xx = linspace(0,maxR,length(c));
rA = interp1(xx,c,r,'linear');

histAVec = zeros(runIter,N);

initAgentId = ceil(armNumA*rand);
agentA = locsA(initAgentId,:); %Initialize current location of Agent A
numVisit = zeros(1,armNumA); %Tally of total number of visits to each location by agentA
numSuccess = zeros(1,armNumA); %Tally of total 
a = 0;c=0;
for u = 1:runIter
    for t = 1:iter
        a = a+1;
        if a == armNumA+1 %Cycle through all locations once reaching max arm number
            a = 1;
        end
        %Keep track of the number of visits at each location
        numVisit(a) = numVisit(a) + 1;
        %Use cofrEst to dictate the successful receipt at location a
        %Use binomial function to calculate the probability that you
        %receive at least 1 success given the weighted probability of
        %success, cofrEst, and the number of times visited
        %numSuccess(a) = numSuccess + (1-binopdf(0,numVisit(a),rA(a)));
        %LTP, 03/27/16: Attempt at implementing Bernoulli-like
        %reward system
        %if ((rand(1) < rA(a)) && (rand(1) > gThresh)) %Use threshold
        succVal = rand(1); %Draw a value at random
        if (succVal < rA(a)) %Don't use threshold
            numSuccess(a) = numSuccess(a) + 1;
        end
    end
    histAVec(u,:) = numSuccess; %Only retain the statistics of the locations
    numSuccess = zeros(1,armNumA); %Tally of total %Reset
end
yA = mean(histAVec,1);
eA = std(histAVec,1,1);
%Aggregate statistics with CofR
%[mean ; std ; C(r); r]
compA = [yA;eA;rA;r];

compA = sortrows(compA',4)';
errorbar(compA(1,:)/runIter,compA(2,:),'rx');hold on
plot([1:N],compA(3,:),'k')