clear all;

x=load('dailyGBP.csv');
%p=(x(:,2)+x(:,3))/2;

p=x(:,6);
len=length(p);


% RSI SIGNAL
rsi=rsindex(p,20);
Srsi1=(rsi<30);
Srsi2=(rsi>70)*2;
Srsi=Srsi1+Srsi2;

%Momentum
mom=tsmom(p);
Smom=(mom<0)+1;


%MACD    settare parametri macd
[macdvec,nineperma]=macd(p);


%mac=macd(p);
%ema9=movavg(p,9,9,'e');
Smacd=(macdvec<nineperma)+1;

%MA_CrossOver
[short,long]=movavg(p,20,200,0);
sMaCO=(short<long)+1;

%PCB
upper=hhigh(p,15);
lower=llow(p,15);
Spcb1=(p>=upper);
Spcb2=(p<=lower)*2;
Spcb=Spcb1+Spcb2;

%CCI
cci=cci(p,15);
Scci2=(cci>100)*2;
Scci1=(cci<-100);
Scci=Scci1+Scci2;

%STOC
h=x(:,4);
l=x(:,5);
stoc=stochosc(h,l,p,14);
sk=stoc(:,1);
sd=stoc(:,2);
Sstoch2=(sk<sd & sd>80 & sk>80)*2;
Sstoch1=(sk>sd & sd<20 & sk<20);

Sstoch=Sstoch1+Sstoch2;
SstochF=zeros(len,1);
min=llow(l,4);
max=hhigh(h,4);
lastSread=0;
for i=1:1:len
    if(Sstoch(i)~=0)
    lastSread=Sstoch(i);
    end
    if(p(i)<=min(i) || p(i)>=max(i))
    lastSread=0;
    end
    SstochF(i)=lastSread;
   
end

% prove di combinazioni di indicatori
stocrsi=(Srsi==SstochF & Srsi==Smacd).*Srsi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate final dataset
state=[Smacd,Smom,sMaCO,Srsi,Scci,Spcb,SstochF];
%%%% OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOo
%%%%%%%%%%%%%%%%%%%%%%%%

final = [state, p];
trainingSet = final(200:len - len / 3 + 200, :);
testSet = final(len - len / 3 + 200:len, :);

%%% Questo header indica quanti valori può avere il corrispondente indicatore
%%% Va cambiato ogni volta a seconda degli indicatori che vengono usati!!!!!!!
header = [2 2 2 3 3 3 3 3];

dlmwrite('trainingSet_d_7.txt', [header; trainingSet]);
dlmwrite('testSet_d_7.txt', [header; testSet]);

state = state(201:end, :);

finalState = zeros(9 * size(state, 1), size(state, 2) + 2);
finalState(:, 1:size(state, 2)) = repmat(state, 9, 1);

actionsCouples = [0 0; 1 1; 2 2; 0 1; 1 0; 0 2; 2 0; 1 2; 2 1]';
currentActions = 0;
a = 1;
episodeDim = size(finalState, 1) / 9;
for i = 0:episodeDim * 9 - 1
    if(mod(i, episodeDim) == 0)
        currentActions = currentActions + 1;
        a = 1;
    end
    finalState(i + 1, end - 1:end) = [actionsCouples(a, currentActions) actionsCouples(3 - a, currentActions)];
    a = 3 - a;
end

rew = zeros(size(finalState, 1), 1);
initialPrice = p(200);
prices = p(201:end);
prices = repmat(prices, 9, 1);

prevPrice = initialPrice;
for i = 1:size(finalState, 1)
    if(mod(i, size(state, 1) + 1) == 0)
        prevPrice = initialPrice;
    end
    
    currentPrice = prices(i);
    diff = currentPrice - prevPrice;
    cost = 0;
    if(finalState(i, end) ~= 0 && finalState(i, end - 1) ~= finalState(i, end))
        cost = 0.0002;
    end
    if(finalState(i, end - 1) == 0)
        rew(i) = 0 - cost;
    end
    if(finalState(i, end - 1) == 1)
        rew(i) = diff - cost;
    end
    if(finalState(i, end - 1) == 2)
        rew(i) = -diff - cost;
    end
    
    prevPrice = currentPrice;
end

sars = [finalState, rew];

states = zeros(size(sars, 1), 1);
for i = 1:size(sars, 1)
    stateN = 0;
    currentState = sars(i, 1:size(finalState, 2) - 1);
    
    for j = 1:length(currentState) - 1
        if(header(j) == 2)
            if(currentState(j) == 1)
                fact = 0;
            else
                fact = 1;
            end
        else
            fact = currentState(j);
        end
        stateN = stateN + fact * prod(header(j + 1:end));
    end
    if(header(end) == 2)
        if(currentState(end) == 1)
            fact = 0;
        else
            fact = 1;
        end
    else
        fact = currentState(end);
    end
    stateN = stateN + fact;

    states(i) = stateN;
end

actions = finalState(:, end);

% Header for Rele sars file
sarsRele = [1 1 1 nan nan];

stateActionRew = [zeros(size(states, 1), 2), states, actions, rew];
episodeDim = size(stateActionRew, 1) / 9;

for i = 0:8
    currentEpisode = stateActionRew(i * episodeDim + 1 : i * episodeDim + episodeDim, :);
    currentEpisode = [currentEpisode; [1 0 currentEpisode(end, 3) nan nan]];

    sarsRele = [sarsRele; currentEpisode];
end

csvwrite('sarsRele.csv', sarsRele);
