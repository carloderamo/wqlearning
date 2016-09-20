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


final=[state,p];
trainingSet=final(100:len-len/3+100,:);
testSet=final(len-len/3+100:len,:);
dlmwrite('trainingSet_d_7.txt',trainingSet);
dlmwrite('testSet_d_7.txt',testSet);

rewards=zeros(len,7);

for ind=1:1:7

% iterate on dataset
prevAction=0;
prevPrice=p(1);
c=0.;
diff=0.;
rew=zeros(len,1);
states=zeros(len,8);
actions=zeros(len,1);
for i=1:1:len
    action=state(i,ind);
    %action=macorsi(i);
    currentPrice=p(i);
    diff=currentPrice-prevPrice;
    c=0;

    if(action~=0 && prevAction~=action)
		c=0.0002;
    end
    if(prevAction==0)
		rew(i)=0-c;
    end
    if(prevAction==1)
		rew(i)=diff-c;
    end
    if(prevAction==2)
        rew(i)=-diff-c;
    end    
     
    prevPrice=currentPrice;
    prevAction=action;
    actions(i)=action;

%     states(i,1:7)=state(i);
%     states(i,8)=actions;
    states=[state,actions];
    rewards(i,ind)=rew(i);
    
end
 

%sars=[states,actions,p,rew];
profit=sum(rew(200:len,:))
% 
% stateMatrix=load('stati.txt');
% 
% 
% 
% 

end

profits=zeros(1,7);
for col=1:1:7
    profits(1,col)=sum(rewards(200:len,col))
%[~,indx]=ismember(a,b,'rows')
end

dlmwrite('rewardsprova_d.txt',[profits;rewards]);

plot(rewards(:,3))

