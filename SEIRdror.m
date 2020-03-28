%%% This is the code for the simulations accompanying this arXiv paper:
%%%
%%% The notation here is consistent with the notation in the paper.

%% model parameters
alpha = 1/11;
%beta = 0.45; %without social awareness
beta = 0.3; %with social awareness
gamma=1/5; 
vecf=linspace(0.01,0.99,99);
%vecf=0.3; % 0 means full compliance, 1 is the opposite
%Tshiftvec=linspace(1,30,100); % length of a shift
Tshiftvec=7; % length of a shift
N=7000000; % population size
t0=60; % time of intervention

%% real data
xData=linspace(9,33,13); % every two days
yData=[7,12,15,21,50,97,126,213,337,677,833,1442,2030]/N; % number of infected

%% 
% first SEIR
y0first=[N-1,1,0,0]; % N-1 susceptible, 1 exposed, 0 infected, 0 recovered

% running the model
ysol = ode45(@(t1,y1) SIRBeq(t1,y1,alpha,beta,gamma,0,sum(y0first(1:4)),0), [0,t0], y0first(1:4));
tB=linspace(0,t0,40)'; % time until intervention
yB=deval(ysol,tB)'; % model prediction for times tB
yT=yB;
tT=tB;
% figure(1);
% hold on;
% plot(tT(1:200),((yT(1:200,3)/N))); % we're interested in the short term (first 200 days)
% plot(xData,((yData)),'o');

for kshift=1:length(Tshiftvec)
Tshift=Tshiftvec(kshift);
fTs=floor(2000/Tshift); % number of shifts in 2000 days
Tsspan=linspace(t0,2000+t0,fTs+1);
for kf=1:length(vecf)
    f=vecf(kf);
    yT=yB;
    tT=tB;
    %yi=[yi_LC,yi_LD,yi_FC,yi_FD]  where yi_## :=[S_##,E_##,I_##,R_##]
    yi=[(1/2)*(1-f)*yB(end,:),(1/2)*(f)*yB(end,:),(1/2)*(1-f)*yB(end,:),(1/2)*(f)*yB(end,:)]; %initial condition separation of lock down
for k=1:fTs
%@% consider changing the names to, say, yInSol and yOutSol !!DONE!!
yOutSol = ode45(@(t1,y1) SEIRoutdoor(t1,y1,alpha,beta,gamma,sum(yi(5:16))), Tsspan(k:(k+1)), yi(5:16));
yInSol = ode45(@(t1,y1) SEIRindoor(t1,y1,alpha,gamma), Tsspan(k:(k+1)), yi(1:4));
t=linspace(Tsspan(k),Tsspan(k+1),40)';
yIn=deval(yInSol,t)';
yOut=deval(yOutSol,t)';
yi=[yOut(end,5:8),yOut(end,9:12),yIn(end,:),yOut(end,1:4)]; % Flip
y=yOut(:,1:4)+yOut(:,5:8)+yOut(:,9:12)+yIn;
yT=[yT;y];
tT=[tT;t];
% disp(k);
end
% figure(1);
% hold on;
% plot(tT,yT(:,1)/N);
% figure(2);
% hold on;
% plot(tT,yT(:,2)/N);
% figure(3);
% hold on;
% plot(tT,yT(:,3)/N);
% plot(tT(1:200),(yT(1:200,3)/N));
% plot(xData,(yData),'o');
% plot(xData,yData,'o');
% figure(4);
% hold on;
% plot(tT,yT(:,4)/N);
 [peak,argpeak]=max(yT(:,3)/N);
 Ipeak(kf)=peak;
 Tpeak(kf)=tT(argpeak);
 Rinf(kf)=yT(end,4)/N;
 disp(kf);
%  [peak,argpeak]=max(yT(:,3)/N);
%  Ipeak(kf,kshift)=peak;
%  Tpeak(kf,kshift)=tT(argpeak);
%  Rinf(kf,kshift)=yT(end,4)/N;
%  disp([kf,kshift])
end
 figure(1);
 hold on;
 plot(vecf,Rinf,'o');
 figure(2);
 hold on;
 plot(vecf,Ipeak,'o');
 figure(3);
 hold on;
 plot(vecf,Tpeak,'o');
%  [peak,argpeak]=max(yT(:,3)/N);
%  Ipeak(kshift)=peak;
%  Tpeak(kshift)=tT(argpeak);
%  Rinf(kshift)=yT(end,4)/N;
%  disp(kshift)
end
%  figure(1);
%  hold on;
%  plot(Tshiftvec,Rinf,'o');
%  figure(2);
%  hold on;
%  plot(Tshiftvec,Ipeak,'o');
%  figure(3);
%  hold on;
%  plot(Tshiftvec,Tpeak,'o');
%  figure(1);
%  hold on;
%  imagesc(vecf,Tshiftvec,Rinf);
%  figure(2);
%  hold on;
%  imagesc(vecf,Tshiftvec,Ipeak);
%   figure(3);
%  hold on;
%  imagesc(vecf,Tshiftvec,Tpeak);
%% SEIR without flip
% t0=400; % time of intervention
% y0first=[N-1,1,0,0]; % N-1 susceptible, 1 exposed, 0 infected, 0 recovered
% ysol = ode45(@(t1,y1) SIRBeq(t1,y1,alpha,beta,gamma,0,sum(y0first(1:4)),0), [0,t0], y0first(1:4));
% tB=linspace(0,t0,1000)'; % time until intervention
% yB=deval(ysol,tB)'; % model prediction for times tB
% yT=yB;
% tT=tB;
% figure(1);
% hold on;
% plot(tT,yT(:,1)/N);
% figure(2);
% hold on;
% plot(tT,yT(:,2)/N);
% figure(3);
% hold on;
% plot(tT,yT(:,3)/N);
% figure(4);
% hold on;
% plot(tT,yT(:,4)/N);