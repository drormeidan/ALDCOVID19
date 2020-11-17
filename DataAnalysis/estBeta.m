% estimate increase rates Fig3
clear;
close all;
load('timeseriescovid19deathsglobal.mat');
load('timeseriescovid19confirmedglobal.mat');

% parameters
countryvec=[138,226,202,137,121,176,83,7,170,10,17,224]; % the order of conutry include in the csv file
lockdownvec=[38,61,43,54,52,49,56,50,54,52,45,46];
Npopvec=[60,328,47,9,83,5.4,51,45,17,8,9,56]*(10^6);
effectPI=3;
effectNI=5;
effectPD=12;
effectND=8;

for k=1:length(countryvec)
country=countryvec(k);
lockdown=lockdownvec(k);
Npop=Npopvec(k);

% estimate beta from I(t)
vecxI=(lockdown-effectNI):(lockdown+effectPI);
vecyI=timeseriescovid19confirmedglobal(country,(lockdown-effectNI):(lockdown+effectPI))/(Npop);
pI=polyfit(vecxI(vecyI>(log(Npop)/Npop)),log(vecyI(vecyI>(log(Npop)/Npop))),1);
rateI(k)=pI(1);

% estimate beta from D(t)
vecxD=(lockdown-effectND):(lockdown+effectPD);
vecyD=timeseriescovid19deathsglobal(country,(lockdown-effectND):(lockdown+effectPD))/(Npop);
pD=polyfit(vecxD(vecyD>0),log(vecyD(vecyD>0)),1);
rateD(k)=pD(1);

firstConfirmed(k)=find(timeseriescovid19confirmedglobal(country,:)>0,1);

figure(k)
% subplot(1,2,1)
hold on;
plot(1:(82-firstConfirmed(k)+1),timeseriescovid19confirmedglobal(country,firstConfirmed(k):82)/Npop,'o');
plot(1:(82-firstConfirmed(k)+1),exp(polyval(pI,firstConfirmed(k):82)));
end
disp("beta: " + mean(rateI) + " +- " + std(rateI));


%% edit figures

for k1=1:length(countryvec)
matCOLOR=[[0.77,0.35,0.07];[0,0,1]];
axes2 = gca(figure(k1));
hold(axes2,'on');
box(axes2,'on');
set(axes2,'FontSize',72,'FontWeight','bold','LineWidth',8,...
    'PlotBoxAspectRatio',[1 1 1]);
plots=get(gca, 'Children');
set(gca,'TickLength',[0.02, 0.001]);
set(gca,'YScale','log');
plots(2).MarkerSize=14;
for k=1:2
plots(k).Color=matCOLOR(3-k,:);
plots(k).LineWidth=6;
end
end