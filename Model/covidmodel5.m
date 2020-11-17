% covid model ver.5
%%% This is the code for the simulations accompanying this arXiv paper:
%%% https://arxiv.org/abs/2004.01453
%%% The notation here is consistent with the notation in the paper.

%%% this code produces:
% Main Figure: Fig4/ FigS2 (SF).
% PWQ: Fig6/ FigS3 (SF).
% panels (o) and (n) in Fig3.

% to produce the figures activte (i.e. comment/uncomment) the relevant code.


%% clear & clc
clear;
clc;

%% parameters
T=150; %time of simulation.
p_NS=0.3; % probability of Asympomatic.
N_H=4000; % Number of households for alpha>0.
% N_H=10000; % Number of households for alpha=0.
% f1vec=linspace(0,0.5,15); %probability of defector.
f1vec=0; %probability of defector.
parameter=3; % the expacted value of outside degree.
PS=[0.78,0.14,0.08]; %[P_M,P_S,P_C].
PS=cumsum(PS);
p_MSC=[1/5,1/4,1/3]; %[r_M,r_S,r_C]
p_H=1/11; %r_H
p_V=1/13; %r_V
p_HR=0.85;
p_VR=0.5;
% HHSprob=[1,0,0,0,0,0]; %distrubation of household size for model without houses
HHSprob=[0.3,0.23,0.23,0.1,0.1,0.04]; %distrubation of household size for model with houses
relalizations=15;
randWeibull=@(lambda,k,m,n) lambda*((-log(1-rand(m,n))).^(1/k)); % generate number from Weibull distrubtion
randExp=@(lambda,m,n) (-1/lambda).*log(1-rand(m,n)); % generate number from Exponetial distrubtion
dt=1/96; % segment of time (in units of days)
steps=ceil(T/dt); %
HQparameter=0.5; % parameter of population wide quarantine.
f2=0;
f3=0;

%% integreted matrix 
% A - adajecny matrix of integrated network. 
% A_H - a.m. of the cliques of households.
% A_O - a.m. of the outside  links.
% household - cell of the members of households.
% node2HH - node -> household(node).
[A,A_H1,A_O,household,node2HH]=BA_generator_fast(N_H,parameter,HHSprob,0); % 0-ER 1-SF
sizeH=zeros(1,N_H); %vector of household sizes.
for k=1:N_H
    sizeH(k)= length(household{k});
end

%% relalization
% 0-UM, 1-HQ, 2-IQ, 3-AQ, 4-FQ.
% T_out/T_in = % The average time of outside/inside meeting in one day.
% p_infect = Probability of integrated link in house.
% HQparameter = the fraction of isolated households in half-quarantine (default= 0.5, change for PWQ)


[N,~]=size(A); %number of nodes
T_outvec=logspace(log10(0.005),log10(0.05),15); % alpha>0
% T_outvec=logspace(log10(0.005),log10(0.05),30); %SF alpha=0;
% T_outvec=logspace(log10(0.005),log10(0.05),30); %ER alpha=0;

p_infectHvec=linspace(0,1,10); %p_infect where A_O is Erdos Renyi.
% p_infectHvec=linspace(0,0.5,15); %p_infect where A_O is Scale Free.

T_outvecC=T_outvec([14,13,13,12,12,11]); % ER alpha>0
%T_outvecC=T_outvec([22,21,19]); % ER alpha=0
%T_outvecC=T_outvec([9,9,8,8,7,7]); % SF alpha>0
%T_outvecC=T_outvec([18,15,14]); % SF alpha=0


T_invecC=ones(1,length(T_outvecC));

% p_infectHvecC=ones(1,length(T_invecC)); %alpha=0
p_infectHvecC=p_infectHvec([2,4,2,4,2,3]); % alpha>0 ER
% p_infectHvecC=p_infectHvec([4,7,4,6,4,5]);% alpha>0 SF

HQparametervec=0.5; % Main Figure: 0.5 PWQ: [0.5,0.6,0.7,0.75]
iNodevec=randi(N,relalizations,1); % the index case node

for kTC=1:length(T_invecC)
    T_out=T_outvecC(kTC);
    T_in=T_invecC(kTC);
    p_infectH=p_infectHvecC(kTC);
    C=temporalmatrix(A_O,dt,T,(1/(2*T_out))); 
    p_in=T_in/(4*12); %inside infection probability.
    p_out=T_out/(12); %outside infection probability.
    
    householdNet=cell(1,N_H); % the adajecncy matrix of infectable links in any household.
    IN_MAT=zeros(N,N); % the adajeccy matrix of indoor infectable links.
    for k=1:N_H
        sizeH(k)= length(household{k});
        IN_H=rand(sizeH(k))<p_infectH;
        IN_H=triu(IN_H,1)+triu(IN_H,1)';
        householdNet{k}=sparse(IN_H);
        IN_MAT(household{k},household{k})=householdNet{k};
    end
    A_H = A_H1 & IN_MAT;

    for kHQparameter=1:length(HQparametervec)
        HQparameter=HQparametervec(kHQparameter); 
        HQmat=zeros(relalizations,ceil(HQparameter*N_H)); % the household are isolated
        HQmatNode=zeros(relalizations,N); % the node are isolated

        for k=(1:relalizations)
            HQmat(k,:)=(randsample(N_H,ceil(HQparameter*N_H)))';
        end

        for k=(1:relalizations)
            HQ=HQmat(k,:);
            for k1=1:length(HQ)
            HQmatNode(k,household{HQ(k1)})=ones(length(household{HQ(k1)}),1);
            end
        end

    for quarentineC=0:4 %Main Figure: 0:4, Fig3: 0, PWQ:1
        y=zeros(steps,12); %the number of nodes in any state.
        yP=y; % the last state of any node.
        logy=zeros(1,steps);
        I_in=0; % total inside infections. 
        I_out=0; % total outside infections.
        betaR=0; % the average beta on relaizations.
        Rsum=0;
        PeakR=0; % The peak of I(t) in relazions
        logPeak=0; % logPeak...
        relVecSS=zeros(steps,relalizations); % summtion on S(t) relaiztions
        relVecES=zeros(steps,relalizations);
        relVecIS=zeros(steps,relalizations);
        relVecHS=zeros(steps,relalizations);
        relVecVS=zeros(steps,relalizations);
        relVecRS=zeros(steps,relalizations);
        relVecDS=zeros(steps,relalizations);
        
        for f1idx=1:length(f1vec)
            f1=f1vec(f1idx);
        for rel=1:relalizations
            QQ=0; % QQ=1 there are qurantine
            Qt=0; % the time of quarantine start 
            Nodef1=rand(1,N)<f1;
            Nodef3=rand(1,N)<f3;

            states=zeros(N,12); % the state of any node
            pstates=zeros(N,12); % the last state of any node
            states(:,1)=ones(N,1);
            iNode=iNodevec(rel);
            Exposed=zeros(N,1); % exposed nodes
            states(iNode,2)= ceil((1/dt)*randWeibull(4.42,1.47,1,1)); %initial condtion
            states(iNode,1)=0;
            Exposed(iNode)=1;
            for step=1:steps
                pstates=states;
                week=floor(mod((step-Qt),4*24*7*2)/(4*24*7)); % week 0 or 1 (for IQ and AQ)

                if (sum(sum(~~pstates(:,6:8)))/N)<(((log(N))/N)) && (QQ==0)
                    quarentine=0;
                else
                    quarentine=quarentineC;
                    if QQ==0
                        Qt=step;    
                    end
                    QQ=1;
                end
                closeHomeI=unique(node2HH(pstates(:,6) | pstates(:,7) | pstates(:,8)));
                
                % Quarantine
                if quarentineC==1 || quarentineC==3 
                    HQ=HQmat(rel,:);
                    HQnode=HQmatNode(rel,:);
                end
                
                switch quarentine
                    case 0
                    closeHomeQ=[];
                    closeHomeQnode=zeros(1,N);
        
                    case 1
                    closeHomeQ=HQ;
                    closeHomeQnode=HQnode;
        
                    case 2
                    if week==0
                        closeHomeQ=1:N_H;
                        closeHomeQnode=ones(1,N);
                    end
        
                    if week==1
                        closeHomeQ=[]; 
                        closeHomeQnode=zeros(1,N);
                    end
                    
                    case 3
                    if week==0
                        closeHomeQ=HQ;
                        closeHomeQnode=HQnode;
                    end
        
                    if week==1
                        vec=ones(1,N_H);
                        vec(HQ)=0;
                        closeHomeQ=find(vec); 
                        closeHomeQnode=~HQnode;
                    end
                        
                    case 4
                        closeHomeQ=1:N_H;
                        closeHomeQnode=ones(1,N);
            
                end
                

      
                states=zeros(N,12);
                % link process 
                % S_O outside infected links, S_H the inside links.  
                if floor(mod(step,4*12*2)/(4*12))==0 && (quarentine~=4)
                    [S_O,S_H]=closeS(C{step},household,sizeH,closeHomeQ,closeHomeI,f1+f2,f3,Nodef1,Nodef3,closeHomeQnode,householdNet);  
                else
                    S_O=sparse(N,N);
                    S_H=A_H;
                end
                
                infects=pstates(:,4) | pstates(:,5) | (pstates(:,6) & rand(N,1)<f3);
                CurExposed1_O=((double(S_O)*double(infects)) & ~Exposed);
                CurExposed1_H=((double(S_H)*double(infects)) & ~Exposed);
                vecE_O=find(CurExposed1_O);
                vecE_H=find(CurExposed1_H);
                vecHELP=rand(1,length(vecE_H));

                vecE1_O=rand(1,length(vecE_O))<p_NS;
                vecE2_O=rand(1,length(vecE_O))<1; % !!!!
                vecE1_H=rand(1,length(vecE_H))<p_NS;
                vecE2_H=vecHELP<p_in; 
                states(vecE_O,1)= (~vecE2_O)'; %S
                states(vecE_H,1)= (~vecE2_H)'; %S
                states(vecE_O,2)= ceil((1/dt)*randWeibull(4.42,1.47,1,length(vecE_O))).*(vecE1_O & vecE2_O); %E_NS
                states(vecE_O,3)= ceil((1/dt)*randWeibull(2.21,1.47,1,length(vecE_O))).*(~vecE1_O & vecE2_O); %E_S
                states(vecE_H,2)= ceil((1/dt)*randWeibull(4.42,1.47,1,length(vecE_H))).*(vecE1_H & vecE2_H); %E_NS
                states(vecE_H,3)= ceil((1/dt)*randWeibull(2.21,1.47,1,length(vecE_H))).*(~vecE1_H & vecE2_H); %E_S
                Exposed= Exposed | states(:,2) | states(:,3);
                I_in=I_in+sum(vecE2_H);
                I_out=I_out+sum(vecE2_O);

                % node process
                states(:,1)=~Exposed; % S
                states(:,11)=pstates(:,11);
                states(:,12)=pstates(:,12);
                states(pstates(:,2)>0,2)=pstates((pstates(:,2)>0),2)-1; % E_NS
                states(pstates(:,2)==1,4)=ceil((1/dt)*randWeibull(7,1.47,1,sum(pstates(:,2)==1)));
                states(pstates(:,3)>0,3)=pstates((pstates(:,3)>0),3)-1; % E_S
                states(pstates(:,3)==1,5)=ceil((1/dt)*randWeibull(3.3,1.47,1,sum(pstates(:,3)==1)));    
                states(pstates(:,4)>0,4)=pstates((pstates(:,4)>0),4)-1; % I_NS
                states(pstates(:,4)==1,11)=ones(sum(pstates(:,4)==1),1);
                states(pstates(:,5)>0,5)=pstates((pstates(:,5)>0),5)-1; % I_PS
                vecPS=find(pstates(:,5)==1);
                for k=1:length(vecPS) 
                   states(vecPS(k),5)=0;
                   gT=find(rand<PS,1);
                   states(vecPS(k),5+gT)=ceil((1/dt)*randExp(p_MSC(gT),1,1));  
                end
                states(pstates(:,6)>0,6)=pstates((pstates(:,6)>0),6)-1; % I_M
                states(pstates(:,7)>0,7)=pstates((pstates(:,7)>0),7)-1; % I_S
                states(pstates(:,8)>0,8)=pstates((pstates(:,8)>0),8)-1; % I_C    
                states(pstates(:,6)==1,11)=ones(sum(pstates(:,6)==1),1); %I_M2R
                states(pstates(:,7)==1,9)=ceil((1/dt)*randExp(p_H,1,sum(pstates(:,7)==1))); %I_S2H
                states(pstates(:,8)==1,10)=ceil((1/dt)*randExp(p_V,1,sum(pstates(:,8)==1))); %I_C2V
    
                states(pstates(:,9)>0,9)=pstates((pstates(:,9)>0),9)-1; % H
                states(pstates(:,10)>0,10)=pstates((pstates(:,10)>0),10)-1; % V
                vecPH=rand(sum(pstates(:,9)==1),1)<p_HR;
                states(pstates(:,9)==1,11)=vecPH; %H2R
                states(pstates(:,9)==1,12)=~vecPH; %H2D    
                vecPV=rand(sum(pstates(:,10)==1),1)<p_VR;
                states(pstates(:,10)==1,11)=vecPV; %V2R
                states(pstates(:,10)==1,12)=~vecPV; %V2D    

                y(step,:)=y(step,:)+sum(~~states);
                
                disp([step/steps,rel,quarentineC,kTC]);
            end
            yRel=y-yP;
            yP=y;
            vecS=yRel(:,1)/N;
            vecE=sum(yRel(:,2:3)')/N;
            vecI=sum(yRel(:,6:8)')/N;
            vecH=yRel(:,9)/N;
            vecV=yRel(:,10)/N;
            vecR=yRel(:,11)/N;
            vecD=yRel(:,12)/N;
            logy=logy+log(vecI);
            [Ipeak,argI]=max(vecI);
            
            % beta calculation
            firs=find(vecI>0,1);
            fir=find(vecI>((log(N))/N),1);
            x2=firs+ceil((argI-firs)/2);
            x1=firs+ceil((argI-firs)/4);
            if x2~=x1
            pV=polyfit((x1:x2)*dt,log(vecI(x1:x2)),1);
            if ~isnan(pV(1)) && ~isinf(pV(1))
            betaR=betaR+pV(1);
            Rsum=Rsum+1;
            end
            end

            
            
            if (Ipeak>0)
            logPeak=logPeak+log(Ipeak);
            PeakR=PeakR+1;
            end
            
%             relVecSS((1:(steps-fir+1)),rel)=vecS(fir:steps); %Fig3
%             relVecES((1:(steps-fir+1)),rel)=vecE(fir:steps); %Fig3
            relVecIS((1:(steps-fir+1)),rel)=vecI(fir:steps);
%             relVecHS((1:(steps-fir+1)),rel)=vecH(fir:steps); %Fig3
%             relVecVS((1:(steps-fir+1)),rel)=vecV(fir:steps); %Fig3
%             relVecRS((1:(steps-fir+1)),rel)=vecR(fir:steps); %Fig3
%             relVecDS((1:(steps-fir+1)),rel)=vecD(fir:steps); %Fig3
            
            deltaD(rel)=yRel(end,12); % Main Figure
            Hpeak(rel)= max(yRel(:,9));
            Vpeak(rel)= max(yRel(:,10));
            
        end
        
%         figure(3)
%         hold on;
%         plot((1:steps)/96,mean(relVecIS(:,relVecIS(1,:)~=0),2),'LineWidth',4);
%         drawnow;

        end
        logy=logy/relalizations;
        Peak(kTC)=exp(logPeak/PeakR);
        betaR=betaR/Rsum;
        y=y/relalizations;

        betavecC(kTC)=betaR;
        I_in=I_in/relalizations;
        I_out=I_out/relalizations;
        fC(kTC)=I_in/(I_in+I_out);

        % plots
        figure(kTC); % Main Figure
        hold on;
        plot((1:steps)/96,mean(relVecIS(:,relVecIS(1,:)~=0),2),'LineWidth',4);
        drawnow;

%     figure(1) % panel (n)
%     hold on;
%     plot((1:steps)/96,mean(relVecSS(:,relVecIS(1,:)~=0),2),'LineWidth',4);
%     plot((1:steps)/96,mean(relVecES(:,relVecIS(1,:)~=0),2),'LineWidth',4);
%     plot((1:steps)/96,mean(relVecIS(:,relVecIS(1,:)~=0),2),'LineWidth',4);
%     plot((1:steps)/96,mean(relVecRS(:,relVecIS(1,:)~=0),2),'LineWidth',4);
%     
%     figure(2) %panel (o)
%     hold on;
%     plot((1:steps)/96,mean(relVecHS(:,relVecIS(1,:)~=0),2),'LineWidth',4);
%     plot((1:steps)/96,mean(relVecVS(:,relVecIS(1,:)~=0),2),'LineWidth',4);
%     plot((1:steps)/96,mean(relVecDS(:,relVecIS(1,:)~=0),2),'LineWidth',4);
    
        % Main Figure
        deltaDf(kTC,quarentineC)=mean(deltaD(relVecIS(1,:)~=0))/N;
        Hpeakf(kTC,quarentineC)=mean(Hpeak(relVecIS(1,:)~=0))/N;
        Vpeakf(kTC,quarentineC)=mean(Vpeak(relVecIS(1,:)~=0))/N;
    
    end
    end

 end



%% edit plots

for k1=1:length(T_invecC)
matCOLOR=[[0.77,0.35,0.07];[0.4,0,0];[0,0.6,0.6];[0,0,0.6];[0,0.2,0.2]];
axes2 = gca(figure(k1));
hold(axes2,'on');
box(axes2,'on');
set(axes2,'FontSize',60,'FontWeight','bold','LineWidth',8,...
    'PlotBoxAspectRatio',[7.5 9 1]);
plots=get(gca, 'Children');
% first figure
%   {'UM','HQ','IQ','AQ','FQ'},'EdgeColor',[1,1,1],'Location','northwest','FontSize',17);
set(gca,'TickLength',[0.02, 0.001]);
axes2.XAxis.MinorTick= 'on';
axes2.YAxis.MinorTick= 'on';
for k=1:5
plots(k).Color=matCOLOR(6-k,:);
plots(k).LineWidth=8;
end
end

%% edit plot-PWQ
% for k1=2:3 % 1 for (a) panel, 2:3 for panels (b)-(c).
% matCOLOR=[[0.4,0,0];[1,0.39,0.28];[0.98,0.5,0.45];[1,0.55,0];[1,0.84,0];[0,0,0.6]];
% axes2 = gca(figure(k1));
% hold(axes2,'on');
% box(axes2,'on');
% set(axes2,'FontSize',42,'FontWeight','bold','LineWidth',6,...
%     'PlotBoxAspectRatio',[3 4 1]);
% plots=get(gca, 'Children');
% % comment for panel (a)
% % legend(plots(6:1), {'50%','60%','70%','75%','80%','AQ'},'EdgeColor',[1,1,1],'Location','northwest','FontSize',17);
% set(gca,'TickLength',[0.04, 0.001]);
% set(gca,'XMinorTick','on');
% plots(1).Color=matCOLOR(6,:);
% plots(1).LineWidth=9;
% for k=2:6
% plots(k).Color=matCOLOR(7-k,:);
% plots(k).LineWidth=6;
% end
% set(gca,'YScale','log')
% end

%% edit panel (n)
% for k1=1
% matCOLOR=[[0,0.47,0.6];[1,0.61,0.07];[0.77,0.35,0.07];[0.66,0.82,0.56]];
% axes2 = gca(figure(k1));
% hold(axes2,'on');
% box(axes2,'on');
% set(axes2,'FontSize',60,'FontWeight','bold','LineWidth',8,...
%     'PlotBoxAspectRatio',[1 1 1]);
% plots=get(gca, 'Children');
% legend(plots(4:1), {'S(t)','E(t)','I(t)','R(t)'},'EdgeColor',[1,1,1],'Location','east','FontSize',60);
% set(gca,'TickLength',[0.02, 0.001]);
% for k=1:4
% plots(k).Color=matCOLOR(5-k,:);
% plots(k).LineWidth=8;
% end
% end

%% edit panel (o)

% for k1=2
% matCOLOR=[[0.4,0,0.2];[0.4,0.35,0.2];[0.65,0.65,0.65]];
% axes2 = gca(figure(k1));
% hold(axes2,'on');
% box(axes2,'on');
% set(axes2,'FontSize',60,'FontWeight','bold','LineWidth',8,...
%     'PlotBoxAspectRatio',[1 1 1]);
% plots=get(gca, 'Children');
%      legend(plots(3:1), {'H(t)','V(t)','D(t)'},'EdgeColor',[1,1,1],'Location','east','FontSize',60);
% set(gca,'TickLength',[0.02, 0.001]);
% for k=1:3
% plots(k).Color=matCOLOR(4-k,:);
% plots(k).LineWidth=8;
% end
% end


%%

for k1=(length(T_invecC)+1):(2*length(T_invecC)) %deltaD
   figure(k1)
   hold on;
   bar(deltaDf(k1-length(T_invecC),2:4)-deltaDf(k1-6,5));
end

for k1=(2*length(T_invecC)+1):(3*length(T_invecC)) %Hpeak
   figure(k1)
   hold on;
   bar(Hpeakf(k1-(2*length(T_invecC)),2:4));
end

for k1=(3*length(T_invecC)+1):(4*length(T_invecC)) %Vpeak
   figure(k1)
   hold on;
   bar(Vpeakf(k1-(3*length(T_invecC)),2:4));
end


%% BARS FIGURE
for k1=(length(T_invecC)+1):(4*length(T_invecC))
matCOLOR=[[0.4,0,0];[0,0.6,0.6];[0,0,0.6]];
% matCOLORb=[[0.8,0,0];[0,1,1];[0.6,0.6,1]];
axes2 = gca(figure(k1));
hold(axes2,'on');
box(axes2,'on');
set(axes2,'FontSize',72,'FontWeight','bold','LineWidth',10,...
    'PlotBoxAspectRatio',[4 3 1]);
plots=get(gca, 'Children');
plots(1).FaceColor='flat';
plots(1).CData=matCOLOR;
plots(1).LineWidth=2;

end