% covid model ver.5
%%% This is the code for the simulations accompanying this arXiv paper:
%%% https://arxiv.org/abs/2004.01453
%%% The notation here is consistent with the notation in the paper.

%%% this code produces:
% Selective Isolation: FigS4

% For description of variables see covidmodel5.m

%% clear & clc
clear;
clc;

%% parameters
T=300; %time of simulation.
p_NS=[0.32;0.25]; % probability of Asympomatic.
N_H=4000; % Number of households for alpha>0.
% N_H=10000; % Number of households for alpha=0.
% f1vec=linspace(0,0.5,15); %probability of defector.
f1vec=0; %probability of defector.
parameter=7.5; % the expacted value of outside degree.
PS=[[0.82,0.12,0.06];[0.6,0.27,0.13]]; %[P_M,P_S,P_C].
PS=cumsum(PS,2);
p_MSC=[[1/5,1/4,1/3];[1/5,1/4,1/3]]; %[r_M,r_S,r_C]
p_H=[1/11;1/11]; %r_H
p_V=[1/13;1/13]; %r_V
p_HR=[0.86;0.79];
p_VR=[0.5;0.5];
% HHSprob=[1,0,0,0,0,0]; %distrubation of household size
HHSprob=[0.3,0.23,0.23,0.1,0.1,0.04]; %distrubation of household size
relalizations=15;
randWeibull=@(lambda,k,m,n) lambda*((-log(1-rand(m,n))).^(1/k));
randExp=@(lambda,m,n) (-1/lambda).*log(1-rand(m,n));
dt=1/96;
f2=0;
f3=0;
steps=ceil(T/dt);
HQparameter=0.5; % parameter of population wide quarantine.
p_old=0.2; % The fraction of old people in our population.

%% integreted matrix 
[A,A_H1,A_O,household,node2HH]=BA_generator_fast(N_H,parameter,HHSprob,0);
sizeH=zeros(1,N_H); %vector of household sizes.
for k=1:N_H
    sizeH(k)= length(household{k});
end

%% relalization
% 0-UM, 1-HQ, 2-IQ, 3-AQ, 4-FQ.
[N,~]=size(A); %number of nodes
T_outvec=logspace(log10(0.005),log10(0.05),15); 
p_infectHvec=linspace(0,1,10); %Probability of integrated link in house where A_O is Erdos Renyi.
T_outvecC=T_outvec(14);
T_invecC=ones(1,1);
p_infectHvecC=p_infectHvec(2);

HQparametervec=0.5;
HQmat=zeros(relalizations,ceil(N_H/2));
HQmatNode=zeros(relalizations,N);
iNodevec=randi(N,relalizations,1);
nodeSpi=rand(N,relalizations)<p_old;

for k=(1:relalizations)
    HQmat(k,:)=(randsample(N_H,ceil(N_H/2)))';
end

for k=(1:relalizations)
    HQ=HQmat(k,:);
    for k1=1:length(HQ)
    HQmatNode(k,household{HQ(k1)})=ones(length(household{HQ(k1)}),1);
    end
end


for kTC=1:length(T_invecC)
    T_out=T_outvecC(kTC);
    T_in=T_invecC(kTC);
    p_infectH=p_infectHvecC(kTC);
    C=temporalmatrix(A_O,dt,T,(1/(2*T_out))); 
    p_in=T_in/(4*12); %inside infection probability.
    p_out=T_out/(12); %outside infection probability.
    
    householdNet=cell(1,N_H);
    IN_MAT=zeros(N,N);
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
        HQmat=zeros(relalizations,ceil(HQparameter*N_H));
        HQmatNode=zeros(relalizations,N);

        for k=(1:relalizations)
            HQmat(k,:)=(randsample(N_H,ceil(HQparameter*N_H)))';
        end

        for k=(1:relalizations)
            HQ=HQmat(k,:);
            for k1=1:length(HQ)
            HQmatNode(k,household{HQ(k1)})=ones(length(household{HQ(k1)}),1);
            end
        end

    for quarentineC=0:3
        for helpM=0:1
        y=zeros(steps,12);
        yP=y;
        logy=zeros(1,steps);
        I_in=0;
        I_out=0;
        betaR=0;
        Rsum=0;
        PeakR=0;
        logPeak=0;
        relVecIS=zeros(steps,relalizations);
        relVecHS=zeros(steps,relalizations);
        relVecVS=zeros(steps,relalizations);
        relVecDS=zeros(steps,relalizations);
        
        for f1idx=1:length(f1vec)
            f1=f1vec(f1idx);
        for rel=1:relalizations
            NodeOld=nodeSpi(:,rel);
            QQ=0;
            Qt=0;
            Nodef1=rand(1,N)<f1;
            Nodef3=rand(1,N)<f3;

            states=zeros(N,12);
            pstates=zeros(N,12);
            states(:,1)=ones(N,1);
            iNode=iNodevec(rel);
            Exposed=zeros(N,1);
            states(iNode,2)= ceil((1/dt)*randWeibull(4.42,1.47,1,1)); %initial condtion 
            states(iNode,1)=0;
            Exposed(iNode)=1;
            for step=1:steps
                pstates=states;
                week=floor(mod((step-Qt),4*24*7*2)/(4*24*7));

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
                if floor(mod(step,4*12*2)/(4*12))==0 && (quarentine~=4)
                    [S_O,S_H]=closeS_spi(C{step},household,sizeH,closeHomeQ,closeHomeI,f1+f2,f3,Nodef1,Nodef3,closeHomeQnode,householdNet,NodeOld,QQ,helpM);  
                else
                    S_O=sparse(N,N);
                    S_H=A_H;
                end
                infects=pstates(:,4) | pstates(:,5) | (pstates(:,6) & rand(N,1)<f3);
                CurExposed1_O=((double(S_O)*double(infects)) & ~Exposed);
                CurExposed1_H=((double(S_H)*double(infects)) & ~Exposed);
                % Young
                vecE_O=find(CurExposed1_O & ~NodeOld);
                vecE_H=find(CurExposed1_H & ~NodeOld);
                vecE1_O=rand(1,length(vecE_O))<p_NS(1);
                vecE2_O=rand(1,length(vecE_O))<1; % !!!!
                vecE1_H=rand(1,length(vecE_H))<p_NS(1);
                vecE2_H=rand(1,length(vecE_H))<p_in; 
                states(vecE_O,1)= (~vecE2_O)'; %S
                states(vecE_H,1)= (~vecE2_H)'; %S
                states(vecE_O,2)= ceil((1/dt)*randWeibull(4.42,1.47,1,length(vecE_O))).*(vecE1_O & vecE2_O); %E_NS
                states(vecE_O,3)= ceil((1/dt)*randWeibull(2.21,1.47,1,length(vecE_O))).*(~vecE1_O & vecE2_O); %E_S
                states(vecE_H,2)= ceil((1/dt)*randWeibull(4.42,1.47,1,length(vecE_H))).*(vecE1_H & vecE2_H); %E_NS
                states(vecE_H,3)= ceil((1/dt)*randWeibull(2.21,1.47,1,length(vecE_H))).*(~vecE1_H & vecE2_H); %E_S
                Exposed= Exposed | states(:,2) | states(:,3);
                I_in=I_in+sum(vecE2_H);
                I_out=I_out+sum(vecE2_O);
                % Old
                vecE_O=find(CurExposed1_O & NodeOld);
                vecE_H=find(CurExposed1_H & NodeOld);
                vecE1_O=rand(1,length(vecE_O))<p_NS(2);
                vecE2_O=rand(1,length(vecE_O))<1; % !!!!
                vecE1_H=rand(1,length(vecE_H))<p_NS(2);
                vecE2_H=rand(1,length(vecE_H))<p_in; 
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
                states(pstates(:,2)==1 & NodeOld,4)=ceil((1/dt)*randWeibull(7,1.47,1,sum(pstates(:,2)==1 & NodeOld)));
                states(pstates(:,2)==1 & ~NodeOld,4)=ceil((1/dt)*randWeibull(7,1.47,1,sum(pstates(:,2)==1 & ~NodeOld)));                
                states(pstates(:,3)>0,3)=pstates((pstates(:,3)>0),3)-1; % E_S
                states(pstates(:,3)==1 & NodeOld,5)=ceil((1/dt)*randWeibull(3.3,1.47,1,sum(pstates(:,3)==1 & NodeOld)));
                states(pstates(:,3)==1 & ~NodeOld,5)=ceil((1/dt)*randWeibull(3.3,1.47,1,sum(pstates(:,3)==1 & ~NodeOld)));
                states(pstates(:,4)>0,4)=pstates((pstates(:,4)>0),4)-1; % I_NS
                states(pstates(:,4)==1,11)=ones(sum(pstates(:,4)==1),1);
                states(pstates(:,5)>0,5)=pstates((pstates(:,5)>0),5)-1; % I_PS
                % Young
                vecPS=find(pstates(:,5)==1 & ~NodeOld);
                for k=1:length(vecPS) 
                   states(vecPS(k),5)=0;
                   gT=find(rand<PS(1,:),1);
                   states(vecPS(k),5+gT)=ceil((1/dt)*randExp(p_MSC(1,gT),1,1));  
                end
                % Old
                vecPS=find(pstates(:,5)==1 & NodeOld);
                for k=1:length(vecPS) 
                   states(vecPS(k),5)=0;
                   gT=find(rand<PS(2,:),1);
                   states(vecPS(k),5+gT)=ceil((1/dt)*randExp(p_MSC(2,gT),1,1));  
                end
                
                states(pstates(:,6)>0,6)=pstates((pstates(:,6)>0),6)-1; % I_M
                states(pstates(:,7)>0,7)=pstates((pstates(:,7)>0),7)-1; % I_S
                states(pstates(:,8)>0,8)=pstates((pstates(:,8)>0),8)-1; % I_C    
                states(pstates(:,6)==1,11)=ones(sum(pstates(:,6)==1),1); %I_M2R
                states(pstates(:,7)==1 & NodeOld,9)=ceil((1/dt)*randExp(p_H(2),1,sum(pstates(:,7)==1 & NodeOld))); %I_S2H
                states(pstates(:,8)==1 & NodeOld,10)=ceil((1/dt)*randExp(p_V(2),1,sum(pstates(:,8)==1 & NodeOld))); %I_C2V
                states(pstates(:,7)==1 & ~NodeOld,9)=ceil((1/dt)*randExp(p_H(1),1,sum(pstates(:,7)==1 & ~NodeOld))); %I_S2H
                states(pstates(:,8)==1 & ~NodeOld,10)=ceil((1/dt)*randExp(p_V(1),1,sum(pstates(:,8)==1 & ~NodeOld))); %I_C2V
    
                states(pstates(:,9)>0,9)=pstates((pstates(:,9)>0),9)-1; % H
                states(pstates(:,10)>0,10)=pstates((pstates(:,10)>0),10)-1; % V
                % Young
                vecPH=rand(sum((pstates(:,9)==1) & ~NodeOld),1)<p_HR(1);
                states(pstates(:,9)==1 & ~NodeOld,11)=vecPH; %H2R
                states(pstates(:,9)==1 & ~NodeOld,12)=~vecPH; %H2D    
                vecPV=rand(sum((pstates(:,10)==1 & ~NodeOld)),1)<p_VR(1);
                states(pstates(:,10)==1 & ~NodeOld,11)=vecPV; %V2R
                states(pstates(:,10)==1 & ~NodeOld,12)=~vecPV; %V2D
                % Old
                vecPH=rand(sum((pstates(:,9)==1) & ~NodeOld),1)<p_HR(1);
                states(pstates(:,9)==1 & ~NodeOld,11)=vecPH; %H2R
                states(pstates(:,9)==1 & ~NodeOld,12)=~vecPH; %H2D    
                vecPV=rand(sum((pstates(:,10)==1 & ~NodeOld)),1)<p_VR(1);
                states(pstates(:,10)==1 & ~NodeOld,11)=vecPV; %V2R
                states(pstates(:,10)==1 & ~NodeOld,12)=~vecPV; %V2D


    
                y(step,:)=y(step,:)+sum(~~states);

                disp([step/steps,rel,helpM,kTC]);
            end
            yRel=y-yP;
            yP=y;
            vecI=sum(yRel(:,6:8)')/N;
            vecD=yRel(:,12)/N;
            vecH=yRel(:,9)/N;
            vecV=yRel(:,10)/N;
            logy=logy+log(vecI);
            [Ipeak,argI]=max(vecI);

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
            
            relVecIS((1:(steps-fir+1)),rel)=vecI(fir:steps);
            relVecHS((1:(steps-fir+1)),rel)=vecH(fir:steps);
            relVecVS((1:(steps-fir+1)),rel)=vecV(fir:steps);
            relVecDS((1:(steps-fir+1)),rel)=vecD(fir:steps);
            
            if quarentineC~=0
            deltaD(rel)=yRel(end,12);
            Hpeak(rel)= max(yRel(:,9));
            Vpeak(rel)= max(yRel(:,10));
            end
        end
        
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
        figure(1+(3*helpM));
        hold on;
        plot((1:steps)/96,mean(relVecHS(:,relVecIS(1,:)~=0),2),'LineWidth',4);
        drawnow;
        figure(2+(3*helpM));
        hold on;
        plot((1:steps)/96,mean(relVecVS(:,relVecIS(1,:)~=0),2),'LineWidth',4);
        drawnow;
        figure(3+(3*helpM));
        hold on;
        plot((1:steps)/96,mean(relVecDS(:,relVecIS(1,:)~=0),2),'LineWidth',4);
        drawnow;

        if quarentineC~=0
        deltaDf(quarentineC,helpM+1)=mean(deltaD)/N;
        Hpeakf(quarentineC,helpM+1)=mean(Hpeak)/N;
        Vpeakf(quarentineC,helpM+1)=mean(Vpeak)/N;
        end
        
        end
    end
    figure(7)
    hold on
    bar(Hpeakf)
    drawnow;

    figure(8)
    hold on
    bar(Vpeakf)
    drawnow;

    
    figure(9)
    hold on
    bar(deltaDf-(2.56*(10^-4)))
    drawnow;

    end
    

 end



%% edit plots

for k1=1:6
    if (k1>4)
        matCOLOR=[[0.77,0.35,0.07];[0.4,0,0];[0,0.6,0.6];[0,0,0.6]];
    else
        matCOLOR=[[1,0.4,0.4];[0.8,0,0];[0,1,1];[0.6,0.6,1]];
    end
axes2 = gca(figure(k1));
hold(axes2,'on');
box(axes2,'on');
set(axes2,'FontSize',40,'FontWeight','bold','LineWidth',8,...
    'PlotBoxAspectRatio',[4 3 1]);
plots=get(gca, 'Children');
%     legend(plots(4:1), {'UM','HQ','IQ','AQ'},'EdgeColor',[1,1,1],'Location','northeast','FontSize',32);
set(gca,'TickLength',[0.02, 0.001]);
for k=1:4
plots(k).Color=matCOLOR(5-k,:);
plots(k).LineWidth=6;
end
end

%% edit bars
for k1=7:9
matCOLOR=[[0.4,0,0];[0,0.6,0.6];[0,0,0.6]];
matCOLORb=[[0.8,0,0];[0,1,1];[0.6,0.6,1]];
axes2 = gca(figure(k1));
hold(axes2,'on');
box(axes2,'on');
set(axes2,'FontSize',40,'FontWeight','bold','LineWidth',8,...
    'PlotBoxAspectRatio',[4 3 1]);
plots=get(gca, 'Children');
plots(1).FaceColor='flat';
plots(1).CData=matCOLOR;
plots(2).FaceColor='flat';
plots(2).CData=matCOLORb;

% xtips1 = plots(1).XEndPoints;
% ytips1 = plots(1).YEndPoints;
% labels1 = string(plots(1).YData);
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
% 
% xtips1 = plots(2).XEndPoints;
% ytips1 = plots(2).YEndPoints;
% labels1 = string(plots(2).YData);
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')

end