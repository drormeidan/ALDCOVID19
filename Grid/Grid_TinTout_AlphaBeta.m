% covid model ver.5
%%% This is the code for the simulations accompanying this arXiv paper:
%%% https://arxiv.org/abs/2004.01453
%%% The notation here is consistent with the notation in the paper.

%%% This code calculate grid of alpha and beta where given T_in and T_out
%%% values.

% For description of variables see covidmodel5.m

%% clear & clc
clear;
clc;

%% parameters
T=100; %time of simulation.
p_NS=0.3; % probability of Asympomatic
N_H=4000; % Number of households.
t0=10;
f1=0;
f2=0;
f3=0;
p_infectHvec=linspace(0,0.5,15);
% p_infectHvec=1; %alpha=0
parameter=7.5;
PS=[0.78,0.14,0.08];
PS=cumsum(PS);
p_MSC=[1/5,1/4,1/3];
p_H=1/11;
p_V=1/13;
p_HR=0.85;
p_VR=0.5;
HHSprob=[0.3,0.23,0.23,0.1,0.1,0.04];
%  HHSprob=[1,0,0,0,0,0];
relalizations=1;
randWeibull=@(lambda,k,m,n) lambda*((-log(1-rand(m,n))).^(1/k));
randExp=@(lambda,m,n) (-1/lambda).*log(1-rand(m,n));
dt=1/96;
steps=ceil(T/dt);

%% integreted matrix 
[A,A_H1,A_O,household,node2HH]=BA_generator_fast(N_H,parameter,HHSprob,1);
sizeH=zeros(1,N_H);
for k=1:N_H
    sizeH(k)= length(household{k});
end

%% relalization
% 0-UM, 1-HQ, 2-IQ, 3-AQ, 4-FQ.
[N,~]=size(A);
   T_invec=1;
 T_outvec=logspace(log10(0.005),log10(0.05),15);

HQmat=zeros(relalizations,ceil(N_H/2));
for k=(1:relalizations)
    HQmat(k,:)=(randsample(N_H,ceil(N_H/2)))';
end

householdNetCell=cell(1,length(p_infectHvec));
A_HCell=cell(1,length(p_infectHvec));
for kIn=1:length(p_infectHvec)
    p_infectH=p_infectHvec(kIn);
    householdNetCell{kIn}=cell(1,N_H);
    IN_MAT=zeros(N,N);
    for k5=1:N_H
        IN_H=rand(sizeH(k5))<p_infectH;
        IN_H=triu(IN_H,1)+triu(IN_H,1)';
        householdNetCell{kIn}{k5}=sparse(IN_H);
        IN_MAT(household{k5},household{k5})=householdNetCell{kIn}{k5};
    end
    A_HCell{kIn} = A_H1 & IN_MAT;
end

fCell=cell(1,5);
betavecCell=cell(1,5);
HTvecCell=cell(1,5);


parfor K=1:15
f=zeros(length(p_infectHvec),length(T_outvec));
betavec=zeros(length(p_infectHvec),length(T_outvec));
HTvec=zeros(length(p_infectHvec),length(T_outvec));

    for kIn=1:length(p_infectHvec)
        A_H=A_HCell{kIn};
        householdNet=householdNetCell{kIn};
    
for kT_in=1:length(T_invec)
for kT_out=1:length(T_outvec)
    
    T_out=T_outvec(kT_out);
    C=temporalmatrix(A_O,dt,T,(1/(2*T_out)));
    T_in=T_invec(kT_in);
    p_in=T_in/(4*12); %inside infection probability.
    p_out=T_out/(12); %outside infection probability.

    for quarentineC=0
        y=zeros(steps,12);
        yP=y;
        I_in=0;
        I_out=0;
        betaR=0;
        Rsum=0;

        for rel=1:relalizations
            states=zeros(N,12);
            pstates=zeros(N,12);
            states(:,1)=ones(N,1);
            iNode=randi(N);
            Exposed=zeros(N,1);
            states(iNode,2)= ceil((1/dt)*randWeibull(4.42,1.47,1,1)); %initial condtion 
            states(iNode,1)=0;
            Exposed(iNode)=1;
            for step=1:steps
                pstates=states;
                week=floor(mod((step-4*24*t0),4*24*7*2)/(4*24*7));
                if step<4*24*t0
                    quarentine=0;
                else
                    quarentine=quarentineC;
                end
                closeHomeI=unique(node2HH(pstates(:,6) | pstates(:,7) | pstates(:,8)));
                % Quarantine
                if quarentineC==1 || quarentineC==3 
                    HQ=HQmat(rel,:);
                end
                switch quarentine
                    case 0
                    closeHomeQ=[];
        
                    case 1
                    closeHomeQ=HQ;
        
                    case 2
                    if week==0
                        closeHomeQ=1:N_H;
                    end
        
                    if week==1
                        closeHomeQ=[]; 
                    end
                    
                    case 3
                    if week==0
                        closeHomeQ=HQ;
                    end
        
                    if week==1
                        vec=ones(1,N_H);
                        vec(HQ)=0;
                        closeHomeQ=find(vec); 
                    end
                        
                    case 4
                        closeHomeQ=1:N_H;
            
                end
      
                states=zeros(N,12);
                % link process
                if floor(mod(step,4*12*2)/(4*12))==0
                    [S_O,S_H]=closeS_grid(C{step},household,sizeH,closeHomeQ,closeHomeI,f1+f2,f3,householdNet);
                else
                    S_O=sparse(N,N);
                    S_H=A_H;
                end
                infects=pstates(:,4) | pstates(:,5) | (pstates(:,6) & rand(N,1)<f3);
                CurExposed1_O=((double(S_O)*double(infects)) & ~Exposed);
                CurExposed1_H=((double(S_H)*double(infects)) & ~Exposed);
                vecE_O=find(CurExposed1_O);
                vecE_H=find(CurExposed1_H);
                vecE1_O=rand(1,length(vecE_O))<p_NS;
                vecE2_O=rand(1,length(vecE_O))<1; % !!!!
                vecE1_H=rand(1,length(vecE_H))<p_NS;
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
                states(pstates(:,2)==1,4)=ceil((1/dt)*randWeibull(7,1.47,1,sum(pstates(:,2)==1)));
                states(pstates(:,3)>0,3)=pstates((pstates(:,3)>0),3)-1; % E_S
                states(pstates(:,3)==1,5)=ceil((1/dt)*randWeibull(3.3,1.47,1,sum(pstates(:,3)==1))); % !CHECK!    
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
                states(pstates(:,10)==1,11)=vecPV; %H2R
                states(pstates(:,10)==1,12)=~vecPV; %H2D    


    
                y(step,:)=y(step,:)+sum(~~states);

            end
            disp([rel,kT_out,kIn,K]);
            yRel=y-yP;
            yP=y;
            vecI=sum(yRel(:,6:8)')/N;
            [~,argI]=max(vecI);
            firs=find(vecI>0,1);
            x2=firs+ceil((argI-firs)/2);
            x1=firs+ceil((argI-firs)/4);
            if x2~=x1
            pV=polyfit((x1:x2)*dt,log(vecI(x1:x2)),1);
            if ~isnan(pV(1))
            betaR=betaR+pV(1);
            Rsum=Rsum+1;
            end
            end

        end
        betaR=betaR/Rsum;
        y=y/relalizations;
        betavec(kIn,kT_out)=betaR;
        I_in=I_in/relalizations;
        I_out=I_out/relalizations;
        f(kIn,kT_out)=I_in/(I_in+I_out);
        HTvec(kIn,kT_out)=mean(Hdiffr(Hdiffr>0));
    end
end
end
    end 
fCell{K}=f;
betavecCell{K}=betavec;
HTvecCell{K}=HTvec;
end
f=zeros(15,15);
fNAN=zeros(15,15);
betavec=zeros(15,15);
betavecNAN=zeros(15,15);


for k1=1:15
    fNAN=fNAN+(~isnan(fCell{k1}) & ~isinf(fCell{k1}) & fCell{k1}>0);
    fCell{k1}(isnan(fCell{k1}) | isinf(fCell{k1}) | fCell{k1}<0)=0;
    f=f+fCell{k1};
    betavecNAN=betavecNAN+(~isnan(betavecCell{k1}) & ~isinf(betavecCell{k1}) & betavecCell{k1}>0);
    betavecCell{k1}(isnan(betavecCell{k1}) | isinf(betavecCell{k1}) | betavecCell{k1}<0)=0;
    betavec=betavec+betavecCell{k1};

end
f=f./fNAN;
betavec=betavec./betavecNAN;

%% plots
% figure(4)
% subplot(1,2,1);
% hold on;
% plot(T_outvec,f','o','LineStyle','-');
% set(gca,'XScale','log')
% 
% subplot(1,2,2);
% hold on;
% plot(T_invec/4,f,'o','LineStyle','-');
% set(gca,'XScale','log')
% 
% 
% figure(5)
% subplot(1,2,1);
% hold on;
% plot(T_outvec,betavec','o','LineStyle','-');
% set(gca,'XScale','log')
% 
% subplot(1,2,2);
% hold on;
% plot(T_invec/4,betavec,'o','LineStyle','-');
% set(gca,'XScale','log')

% figure(6)
% subplot(1,2,1);
% hold on;
% plot(T_outvec,HTvec'/96,'o','LineStyle','-');
% set(gca,'XScale','log')
% 
% subplot(1,2,2);
% hold on;
% plot(T_invec/4,HTvec/96,'o','LineStyle','-');
% set(gca,'XScale','log')