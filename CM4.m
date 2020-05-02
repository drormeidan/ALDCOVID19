function [yT,tT] = CM4(dt,t0,Tshift,Tweekend,f,rate,para1,para2,prob,NumberOfFlips,N,g1Choose,betaI,state,RQp)
%%% This is the code for the simulations accompanying this arXiv paper:
%%% https://arxiv.org/abs/2004.01453
%%% The notation here is consistent with the notation in the paper.

% Variables:
% yB holds the different states / compartment:
% S=yB(1),E_NS=yB(2),E_S=yB(3),I_NS=yB(4),I_PS=yB(5),I_M=yB(6),I_S=yB(7),I_C=yB(8),H=yB(9),V=yB(10),R=yB(11),D=yB(12)
% dP=[E_NS+,E_S+]
% prob=[p_NS,p_M,p_S,p_C,p_HR,p_HD,p_VR,p_VD]
% para1/2=[P1,P2,P3,P4];
% rate=[beta,r_MR,r_SH,r_CV,r_HR,r_HD,r_VR,r_VD]

% IMPORTANT NOTE: Tweekend<Tshift.
% state=0: Alterting quarantine.
% state=1: Intermitet quarantine.
% state=2: Half Quarantine (HQ).
% state=3: Full quarantine.
% state=4: Unmitigated.

% Flip variables
% S_LC=y(1),E_NS_LC=y(2),E_S_LC=y(3),I_NS_LC=y(4),I_PS_LC=y(5),I_M_LC=y(6),I_S_LC=y(7),I_C_LC=y(8),H_LC=y(9),V_LC=y(10),R_LC=y(11),D_LC=y(12)
% S_FC=y(13),E_NS_FC=y(14),E_S_FC=y(15),I_NS_FC=y(16),I_PS_FC=y(17),I_M_FC=y(18),I_S_FC=y(19),I_C_FC=y(20),H_FC=y(21),V_FC=y(22),R_FC=y(23),D_FC=y(24)
% S_D=y(25),E_NS_D=y(26),E_S_D=y(27),I_NS_D=y(28),I_PS_D=y(29),I_M_D=y(30),I_S_D=y(31),I_C_D=y(32),H_D=y(33),V_D=y(34),R_D=y(35),D_D=y(36)

InterT = t0/dt;
if (state == 4)
    t0 = t0+(NumberOfFlips*2*Tshift);
    NumberOfFlips = 0;
end

y0first = [N-2,1,1,0,0,0,0,0,0,0,0,0]; % N-1 susceptible, 1 ExposedNS, 1 ExposedS.
tB = (0:dt:(t0-dt))'; % time until intervention
yB = y0first;

if (g1Choose == 0)
    % preparing lognormal distribution
    g1 = @(x,mu,sigma) (1./(sqrt(2*pi)*sigma*x)).*exp(-((log(x)-mu).^2)/(2*(sigma.^2))); %lognormal
elseif (g1Choose == 1)
    % preparing Weibull distribution
    g1 = @(x,k,lambda) (k/lambda).*((x/lambda).^(k-1)).*(exp(-((x/lambda).^k))); %weibull
end

% Discretize the g function.
t = (dt:dt:20)';
gvec = zeros(length(t),4);
for kg = 1:4
    gvec(:,kg) = g1(t,para1(kg),para2(kg))./sum(g1(t,para1(kg),para2(kg))*dt);
end
decay = zeros(1,4);

% Plus variables (see SI)
dP = zeros(length(t),2);
dP = [dP;[y0first(2)/dt,y0first(3)/dt]];
atB = length(tB)+length(t);
plusm = [1,1,2,2];

for tidx = 2:length(tB)
    % introduce
    infects = yB(tidx-1,4)+yB(tidx-1,5);
    indexP = length(t)+tidx; 
    for kg = 1:4
        decay(kg) = sum(dP((indexP-1):-1:(indexP-length(gvec(:,kg))),plusm(kg)).*(gvec(:,kg))*dt);
    end
    
    if ((tidx>InterT) && (state == 4))
        rate(1) = betaI;
    end
    if (mod(tB(tidx),5) == 0)
        disp(tB(tidx));
    end
    
    % Euler
    yB(tidx,1) = yB(tidx-1,1)-dt*(rate(1)*(yB(tidx-1,1)/N)*infects); %S
    yB(tidx,2) = yB(tidx-1,2)+dt*((prob(1)*rate(1)*(yB(tidx-1,1)/N)*infects)-decay(1)); %E_NS
    yB(tidx,3) = yB(tidx-1,3)+dt*(((1-prob(1))*rate(1)*(yB(tidx-1,1)/N)*infects)-decay(3)); %E_S
    yB(tidx,4) = yB(tidx-1,4)+dt*(decay(1)-decay(2)); %I_NS
    yB(tidx,5) = yB(tidx-1,5)+dt*(decay(3)-decay(4)); %I_PS
    yB(tidx,6) = yB(tidx-1,6)+dt*((prob(2)*decay(4))-(rate(2)*yB(tidx-1,6))); %I_M
    yB(tidx,7) = yB(tidx-1,7)+dt*((prob(3)*decay(4))-(rate(3)*yB(tidx-1,7))); %I_S
    yB(tidx,8) = yB(tidx-1,8)+dt*((prob(4)*decay(4))-(rate(4)*yB(tidx-1,8))); %I_C
    yB(tidx,9) = yB(tidx-1,9)+dt*((rate(3)*yB(tidx-1,7))-(prob(5)*rate(5)*yB(tidx-1,9))-(prob(6)*rate(6)*yB(tidx-1,9))); %H
    yB(tidx,10) = yB(tidx-1,10)+dt*((rate(4)*yB(tidx-1,8))-(prob(7)*rate(7)*yB(tidx-1,10))-(prob(8)*rate(8)*yB(tidx-1,10))); %V
    yB(tidx,11) = yB(tidx-1,11)+dt*(decay(2)+(rate(2)*yB(tidx-1,6))+(prob(5)*rate(5)*yB(tidx-1,9))+(prob(7)*rate(7)*yB(tidx-1,10))); %R
    yB(tidx,12) = yB(tidx-1,12)+dt*((prob(6)*rate(6)*yB(tidx-1,9))+(prob(8)*rate(8)*yB(tidx-1,10))); %D

    % PlusDef
    dP(indexP,:) = [(prob(1)*rate(1)*(yB(tidx-1,1)/N)*infects),((1-prob(1))*rate(1)*(yB(tidx-1,1)/N)*infects)];
end
yT = yB;
tT = tB;

% Intervention
if (state == 0 || state == 1)
    y(1,:) = [(1/2)*(1-f)*yB(end,:),(1/2)*(1-f)*yB(end,:),(f)*yB(end,:)];
    dPLC = (1/2)*(1-f)*dP;
    dPFC = (1/2)*(1-f)*dP;
    dPD = (f)*dP;
end
if (state == 2)
    y(1,:) = [(RQp)*(1-f)*yB(end,:),(1-RQp)*(1-f)*yB(end,:),(f)*yB(end,:)];
    dPLC = (RQp)*(1-f)*dP;
    dPFC = (1-RQp)*(1-f)*dP;
    dPD = (f)*dP;
end
if (state == 3)
    state = 2; % State 3 (FQ) is implemented via State 2 (HQ) with cohot 1 being everyone
    y(1,:) = [(1-f)*yB(end,:),zeros(1,12),(f)*yB(end,:)];

    dPLC = (1-f)*dP;
    dPFC = 0*dP;
    dPD = (f)*dP;
end

decayLC = zeros(1,4);
decayFC = zeros(1,4);
decayD = zeros(1,4);
rate(1) = betaI;

for kFlips = 1:NumberOfFlips

    t = ((t0+((kFlips-1)*2*Tshift)):dt:(t0+((kFlips)*2*Tshift)-dt))'; %The time vector of two flips.

    dPLC(atB+((kFlips-1)*length(t))+1,:) = [0,0];
    dPFC(atB+((kFlips-1)*length(t))+1,:) = [0,0];
    dPD(atB+((kFlips-1)*length(t))+1,:) = [0,0];

    for tidx = 2:floor(length(t)/2) 

        indexP = atB+((kFlips-1)*length(t))+tidx;
        for kg = 1:4
            decayLC(kg) = sum(dPLC((indexP-1):-1:(indexP-length(gvec(:,kg))),plusm(kg)).*(gvec(:,kg))*dt);
            decayFC(kg) = sum(dPFC((indexP-1):-1:(indexP-length(gvec(:,kg))),plusm(kg)).*(gvec(:,kg))*dt); 
            decayD(kg) = sum(dPD((indexP-1):-1:(indexP-length(gvec(:,kg))),plusm(kg)).*(gvec(:,kg))*dt); 
        end

        %indoor
        if ((state == 1) && (tidx<(floor(length(t)/2)-floor(Tweekend/dt))))
            infects = (y(tidx-1,4)+y(tidx-1,5)+y(tidx-1,28)+y(tidx-1,29)+y(tidx-1,16)+y(tidx-1,17))/N;
            y(tidx,1) = y(tidx-1,1)-dt*(rate(1)*(y(tidx-1,1))*infects); %S
            y(tidx,2) = y(tidx-1,2)+dt*((prob(1)*rate(1)*(y(tidx-1,1))*infects)-decayLC(1)); %E_NS
            y(tidx,3) = y(tidx-1,3)+dt*(((1-prob(1))*rate(1)*(y(tidx-1,1))*infects)-decayLC(3)); %E_S
            dPLC(indexP,:) = [(prob(1)*rate(1)*(y(tidx-1,1))*infects),((1-prob(1))*rate(1)*(y(tidx-1,1))*infects)];
        else
            y(tidx,1) = y(tidx-1,1); %S
            y(tidx,2) = y(tidx-1,2)+dt*(-decayLC(1)); %E_NS
            y(tidx,3) = y(tidx-1,3)+dt*(-decayLC(3)); %E_S    
            dPLC(indexP,:) = [0,0];
        end
        y(tidx,4) = y(tidx-1,4)+dt*(decayLC(1)-decayLC(2)); %I_NS
        y(tidx,5) = y(tidx-1,5)+dt*(decayLC(3)-decayLC(4)); %I_PS
        y(tidx,6) = y(tidx-1,6)+dt*((prob(2)*decayLC(4))-(rate(2)*y(tidx-1,6))); %I_M
        y(tidx,7) = y(tidx-1,7)+dt*((prob(3)*decayLC(4))-(rate(3)*y(tidx-1,7))); %I_S
        y(tidx,8) = y(tidx-1,8)+dt*((prob(4)*decayLC(4))-(rate(4)*y(tidx-1,8))); %I_C
        y(tidx,9) = y(tidx-1,9)+dt*((rate(3)*y(tidx-1,7))-(prob(5)*rate(5)*y(tidx-1,9))-(prob(6)*rate(6)*y(tidx-1,9))); %H
        y(tidx,10) = y(tidx-1,10)+dt*((rate(4)*y(tidx-1,8))-(prob(7)*rate(7)*y(tidx-1,10))-(prob(8)*rate(8)*y(tidx-1,10))); %V
        y(tidx,11) = y(tidx-1,11)+dt*(decayLC(2)+(rate(2)*y(tidx-1,6))+(prob(5)*rate(5)*y(tidx-1,9))+(prob(7)*rate(7)*y(tidx-1,10))); %R
        y(tidx,12) = y(tidx-1,12)+dt*((prob(6)*rate(6)*y(tidx-1,9))+(prob(8)*rate(8)*y(tidx-1,10))); %D

        %outdoor
        if (tidx<(floor(length(t)/2)-floor(Tweekend/dt)) || state == 2)
            if (state ~= 1)
                infects = (y(tidx-1,16)+y(tidx-1,17)+y(tidx-1,28)+y(tidx-1,29))/N;
                y(tidx,13) = y(tidx-1,13)-dt*(rate(1)*(y(tidx-1,13))*infects); %S
                y(tidx,14) = y(tidx-1,14)+dt*((prob(1)*rate(1)*(y(tidx-1,13))*infects)-decayFC(1)); %E_NS
                y(tidx,15) = y(tidx-1,15)+dt*(((1-prob(1))*rate(1)*(y(tidx-1,13))*infects)-decayFC(3)); %E_S
                dPFC(indexP,:) = [(prob(1)*rate(1)*(y(tidx-1,13))*infects),((1-prob(1))*rate(1)*(y(tidx-1,13))*infects)];
            else
                infects = (y(tidx-1,4)+y(tidx-1,5)+y(tidx-1,28)+y(tidx-1,29)+y(tidx-1,16)+y(tidx-1,17))/N;    
                y(tidx,13) = y(tidx-1,13)-dt*(rate(1)*(y(tidx-1,13))*infects); %S
                y(tidx,14) = y(tidx-1,14)+dt*((prob(1)*rate(1)*(y(tidx-1,13))*infects)-decayFC(1)); %E_NS
                y(tidx,15) = y(tidx-1,15)+dt*(((1-prob(1))*rate(1)*(y(tidx-1,13))*infects)-decayFC(3)); %E_S
                dPFC(indexP,:) = [(prob(1)*rate(1)*(y(tidx-1,13))*infects),((1-prob(1))*rate(1)*(y(tidx-1,13))*infects)];        
            end    
        else
            if f>0
                infects = 0;
            end
            y(tidx,13) = y(tidx-1,13); %S
            y(tidx,14) = y(tidx-1,14)+dt*(-decayFC(1)); %E_NS
            y(tidx,15) = y(tidx-1,15)+dt*(-decayFC(3)); %E_S    
            dPFC(indexP,:) = [0,0];

        end
        y(tidx,16) = y(tidx-1,16)+dt*(decayFC(1)-decayFC(2)); %I_NS
        y(tidx,17) = y(tidx-1,17)+dt*(decayFC(3)-decayFC(4)); %I_PS
        y(tidx,18) = y(tidx-1,18)+dt*((prob(2)*decayFC(4))-(rate(2)*y(tidx-1,18))); %I_M
        y(tidx,19) = y(tidx-1,19)+dt*((prob(3)*decayFC(4))-(rate(3)*y(tidx-1,19))); %I_S
        y(tidx,20) = y(tidx-1,20)+dt*((prob(4)*decayFC(4))-(rate(4)*y(tidx-1,20))); %I_C
        y(tidx,21) = y(tidx-1,21)+dt*((rate(3)*y(tidx-1,19))-(prob(5)*rate(5)*y(tidx-1,21))-(prob(6)*rate(6)*y(tidx-1,21))); %H
        y(tidx,22) = y(tidx-1,22)+dt*((rate(4)*y(tidx-1,20))-(prob(7)*rate(7)*y(tidx-1,22))-(prob(8)*rate(8)*y(tidx-1,22))); %V
        y(tidx,23) = y(tidx-1,23)+dt*(decayFC(2)+(rate(2)*y(tidx-1,18))+(prob(5)*rate(5)*y(tidx-1,21))+(prob(7)*rate(7)*y(tidx-1,22))); %R
        y(tidx,24) = y(tidx-1,24)+dt*((prob(6)*rate(6)*y(tidx-1,21))+(prob(8)*rate(8)*y(tidx-1,22))); %D


        y(tidx,25) = y(tidx-1,25)-dt*(rate(1)*(y(tidx-1,25))*infects); %S
        y(tidx,26) = y(tidx-1,26)+dt*((prob(1)*rate(1)*(y(tidx-1,25))*infects)-decayD(1)); %E_NS
        y(tidx,27) = y(tidx-1,27)+dt*(((1-prob(1))*rate(1)*(y(tidx-1,25))*infects)-decayD(3)); %E_S
        y(tidx,28) = y(tidx-1,28)+dt*(decayD(1)-decayD(2)); %I_NS
        y(tidx,29) = y(tidx-1,29)+dt*(decayD(3)-decayD(4)); %I_PS
        y(tidx,30) = y(tidx-1,30)+dt*((prob(2)*decayD(4))-(rate(2)*y(tidx-1,30))); %I_M
        y(tidx,31) = y(tidx-1,31)+dt*((prob(3)*decayD(4))-(rate(3)*y(tidx-1,31))); %I_S
        y(tidx,32) = y(tidx-1,32)+dt*((prob(4)*decayD(4))-(rate(4)*y(tidx-1,32))); %I_C
        y(tidx,33) = y(tidx-1,33)+dt*((rate(3)*y(tidx-1,31))-(prob(5)*rate(5)*y(tidx-1,33))-(prob(6)*rate(6)*y(tidx-1,33))); %H
        y(tidx,34) = y(tidx-1,34)+dt*((rate(4)*y(tidx-1,32))-(prob(7)*rate(7)*y(tidx-1,34))-(prob(8)*rate(8)*y(tidx-1,34))); %V
        y(tidx,35) = y(tidx-1,35)+dt*(decayD(2)+(rate(2)*y(tidx-1,30))+(prob(5)*rate(5)*y(tidx-1,33))+(prob(7)*rate(7)*y(tidx-1,34))); %R
        y(tidx,36) = y(tidx-1,36)+dt*((prob(6)*rate(6)*y(tidx-1,33))+(prob(8)*rate(8)*y(tidx-1,34))); %D

        dPD(indexP,:) = [(prob(1)*rate(1)*(y(tidx-1,25))*infects),((1-prob(1))*rate(1)*(y(tidx-1,25))*infects)];
    end
    
    % and Flip!
    for tidx = (floor(length(t)/2)+1):length(t) 
        indexP = atB+((kFlips-1)*length(t))+tidx;
        for kg = 1:4
            decayLC(kg) = sum(dPLC((indexP-1):-1:(indexP-length(gvec(:,kg))),plusm(kg)).*(gvec(:,kg))*dt); 
            decayFC(kg) = sum(dPFC((indexP-1):-1:(indexP-length(gvec(:,kg))),plusm(kg)).*(gvec(:,kg))*dt); 
            decayD(kg) = sum(dPD((indexP-1):-1:(indexP-length(gvec(:,kg))),plusm(kg)).*(gvec(:,kg))*dt); 
        end

    
    %indoor
    if (state == 2)
        infects = (y(tidx-1,16)+y(tidx-1,17)+y(tidx-1,28)+y(tidx-1,29))/N;
        y(tidx,13) = y(tidx-1,13)-dt*(rate(1)*(y(tidx-1,13))*infects); %S
        y(tidx,14) = y(tidx-1,14)+dt*((prob(1)*rate(1)*(y(tidx-1,13))*infects)-decayFC(1)); %E_NS
        y(tidx,15) = y(tidx-1,15)+dt*(((1-prob(1))*rate(1)*(y(tidx-1,13))*infects)-decayFC(3)); %E_S
        dPFC(indexP,:) = [(prob(1)*rate(1)*(y(tidx-1,13))*infects),((1-prob(1))*rate(1)*(y(tidx-1,13))*infects)];        
    else
        y(tidx,13) = y(tidx-1,13); %S
        y(tidx,14) = y(tidx-1,14)+dt*(-decayFC(1)); %E_NS
        y(tidx,15) = y(tidx-1,15)+dt*(-decayFC(3)); %E_S
        dPFC(indexP,:) = [0,0];    
    end
    y(tidx,16) = y(tidx-1,16)+dt*(decayFC(1)-decayFC(2)); %I_NS
    y(tidx,17) = y(tidx-1,17)+dt*(decayFC(3)-decayFC(4)); %I_PS
    y(tidx,18) = y(tidx-1,18)+dt*((prob(2)*decayFC(4))-(rate(2)*y(tidx-1,18))); %I_M
    y(tidx,19) = y(tidx-1,19)+dt*((prob(3)*decayFC(4))-(rate(3)*y(tidx-1,19))); %I_S
    y(tidx,20) = y(tidx-1,20)+dt*((prob(4)*decayFC(4))-(rate(4)*y(tidx-1,20))); %I_C
    y(tidx,21) = y(tidx-1,21)+dt*((rate(3)*y(tidx-1,19))-(prob(5)*rate(5)*y(tidx-1,21))-(prob(6)*rate(6)*y(tidx-1,21))); %H
    y(tidx,22) = y(tidx-1,22)+dt*((rate(4)*y(tidx-1,20))-(prob(7)*rate(7)*y(tidx-1,22))-(prob(8)*rate(8)*y(tidx-1,22))); %V
    y(tidx,23) = y(tidx-1,23)+dt*(decayFC(2)+(rate(2)*y(tidx-1,18))+(prob(5)*rate(5)*y(tidx-1,21))+(prob(7)*rate(7)*y(tidx-1,22))); %R
    y(tidx,24) = y(tidx-1,24)+dt*((prob(6)*rate(6)*y(tidx-1,21))+(prob(8)*rate(8)*y(tidx-1,22))); %D


    %outdoor
    if (tidx<(length(t)-floor(Tweekend/dt)) && state == 0)
        infects = (y(tidx-1,4)+y(tidx-1,5)+y(tidx-1,28)+y(tidx-1,29))/N;
        y(tidx,1) = y(tidx-1,1)-dt*(rate(1)*(y(tidx-1,1))*infects); %S
        y(tidx,2) = y(tidx-1,2)+dt*((prob(1)*rate(1)*(y(tidx-1,1))*infects)-decayLC(1)); %E_NS
        y(tidx,3) = y(tidx-1,3)+dt*(((1-prob(1))*rate(1)*(y(tidx-1,1))*infects)-decayLC(3)); %E_S
        dPLC(indexP,:) = [(prob(1)*rate(1)*(y(tidx-1,1))*infects),((1-prob(1))*rate(1)*(y(tidx-1,1))*infects)];
    else
        if (state ~= 2)
            if f>0 && (tidx<(length(t)-floor(Tweekend/dt)))
                infects = (y(tidx-1,28)+y(tidx-1,29))/N;
            else
                infects = 0;
            end
        else
            infects = (y(tidx-1,16)+y(tidx-1,17)+y(tidx-1,28)+y(tidx-1,29))/N;        
        end
        y(tidx,1) = y(tidx-1,1); %S
        y(tidx,2) = y(tidx-1,2)+dt*(-decayLC(1)); %E_NS
        y(tidx,3) = y(tidx-1,3)+dt*(-decayLC(3)); %E_S    
        dPLC(indexP,:) = [0,0];    
    end
    y(tidx,4) = y(tidx-1,4)+dt*(decayLC(1)-decayLC(2)); %I_NS
    y(tidx,5) = y(tidx-1,5)+dt*(decayLC(3)-decayLC(4)); %I_PS
    y(tidx,6) = y(tidx-1,6)+dt*((prob(2)*decayLC(4))-(rate(2)*y(tidx-1,6))); %I_M
    y(tidx,7) = y(tidx-1,7)+dt*((prob(3)*decayLC(4))-(rate(3)*y(tidx-1,7))); %I_S
    y(tidx,8) = y(tidx-1,8)+dt*((prob(4)*decayLC(4))-(rate(4)*y(tidx-1,8))); %I_C
    y(tidx,9) = y(tidx-1,9)+dt*((rate(3)*y(tidx-1,7))-(prob(5)*rate(5)*y(tidx-1,9))-(prob(6)*rate(6)*y(tidx-1,9))); %H
    y(tidx,10) = y(tidx-1,10)+dt*((rate(4)*y(tidx-1,8))-(prob(7)*rate(7)*y(tidx-1,10))-(prob(8)*rate(8)*y(tidx-1,10))); %V
    y(tidx,11) = y(tidx-1,11)+dt*(decayLC(2)+(rate(2)*y(tidx-1,6))+(prob(5)*rate(5)*y(tidx-1,9))+(prob(7)*rate(7)*y(tidx-1,10))); %R
    y(tidx,12) = y(tidx-1,12)+dt*((prob(6)*rate(6)*y(tidx-1,9))+(prob(8)*rate(8)*y(tidx-1,10))); %D


    y(tidx,25) = y(tidx-1,25)-dt*(rate(1)*(y(tidx-1,25))*infects); %S
    y(tidx,26) = y(tidx-1,26)+dt*((prob(1)*rate(1)*(y(tidx-1,25))*infects)-decayD(1)); %E_NS
    y(tidx,27) = y(tidx-1,27)+dt*(((1-prob(1))*rate(1)*(y(tidx-1,25))*infects)-decayD(3)); %E_S
    y(tidx,28) = y(tidx-1,28)+dt*(decayD(1)-decayD(2)); %I_NS
    y(tidx,29) = y(tidx-1,29)+dt*(decayD(3)-decayD(4)); %I_PS
    y(tidx,30) = y(tidx-1,30)+dt*((prob(2)*decayD(4))-(rate(2)*y(tidx-1,30))); %I_M
    y(tidx,31) = y(tidx-1,31)+dt*((prob(3)*decayD(4))-(rate(3)*y(tidx-1,31))); %I_S
    y(tidx,32) = y(tidx-1,32)+dt*((prob(4)*decayD(4))-(rate(4)*y(tidx-1,32))); %I_C
    y(tidx,33) = y(tidx-1,33)+dt*((rate(3)*y(tidx-1,31))-(prob(5)*rate(5)*y(tidx-1,33))-(prob(6)*rate(6)*y(tidx-1,33))); %H
    y(tidx,34) = y(tidx-1,34)+dt*((rate(4)*y(tidx-1,32))-(prob(7)*rate(7)*y(tidx-1,34))-(prob(8)*rate(8)*y(tidx-1,34))); %V
    y(tidx,35) = y(tidx-1,35)+dt*(decayD(2)+(rate(2)*y(tidx-1,30))+(prob(5)*rate(5)*y(tidx-1,33))+(prob(7)*rate(7)*y(tidx-1,34))); %R
    y(tidx,36) = y(tidx-1,36)+dt*((prob(6)*rate(6)*y(tidx-1,33))+(prob(8)*rate(8)*y(tidx-1,34))); %D

    dPD(indexP,:) = [(prob(1)*rate(1)*(y(tidx-1,25))*infects),((1-prob(1))*rate(1)*(y(tidx-1,25))*infects)];
        
    end
    
    % we chose not to prealocate yT and tT in favor of simplicity and readability
    yT = [yT;y(:,1:12)+y(:,13:24)+y(:,25:36)];
    tT = [tT;t];
    y = y(end,:);
    disp([kFlips,f,t0,state,RQp]);
end % end of the flip loop

end % funtion

