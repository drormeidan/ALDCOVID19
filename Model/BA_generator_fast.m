function [ A,A_H,A_O,Household,Node2HH] = BA_generator_fast( NumberOfHouslehold,parameter,HHSprob,real)
% generate network with households.. 
% HHSprob=[p1,p2,p3,p4,p5,p6] 
vecP=cumsum(HHSprob);
HHvec=zeros(1,NumberOfHouslehold);
Household=cell(1,NumberOfHouslehold);
for k=1:NumberOfHouslehold
    HHvec(k)=find(rand<vecP,1);
end
N=sum(HHvec);
Node2HH=zeros(1,N);
csHHvec=cumsum(HHvec);
csHHvec=[0,csHHvec];
vecHH=randperm(N);
vecDeg=floor(PLVgenerator(3,1,N,2*parameter));
if real==1
    A1=cm_net(vecDeg);
else
    kavgvec=15;
    p=kavgvec/(N-1);
%     A1=zeros(N);
%     for k=1:N-1
%     A1=A1+diag((rand(1,N-k)<p),k);
%     disp(k)
%     end
% A1=A1'+A1;
A1=rand(N);
A1=A1<p;
A1=triu(A1,1)+triu(A1,1)';

end
A2=zeros(N);
for k=1:(NumberOfHouslehold)
    vecHouse=vecHH((csHHvec(k)+1):csHHvec(k+1));
    Node2HH(vecHouse)=k*ones(1,length(vecHouse));
    Household{k}=vecHouse;
    A2(vecHouse,vecHouse)=ones(HHvec(k));
    disp(k);
end
A=sparse((A1 | A2) & ~eye(N));
A_O=sparse(A1);
A_H=sparse(A2 & ~eye(N));
end
