function [ CellMatrix ] = temporalmatrix( A,dt,T,tauavg)

randExp=@(lambda,m,n) (-1/lambda).*log(1-rand(m,n));
[N,~]=size(A);
links=find(triu(A,1));
CellMatrix=cell(1,ceil(T/dt));
for k=(1:ceil(T/dt))
    CellMatrix{k}=sparse(N,N);
end

for k=1:length(links)
    vecT=ceil((1/dt)*mod((cumsum(randExp((1/tauavg),1,ceil(T/tauavg)))),T));
    [s,t]=ind2sub([N,N],links(k));
    for k1=1:length(vecT)
        CellMatrix{vecT(k1)}(s,t)=1;
        CellMatrix{vecT(k1)}(t,s)=1;
    end
    disp(k/length(links))
end

%temporal clique
for k=1:ceil(T/dt)
    vecC=find(sum(CellMatrix{k})>1);
    for k1=1:length(vecC)
        vecC1=find(CellMatrix{k}(:,vecC(k1)));
        CellMatrix{k}(vecC1,vecC1)=~eye(length(vecC1));
    end
end
end

