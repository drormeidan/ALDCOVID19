function [ S_O,S_H ] = closeS_grid( S, household,sizeH, closeHomeQ,closeHomeI,fQ,fI,HouseholdNet)

[source,target,~]=find(S);
[N,~]=size(S);
% isolate
closeNodeQ=zeros(1,N);
for k=1:length(closeHomeQ)
    closeNodeQ(household{closeHomeQ(k)})=rand(1,length(household{closeHomeQ(k)}))>fQ;
end

closeNodeI=zeros(1,N);
for k=1:length(closeHomeI)
    closeNodeI(household{closeHomeI(k)})=rand(1,length(household{closeHomeI(k)}))>fI;
end

closeNode=find(closeNodeI | closeNodeQ);
for k=1:length(closeNode)
    node=closeNode(k);
    vec=(node==source) | (node==target);
    if (sum(vec)>0)
        source(vec)=[];
        target(vec)=[];
    end
end

S_O=sparse(source,target,ones(1,length(source)),N,N);
% households
closeHome= unique([closeHomeQ,closeHomeI]);
sizeCH=sizeH(closeHome);
sourceH=zeros(1,sum(sizeCH.*(sizeCH-1)));
targetH=zeros(1,sum(sizeCH.*(sizeCH-1)));
flag=1;
for k=1:length(closeHome)
    [s,t]=find(HouseholdNet{closeHome(k)});
    for k1=1:length(s)
    sourceH(flag)=household{closeHome(k)}(s(k1));
    targetH(flag)=household{closeHome(k)}(t(k1));
    flag=flag+1;
    end
end
sourceH=sourceH(sourceH~=0);
targetH=targetH(targetH~=0);
if flag>1
    S_H=sparse(sourceH,targetH,ones(1,length(sourceH)),N,N);
else
    S_H=sparse(N,N);
end


end

