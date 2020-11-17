function [ S_O,S_H ] = closeS( S, household,sizeH, closeHomeQ,closeHomeI,fQ,fI,Nodef1,Nodef3,closeHomeQnode,HouseholdNet)

[source,target,~]=find(S);
[N,~]=size(S);
% isolate
closeNodeQ= closeHomeQnode & ~Nodef1;

closeNodeI=zeros(1,N);
for k=1:length(closeHomeI)
    closeNodeI(household{closeHomeI(k)})=ones(1,length(household{closeHomeI(k)}));
end
closeNodeI= closeNodeI & ~Nodef3;

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

