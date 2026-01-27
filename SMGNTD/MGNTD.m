%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions

function [core,var,loss,timerun]=MGNTD(core,var,ngmar,maxiteropt,stopindex,r,kiter)
%% initialization algorithm
loss=[];
timerun=[0];
num=length(size(ngmar));
LK=zeros(1,num);
LCK=0;
L=ones(1,num);
LC=1;
tk=1;


varK=var;
coreK=core;
wk=(tk-1)/(tk);

returnloss=norm(tensor(ngmar));
alpha=[];

[~,~,~,~,ngmar1]=readfile(kiter);
Lapktemp=LLaplace(ngmar1);
Lapk{1}=Lapktemp{end};

Lapktemp=LLaplace(ngmar);
Lapk{2}=Lapktemp{end};
clearvars Lapktemp

for j=1:length(Lapk)
    alpha(j)=1;
end
loss(1)=computeloss(ngmar,core,var,alpha,Lapk);


rate=[0];
t1=clock;


for i=1:maxiteropt
%% update parameters
fprintf("%d\n",i);

LtempK=LK;
LCtempK=LCK;
vartempK=varK;
coretempK=coreK;    
varK=var;
coreK=core;
LK=L;
LCK=LC;    
[core,var,LC,L,btc,bt]=MGNTDupdate(core,var,coretempK,vartempK,num,ngmar,r,wk,L,LC,LtempK,LCtempK,Lapk,alpha);
loss(i+1)=computeloss(ngmar,core,var,alpha,Lapk);


%% Judging whether to extrapolate
if(loss(i+1)>=loss(i))
    var=varK;
    core=coreK;
    L=LK;
    LC=LCK;
    [core,var,LC,L]=MGNTDupdate(core,var,coretempK,vartempK,num,ngmar,r,0,L,LC,LtempK,LCtempK,Lapk,alpha);
    loss(i+1)=computeloss(ngmar,core,var,alpha,Lapk);
end


%% Check if termination condition is met
fprintf("MGNTD\n");
fprintf("nonzero:%d\n",nnz(core));  

check1=norm(tensor(core));
check2=norm(tensor(coreK));
for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}));  
    check1=check1+norm(var{j}-varK{j},'fro');
    check2=check2+norm(varK{j},'fro');
end

bts{i}=[btc,bt];
t2=clock;
timerun(i+1)=etime(t2,t1);
% Res=abs(loss(i+1)-loss(i));
Res=check1/check2;
fprintf("cri：%d\n",Res);
stop=stopcheck(Res,timerun,stopindex);
if(stop==1)
    fprintf("Number of terminations：%d\n",i);
    pause(4);
    break;
end

tk=(1+sqrt(1+4*tk^2))/2;
wk=(tk-1)/(tk);



end
end


function loss=computeloss(ngmar,core,var,alpha,Lapk)
    loss=compute(core,var,ngmar);
    num=length(var);
    for i=1:length(Lapk)
    loss=loss+alpha(i)/2*trace(var{num}'*Lapk{i}*var{num});
    end
end
