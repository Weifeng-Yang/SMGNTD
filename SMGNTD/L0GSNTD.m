%%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions
function [core,var,loss,timerun]=L0GSNTD(core,var,ngmar,maxiteropt,Rdims,stopindex,r,alphat,btmax)
%% initialization algorithm
loss=[];
timerun=[0];
rate=0;
num=length(size(ngmar));
percent=[0.4,0.6,0.5,1];
percore=0.4;
for i=1:length(size(ngmar))
    a(i)=round(Rdims(i)*size(ngmar,i)*percent(i));
end
aa=a;
coreaa=ceil(prod(size(core))*percore);

Lapk=LLaplace(ngmar);
LK=zeros(1,num);
LCK=0;
L=ones(1,num);
LC=1;



for i=1:num-1
    alpha(i)=0;
end


alpha(num)=alphat;



for i=1:num
    varze{i}=zeros(size(var{i}));
end
coreze=zeros(size(core));




bts=[];
varK=var;
coreK=core;
tk=1;
tks=1;
bt=(tk-1)/tks;
wk=(tk-1)/(tk);

returnloss=norm(tensor(ngmar));
loss(1)=computeloss(ngmar,core,var,alpha,Lapk);


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


[core,var,LC,L,btc,bt]=L0GSNTDparaupdate(core,var,coreze,varze,num,ngmar,coreaa,aa,r,coretempK,vartempK,wk,L,LC,LtempK,LCtempK,alpha,Lapk,btmax);
loss(i+1)=computeloss(ngmar,core,var,alpha,Lapk);
%% Judging whether to extrapolate    
check=0;
for j=1:num
    check=check+norm(var{j}-varK{j},'fro')^2;
end
check=check+norm(tensor(core-coreK))^2;
rho=min([L,LC])/(1e+10);
if(loss(i+1)>=loss(i)-rho*check)
    var=varK;
    core=coreK;
    [core,var,LC,L,btc,bt]=L0GSNTDparaupdate(core,var,coreze,varze,num,ngmar,coreaa,aa,r,coretempK,vartempK,0,L,LC,LtempK,LCtempK,alpha,Lapk,btmax);
    loss(i+1)=computeloss(ngmar,core,var,alpha,Lapk);
end





%% Check if termination condition is met
fprintf("L0GSNTD\n");

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
Res=abs(loss(i+1)-loss(i))/returnloss;
% Res=check1/check2;
fprintf("cri：%d\n",Res);
stop=stopcheck(Res,timerun,stopindex);
if(stop==1)
    fprintf("Number of terminations：%d\n",i);
    pause(1);
    break;
end

tk=(1+sqrt(1+4*tk^2))/2;
wk=(tk-1)/(tks);
tks=tk;

end

end



function loss=computeloss(ngmar,core,var,alpha,Lapk)
    loss=compute(core,var,ngmar);
    for i=1:length(size(ngmar))
    loss=loss+alpha(i)/2*trace(var{i}'*Lapk{i}*var{i});
    end
end