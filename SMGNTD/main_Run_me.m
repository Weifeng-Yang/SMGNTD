clearvars -except 
clc
warning('off');
%% Parameter.
%   index     : The dataset to be used, when index=1, use FashionMnist dataset
%               when index=2, use pixraw10P dataset. 
%   r         : Step factor
%   maxiteropt: Maximum iteration alloted to the method
%   trigger   : Whether to enable the indicator array of each method, where
%               when 1∈trigger, enable the 𝓁1-SNTD method
%               when 2∈trigger, enable the ARTD method
%               when 3∈trigger, enable the GSNTD method
%               when 4∈trigger, enable the 𝓁0-GSNTD method
%               when 5∈trigger, enable the MGNTD method
%               when 6∈trigger, enable the AMGRNTD method
%               when 7∈trigger, enable the SMGNTD method 
%   lamda     ：The sparse regularization parameters
%   alphat    : The graph regularization parameters
%   stopindex : The indicator of the stop condition.  
%               To set the specific termination condition, see the 'stopcheck' function for details.  
%               The default termination condition is: ϵ<1e-5 or maxiteropt>8000 
%% Display
%   nonzero   ：The number of non-zero elements in each component.
%   Rel       : The difference in the variable value between two iterations.

%% Parameter settings
rng('shuffle')
index=[1];
r=1.01;
maxiteropt=8000;
trigger=[1,2,3,4,5,6,7];
alphat=[1,1,1];
lamda=[1,1,1];
stopindex=4;




%% Select dataset
[ngmar,R,Rdims,y]=readfile(index);
num=length(size(ngmar));
N=R;
X1=double(tenmat(ngmar,length(size(ngmar))))';






for j=1:10
%% Init
var=[];
core=[];
for i=1:num
    var{i}=rand(size(ngmar,i),Rdims(i));
end
core=tenrand(Rdims);
core=tensor(core);






%% Solving
for i=1:length(trigger)
[datas{i},varss{i}]=ALGOchoose(core,var,ngmar,maxiteropt,Rdims,trigger(i),stopindex,r,alphat,lamda,index);
end
vart=var;
vart{num+1}=core;
datas{length(trigger)+1}=varss;
datas{length(trigger)+2}=vart;
for i=1:length(trigger)
    vars1=datas{length(trigger)+1};
    vars11=vars1{i};
    vartemp=vars11{end};

    if(trigger(i)==3)
        cluster=ttm(tensor(vars11{1}),vartemp{end},num);
        temp=double(tenmat(cluster,num));
        [acc(j,i),rdx(j,i),NMIs(j,i)]=clustermeans(temp,N,y);

    else
         [acc(j,i),rdx(j,i),NMIs(j,i)]=clustermeans(vartemp{end},N,y);
    end

end
    datass{j}=datas;

end

%% Display results
accmean=mean(acc)
rdxmean=mean(rdx)
nmimean=mean(NMIs)
plt0=barplot(trigger,acc,rdx,NMIs); plt0=barplot(trigger,acc,rdx,NMIs)


