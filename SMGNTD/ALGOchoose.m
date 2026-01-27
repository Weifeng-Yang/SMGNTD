
% Input.
% var,core    : initial matrix and core tensor
% ngmar       : decomposed tensor
% The remaining parameters are explained the same as the 'main_Run_me' function

% Output.
% cores, vars : Decomposition matrix and core tensor resulting from the final iterative result
% loss:       : Array of loss functions generated during iteration
% tr:         : Runtime array during iteration

function [data,varss]=ALGOchoose(core,var,ngmar,maxiteropt,Rdims,flag,stopindex,r,alphat,lamda,index)

if(flag==1)
[cores,vars,loss,tr]=SNTD(core,var,ngmar,maxiteropt,stopindex,r,0.5);
varss{1}=cores;
varss{2}=vars;
lossdata=loss;
trdata=tr;
 
  
 
  

elseif(flag==2) 
    
[cores,vars,loss,tr]=ARTD(core,var,ngmar,maxiteropt,stopindex,r,0.5);
varss{1}=cores;
varss{2}=vars;
lossdata=loss;
trdata=tr;
 
  

elseif(flag==3) 
    
[cores,vars,loss,tr]=GSNTD(core,var,ngmar,maxiteropt,Rdims,stopindex,r,1.7,10);
varss{1}=cores;
varss{2}=vars;
lossdata=loss;
trdata=tr;
 
 
  

elseif(flag==4) 
[cores,vars,loss,tr]=L0GSNTD(core,var,ngmar,maxiteropt,Rdims,stopindex,r,1,0.99);  
varss{1}=cores;
varss{2}=vars;
lossdata=loss;
trdata=tr;



elseif(flag==5) 
[cores,vars,loss,tr]=MGNTD(core,var,ngmar,maxiteropt,stopindex,r,index);
varss{1}=cores;
varss{2}=vars;
lossdata=loss;
trdata=tr;

elseif(flag==6) 
[cores,vars,loss,tr]=AMGRNTD(core,var,ngmar,maxiteropt,stopindex,r);
varss{1}=cores;
varss{2}=vars;
lossdata=loss;
trdata=tr;




elseif(flag==7) 
[cores,vars,loss,tr]=SMGNTD(core,var,ngmar,maxiteropt,stopindex,r,alphat,lamda,0.99);
varss{1}=cores;
varss{2}=vars;
lossdata=loss;
trdata=tr;



 
  
end


data{1}=lossdata;
data{2}=trdata;






end



