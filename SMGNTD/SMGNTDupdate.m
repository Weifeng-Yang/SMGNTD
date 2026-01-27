

function [core,var,LC,L,btc,bt]=SMGNTDupdate(core,var,coreK,varK,num,ngmar,r,wk,L,LC,LK,LCK,Lapk,alpha,lamda,btmax)

    btc=min(wk,btmax);
    core=core+btc*(core-coreK);
    [V,LC]=gradcore(core,var,ngmar,r,num);
    core= PROXL1(V,lamda(1),1/(LC*r));
    
    for j=1:num
    bt(j)=min(wk,btmax);
    var{j}=var{j}+bt(j)*(var{j}-varK{j});
    [V,L(j)]=gradMGNTD(core,var,ngmar,r,j,num,Lapk,alpha);
    if(j==num)
        var{j}=PROXn1(V);
    else
        var{j}= PROXL1(V,lamda(j+1),1/(L(j)*r));
    end
    end

end


function [V,L]=gradMGNTD(core,var,ngmar,r,n,num,Lapk,alpha)
core=tensor(core);
index=1:num;
index(n)=[];
coreg=ttm(core,var,index);
tempB=double(tenmat(coreg,n));
temp=tempB*tempB';
Xn=double(tenmat(ngmar,n));
if(n==num)
    Laptemp=alpha(1)*Lapk{1}*var{n};
    L=alpha(1)*norm(Lapk{1},'fro');
    for j=2:length(Lapk)
        Laptemp=Laptemp+alpha(j)*Lapk{j}*var{n};
        L=L+alpha(j)*norm(Lapk{j},'fro');
    end
    U=var{n}*temp-Xn*tempB'+Laptemp;
    L=norm(temp,'fro')+L;
    
else
    U=var{n}*temp-Xn*tempB';
    L=norm(temp,'fro');
end
V=var{n}-1/(r*L)*U;
end

function x=PROXn1(x)
x=double(x);
x(x<0)=0;
end


