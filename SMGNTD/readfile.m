function [ngmar,R,Rdims,tlabel,ngmar1]=readfile(i)
  if(i==1) 
   E=load('.\Data\FashionMnist.mat');
   ngmar= double(E.fea);
   ngmart=ngmar;
   sizea=sqrt(size(ngmar,1));
   ngmar=reshape(ngmar,sizea,sizea,prod(size(ngmar))/sizea^2);
   ngmar1=ngmar;
   label= double(E.gnd);


    M1=tenmat(ngmar,3);
    M11=double(M1)';
    M11=normalize(M11,'range');
    M11=M11';
    rdims =M1.rdims;% Dimensions that were mapped to the rows.
    cdims =M1.cdims;% Dimensions that were mapped to the columns.
    tsize = size(ngmar);
    M10=tenmat(M11,rdims,cdims,tsize);
%     M10=tenmat(M11,length(size(ngmar)),cdims-1,fliplr(tsize));
    ngmar=tensor(M10);


   tlabel=double(label);
    R=length(unique(tlabel));


      if(find(tlabel==0)~=0)
            tlabel=tlabel(:,1)+1;
       end
       

%        cluster = kmeans(M11, length(unique(label))); 
%      NMIs(1)=NMI(label,cluster);

    %     sum1=0;
%     for j=1:length(size(ngmar))
%         sum1=sum1+norm(ngmar(:,:,j)-tensor(squeeze(ngmar1(j,:,:))));
%     end

    len=100;
    labelu=unique(tlabel);
   
    label=[];
    id=[];
     for j=1:R
         temp=find(tlabel==labelu(j));
         id=[id;temp(1:len)];
         label=[label,ones(1,len)*j];
     end
    
     tlabel=tlabel(id);
     ngmart=ngmar(:,:,id);
     ngmar=ngmart;
     ngmar1=ngmar1(:,:,id);


     tsize=size(ngmar);
     Rdims=[ceil(tsize(1:end-1)/2),R];






  elseif(i==2)   
   E=load('.\Data\pixraw10P.mat');
   ngmar= double(E.X);
   ngmart=ngmar;
   ngmar=reshape(ngmar',100,100,size(ngmar,1));
   ngmar1=ngmar;
   label= double(E.Y);

    M1=tenmat(ngmar,3);
    M11=double(M1)';
    M11=normalize(M11,'range');
    M11=M11';
    rdims =M1.rdims;% Dimensions that were mapped to the rows.
    cdims =M1.cdims;% Dimensions that were mapped to the columns.
    tsize = size(ngmar); %或M1.tsize % Size of the original tensor.
    M10=tenmat(M11,rdims,cdims,tsize);
    ngmar=tensor(M10);

   tlabel=double(label);
    R=length(unique(tlabel));

     tsize=size(ngmar);
     Rdims=[ceil(tsize(1:end-1)/2),R];


  end



end


