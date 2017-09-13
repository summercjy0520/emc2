function [SSresult]=CrossEmerge(SS,index,kuser)
Kcurr=size(SS,2);
% SSresult=cell(Kcurr-kuser,1);
% A=1:Kcurr;
f=1;
while Kcurr>kuser
    DSC=zeros(Kcurr,1);
    for i=1:Kcurr
        SSt=SS;
%         ii=rhomaxindex(index,ordrho);
        SSt(i)=[];
        Jsc=0;
        for j=1:Kcurr-2
            for h=j+1:Kcurr-1
                S1=SSt{j};
                S2=SSt{h};
                q=size(S1,2);

                Jsc=Jsc+ComputeCS(S1,S2,q);
%                 Dsc=Dsc-log(Jsc);
            end            
        end
        DSC(i)=-log((2/((Kcurr-1)*(Kcurr-2)))*Jsc);       
    end
    [~,d]=max(DSC); %min(JSC) 
    
    p=0.1;
    zz=size(SS,2);
    SSt=SS;
    SSt(d)=[];
    ind=index;
    ind(d)=[];
    pdist=zeros(size(index,1)-1,1);
    for i=1:zz-1
        pdist(i)=computepd(SS{d},SSt{i},p);
    end 
    [~,aa]=min(pdist);
    xx=find(index==ind(aa));
    SS{xx}=[SSt{aa};SS{d}];
    SS(d)=[];  
    index(d)=[];
    
    f=f+1;
    Kcurr=Kcurr-1;
end

SSresult=SS;
end


function dist=computepd(A,B,p)
k=ceil(min(size(A,1)*p,size(B,1)));
stack=inf*ones(k,1);
for i=1:size(A,1)
    distance=pdist2(B,A(i,:));
    if max(stack)>min(distance)
        [~,ind]=max(stack);
        stack(ind)=min(distance);
    end
end
dist=mean(stack);
end