function [node]=deltaEntropyD(Yra,rho)
[Yra,~]=sort(Yra,'descend');
Y=zeros(size(Yra,1)-1,1);
% Y(1)=(Yra(1)-Yra(2))/2;
for i=2:size(Yra,1)
       Y(i-1)=abs(Yra(i)*rho(i)-Yra(1)*rho(1));
%     Y(i-1)=abs(Yra(i)-Yra(1))*rho(i);
end
% Y=sort(Yra);
% Y=pdist2([Yra,rho],[max(Yra),max(rho)],'cosine');
% Y=1-Y;
maxgain=inf; %-inf
diff=zeros(size(Y,1),1);
% entropy0=computeEntropy(Y);
for i=2:size(Y,1)
    gain=computeConditionalEntropy(Y,i);
    if(gain<=maxgain)
        maxgain=gain;
        diff(i)=gain;%DiffEntropy(Y,i);
    else
        diff(i)=gain;
        node=i;
%         break;
    end
end
% nn=size(find(diff==0),1);
% diff(ismember(diff,0,'rows'))=[];
% if isempty(diff)
%     node=1;
% else
%     [~,node]=min(diff);
%     node=node+nn;
% end
% [c,~]=max(diff);
% K=entropy0-diff(node);
% K=diff(node);
[~,node]=min(diff(2:end));
end

function entropy=computeEntropyL(dataset)
entropy=0;
AU=size(dataset,1);
for i=2:AU
%     if dataset(i)>0
    if dataset(i-1)-dataset(i)>0
         prob=(dataset(i)-dataset(i-1))/(dataset(i)-dataset(1));
%         prob=(dataset(i)-dataset(i-1))/(dataset(i));%(dataset(end)-dataset(1)); 
        entropy=entropy+(-1)*prob*log(prob)/log(2);
    else
        entropy=entropy+0;
    end
end
end

function entropy=computeEntropyR(dataset)
entropy=0;
AU=size(dataset,1);
for i=2:AU
%     if dataset(i)>0
    if dataset(i)-dataset(i-1)>0
        prob=(dataset(i)-dataset(i-1))/(dataset(i)-dataset(1));
%          prob=(dataset(i)-dataset(i-1))/(dataset(end)-dataset(i));
%           prob=(dataset(i)-dataset(i-1))/(dataset(i));%(dataset(end)-dataset(1)); 
        entropy=entropy+(-1)*prob*log(prob)/log(2);
    else
        entropy=entropy+0;
    end
end
end

function conditionalentropy=computeConditionalEntropy(dataset,m)
% conditionalentropy=0;
probx=m/size(dataset,1);
proby=(size(dataset,1)-m+1)/size(dataset,1);

dataset1=dataset(1:m);
dataset2=dataset(m+1:end);

probEntropyx=computeEntropyL(dataset1);
probEntropyy=computeEntropyR(dataset2);

conditionalentropy=probx*probEntropyx+proby*probEntropyy;
end
