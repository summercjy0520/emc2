function [node]=emc2EntropyDD1(Y)
%%% pairdist computing
[Y,KT]=sort(Y);
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
%         node=i;
%         break;
    end
end
[~,node]=min(diff(2:end));
% K=entropy0-diff(node);
% K=diff(node);
% nn=size(find(diff==0),1);
% diff(ismember(diff,0,'rows'))=[];
% if isempty(diff)
%     node=1;
% else
%     [~,node]=min(diff);
%     node=node+nn;
% end
end

function entropy=computeEntropyL(dataset)
entropy=0;
AU=size(dataset,1);
for i=2:AU
%     if dataset(i)>0
    if dataset(i)-dataset(i-1)>0
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
