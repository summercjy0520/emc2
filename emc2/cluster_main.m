%%%%%% maybe the last banben
%%%%%%%%%%%%%%  add the points from area 1,and then merge the samll areas.
%%%%%%%%%%%rhok > number(SS{i})

%%%%%%%%%%%%%cross-entropy

clear all
close all
clc
disp('The only input needed is a distance matrix file')
disp('The format of this file should be: ')
disp('Column 1: id of element i')
disp('Column 2: id of element j')
disp('Column 3: dist(i,j)')

disp('Reading input distance matrix')

TXT= importdata('C:\Documents and Settings\Administrator\×ÀÃæ\cjy\2017\7\data\Jain.txt');
S=TXT(:,1:end-1);
xx=Computedist(S);
% load('S1.mat');

ND=max(xx(:,2));
NL=max(xx(:,1));
if (NL>ND)
  ND=NL;
end
N=size(xx,1);
for i=1:ND
  for j=1:ND
    dist(i,j)=0;
  end
end
for i=1:N
  ii=xx(i,1);
  jj=xx(i,2);
  dist(ii,jj)=xx(i,3);
  dist(jj,ii)=xx(i,3);
end
percent=45.0;
fprintf('average percentage of neighbours (hard coded): %5.6f\n', percent);

position=round(N*percent/100);
sda=sort(xx(:,3));
dc=sda(position)*1;

fprintf('Computing Rho with Cut-off kernel of radius: %12.6f\n', dc);

tic;
rho=zeros(ND,1);
rhoK=zeros(ND,1);
for i=1:ND
    [rho(i),rhoK(i)]=emc2EntropyD(dist(i,:)');
%     rho(i)=minEntropyPDFs(dd,size(S,2));
end

maxd=max(max(dist));

[rho_sorted,ordrho]=sort(rho,'descend');
delta(ordrho(1))=-1.;
nneigh(ordrho(1))=ordrho(1);%ordrho(1);

for ii=2:ND
   delta(ordrho(ii))=maxd;
   for jj=1:ii-1
     if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
        delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
        nneigh(ordrho(ii))=ordrho(jj);
     end
   end
end
% delta(ordrho(1))=max(max(dist));

delta(ordrho(1))=max(delta(:));
disp('Generated file:DECISION GRAPH')
disp('column 1:Density')
disp('column 2:Delta')

fid = fopen('DECISION_GRAPH', 'w');
for i=1:ND
   fprintf(fid, '%6.2f %6.2f\n', rho(i),delta(i));
end

disp('Select a rectangle enclosing cluster centers')
% scrsz = get(0,'ScreenSize');
% figure('Position',[6 72 scrsz(3)/4. scrsz(4)/1.3]);
for i=1:ND
  ind(i)=i;
  gamma(i)=rho(i)*delta(i);
end


cl=ones(1,size(S,1))*(-1);
[~,ordgamma]=sort(gamma);
rrho=mapminmax(rho',0,1);
ddelta=mapminmax(delta,0,1);
ggama=mapminmax(gamma,0,1);
subplot(1,1,1)
tt=plot(rrho(:),ddelta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
title ('Decision Graph','FontSize',15.0)
xlabel ('\rho')
ylabel ('\delta')
% rrho=rho';
% ddelta=delta;
kuser=2;

% [kk]=deltaEntropyD(ddelta',(rrho').^(1/1));
% [kk]=deltaEntropyD(ddelta',(rrho).^(1/1));
[kk1]=emc2EntropyDD1(ddelta');
ddelta0=ddelta;
[ddelta,~]=sort(ddelta);
[rrhoY,kk2]=emc2EntropyDD(rrho');
% index=ordgamma(end-kuser+1:end)';
index=[];
for i=1:kuser
    index=[index;ordgamma(end-kuser+i)];
    index=[index;nneigh(ordgamma(end-kuser+i))];
end
TA2=intersect( find(rrho>rrhoY(kk2)), find(ddelta0>ddelta(kk1)));% zhui su dian  area
TA3=intersect( find(rrho>rrhoY(kk2)), find(ddelta0<ddelta(kk1)));% zhui su dian  area
TA4=intersect( find(rrho<rrhoY(kk2)), find(ddelta0<ddelta(kk1)));% zhui su dian  area
TA2=TA2';
TA3=TA3';
TA4=TA4';

for i=1:ND
    if ddelta(ordrho(i))>ddelta(kk1) && rrho(ordrho(i))>rrhoY(kk2)&& any(ismember(unique(TA3),nneigh(ordrho(i)),'rows'))==0 && any(ismember(unique(TA2),nneigh(ordrho(i)),'rows'))==0 %TA2 
        index=[index;ordrho(i)];
        TA3=[TA3;ordrho(i)];
        TA3=[TA3;nneigh(ordrho(i))];
    end
end

for i=1:ND
    if ddelta0(ordrho(i))>ddelta(kk1) && rrho(ordrho(i))<rrhoY(kk2) && any(ismember(unique(TA4),nneigh(ordrho(i)),'rows'))>0
        index=[index;ordrho(i)];
%         index=[index;nneigh(ordrho(i))];
    end
end


index=unique(index);
kk=size(index,1);
ordrho0=ordrho;

for i=1:kk
    cl(index(i))=i;
    ordrho(find(ordrho==index(i)))=[];
end
icl=index';
    
ordrho=[icl';ordrho];
NCLUST=kk;
Nnneigh=nneigh;
for ii=2:ND
   delta(ordrho(ii))=maxd;
   for jj=1:ii-1
     if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
        delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
        Nnneigh(ordrho(ii))=ordrho(jj);
     end
   end
end
%assignation
for i=1:ND
  if (cl(ordrho(i))==-1)
    cl(ordrho(i))=cl(Nnneigh(ordrho(i)));
  end
end

kin=0;
while (1)
for i=1:kk
    A=find(cl==i);
    if size(A,2)<=rhoK(index(i)) && size(A,2)>0
        for f=1:size(A,2)
            cl(A(f))=cl(nneigh(A(f)));
        end
%         [aa,~]=findnearestcore(index,i,S);
%         for f=1:size(A,2)
%             cl(A(f))=aa;
%         end
    else
        kin=kin+1;
    end    
end
if kin>=kk-1
    break;
else
    kin=0;
end
end

pin=unique(cl);
kk=size(pin,2);
SS=cell(1,kk);
for i=1:kk
    A=find(cl==pin(i));
    S1=zeros(size(A,2),size(S,2));
    for j=1:size(A,2)
        S1(j,:)=S(A(j),:);
    end
    SS(i)={S1};
end

% SS(cellfun(@isempty,SS))=[];

index=index(pin);
SSresult=CrossEmerge(SS,index,kuser);
toc;

ll=[];
re=[];
for j=1:kuser
    tt=zeros(ND,1);
    rr=zeros(ND,1);
    for i=1:ND
        if any(ismember(SSresult{j},S(i,:),'rows'))>0
            rr(i)=j;
            tt(i)=TXT(i,end);
        end
    end
    rr(ismember(rr,0,'rows'))=[];
    tt(ismember(tt,0,'rows'))=[];
    ll=[ll;tt];
    re=[re;rr];
end

cl=ll';
cmap=colormap;
 figure;
for i=1:NCLUST
   ic=int8((i*64.)/(NCLUST*1.)+1);%int8((i*64.)/(NCLUST*1.);
   subplot(2,1,1)
   hold on
   if ic>64
       ic=64;
   end
   plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end
subplot(2,1,2)
disp('Performing 2D nonclassical multidimensional scaling')
Y1 = mdscale(dist, 2, 'criterion','metricsstress');
plot(Y1(:,1),Y1(:,2),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
title ('2D Nonclassical multidimensional scaling','FontSize',15.0)
xlabel ('X')
ylabel ('Y')
for i=1:ND
 A(i,1)=0.;
 A(i,2)=0.;
end
figure;
for i=1:NCLUST
  nn=0;
  ic=int8((i*64.)/(NCLUST*1.)+1);
  for j=1:ND
    if (cl(j)==i)
      nn=nn+1;
      A(nn,1)=Y1(j,1);
      A(nn,2)=Y1(j,2);
    end
  end
  hold on
  if ic>64
       ic=64;
  end
  plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end

faa = fopen('E:\cjy\2017\CLUSTER_ASSIGNATION.txt', 'w');
disp('Generated file:CLUSTER_ASSIGNATION')
disp('column 1:element id')
disp('column 2:cluster assignation without halo control')
disp('column 3:cluster assignation with halo control')
for i=1:ND
   fprintf(faa, '%i %i %i\n',i,cl(i),halo(i));
end
