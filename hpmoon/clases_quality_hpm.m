function [f,err,ind] = clases_quality_hpm(M, V);
global tcp_con_feat;
global labels;
global selec;
% Function to evaluate the performance for the given input vector
% x of selected features. x is an array of selected decision variables (0's and 1's).
% M is the number of objective functions and
% V is the number of decision variables. 
% By Andres Ortiz and Julio Ortega

g = sprintf('Number of selected features: ');
% Obtain the number of selected features
FF = input(g); % Number of features

x=zeros(1,V);
for i = 1 : FF
    g = sprintf('Feature:%d: ',i);
    % Obtain the indices of the selected features
    kk = input(g); % Number of features;
    x(kk)=round(1);
end

g = sprintf('Number of clases: ');
% Obtain the number of selected features
Nclas = input(g); % Number of clases

% Data file is loaded in matrix tcp_con_feat 
tcp_con_feat=importdata('datos.mat');
classes=importdata('clases.mat'); 
data_test=importdata('datos_test.mat');
classes_test=importdata('clases_test.mat');
K=M+V; % Objectives plus decision variables (= features)

selec=zeros(1,K);
%fprintf('Size of selec: %d...\n\n',size(selec,2));
fvect=zeros(1,size(tcp_con_feat,2));
for i = 1 : size(tcp_con_feat,2)
     if x(i)==1
         fvect(i)=1; % used features of tcp_con_feat are the components set to one in fvect
     end
  end

tcp_con_feat2(:,find(find(fvect)))=tcp_con_feat(:,find(fvect));
data_test2(:,find(find(fvect)))=data_test(:,find(fvect));
%fprintf('Selected features: %d...\n\n',size(tcp_con_feat2,2));
%

tcp_con_feat_n=som_normalize(tcp_con_feat2,'var');    % Vector of features is normalized
data_test_n=som_normalize(data_test2,'var');

% Sustituido: sMap = som_randinit(tcp_con_feat_n,'msize',[10 10],'rect'); % SOM map definition and random intitialization
sMap=som_lininit(tcp_con_feat_n,'hexa','msize',[10 10]);

subplot(2,1,1);
% SOM training with the selected features
% sMap = som_make(tcp_con_feat_n,'msize',[10 10],'ep','training','long');
% Sustituido: sMap = som_seqtrain(sMap,tcp_con_feat_n,'radius',[5 1],'trainlen',60,'tracking',3);
sMap = som_batchtrain(sMap,tcp_con_feat_n,'radius_ini',5,'radius_fin',1,'trainlen',60,'tracking',3);

som_show(sMap);
hold off;
%&&&& subplot(1,1,1)

% Quality and topology characteristics evaluated by som_quality
[quality,topology]=som_quality(sMap,tcp_con_feat_n);

% ------------------------------------------------------------------
% Evaluation of IDB (Index of Davies-Bouldin) from the sMap (without
% k-means clustering) from Voronoi sets
[BMU,QERR]=som_bmus(sMap,tcp_con_feat_n);
[hits] = som_hits(sMap,tcp_con_feat_n);

S(1:size(hits,1))=0;
centro(1:size(hits,1),1:size(tcp_con_feat_n,2))=0;
B(1:size(hits,1),1:size(hits,1))=0;
IDB=0;
kk=0;


% ---------------- Classes assignment to units
uni_class(1:max(BMU),1:Nclas)=0;
for i=1:size(tcp_con_feat_n,1)
    uni_class(BMU(i),classes(i))=uni_class(BMU(i),classes(i))+1;
end
%  sum(sum(uni_class));

[class_weight class_estimated]=max(uni_class,[],2);
for i=1:size(uni_class,2)
    bi_class(:,i)=uni_class(:,size(uni_class,2)-i+1);
end    
[biclass_weight biclass_estimated]=max(bi_class,[],2);

busca_class(1:size(uni_class,1),1)=0;
for i=1:size(uni_class,1)
    if (class_weight(i)==0)
        class_estimated(i)=0;
        busca_class(i)=1;
    else
        if (class_estimated(i)~=(size(uni_class,2)-biclass_estimated(i)+1))
           busca_class(i)=1;
        end
    end
end

units_assigned=find(class_estimated);
units_toassign=find(busca_class);

for i=1:size(units_toassign,1)
    for j=1:size(units_assigned)
        distance=norm(sMap.codebook(units_toassign(i),:)-sMap.codebook(units_assigned(j),:));
        if (j==1)
            MIN=distance;
            class_estimated(units_toassign(i))=class_estimated(units_assigned(j));
        else
            if (distance<MIN)
                MIN=distance;
                class_estimated(units_toassign(i))=class_estimated(units_assigned(j));
            end
        end 
    end  
end
% ----------------------- End of classes assignment to units

for i=1:size(BMU,1)
        for j=1:size(tcp_con_feat_n,2)
            centro(BMU(i),j)=centro(BMU(i),j)+(tcp_con_feat_n(i,j)/hits(BMU(i)));
        end
end

%centroids of clusters 
subplot(1,1,1)
hold on;
for i=1:size(tcp_con_feat_n,1)
    if (classes(i)==1) 
        plot(tcp_con_feat_n(i,1),tcp_con_feat_n(i,2),'ok');
    end    
    if (classes(i)==2) 
        plot(tcp_con_feat_n(i,1),tcp_con_feat_n(i,2),'og');
    end    
    if (classes(i)==3) 
        plot(tcp_con_feat_n(i,1),tcp_con_feat_n(i,2),'ob');
    end    
    if (classes(i)==4) 
        plot(tcp_con_feat_n(i,1),tcp_con_feat_n(i,2),'*k');
    end   
end
plot(centro(:,1),centro(:,2),'+r');

for i=1:size(sMap.codebook,1) % units of the map
            for ii=1:size(tcp_con_feat_n,1) % input patterns 
                if (BMU(ii)==i)
                    for jj=1:size(centro,2)
                        vector(jj)=tcp_con_feat_n(ii,jj)-centro(i,jj);
                    end    
                    S(i)=S(i)+(norm(vector)/hits(i));
                end
            end
end

for i=1:size(sMap.codebook,1)
            for ii=1:size(sMap.codebook,1) 
                    if ((hits(i)~=0)&&(hits(ii)~=0))  % OJO i se ha cambiado a ii 
                        B(i,ii)=norm(centro(i)-centro(ii));
                    end
            end
end

for i=1:size(sMap.codebook,1)
            R(i)=0;
            for ii=1:size(sMap.codebook,1) 
                    if (i~=ii)
                        if (B(i,ii)~=0)
                            SS=(S(i)+S(ii))/B(i,ii);
                            if (R(i)<SS)
                                R(i)=SS;
                            end
                        end    
                    end   
            end
end

for i=1:size(sMap.codebook,1)
            IDB=IDB+R(i);
            if (R(i)~=0)
                kk=kk+1;
            end
end
IDB=IDB/kk;
% ------------------------------------------- IDB has been calculated

% ------------------ Classification performance (labeled patterns) 
correct=0;
incorrect=0;
TN(1:Nclas)=0;
TP(1:Nclas)=0;
FN(1:Nclas)=0;
FP(1:Nclas)=0;
ACCUR(1:Nclas)=0;
SENS(1:Nclas)=0;
SPEC(1:Nclas)=0;
JACCARD(1:Nclas)=0;
MatrixKappa(1:Nclas,1:Nclas)=0;
NKappa(1:Nclas)=0;
MKappa(1:Nclas)=0;
ProbClases(1:Nclas)=0;
ProbClasif(1:Nclas)=0;
%Total=size(tcp_con_feat_n,1);
Total=size(data_test_n,1);
Po=0;
Pc=0;

[BMU_class,QERR_class]=som_bmus(sMap,data_test_n);
% Kappa coefficient computation
% for i=1:size(tcp_con_feat_n,1)
for i=1:size(data_test_n,1)
    if (class_estimated(BMU_class(i))==classes_test(i))
        correct=correct+1;
    else
        incorrect=incorrect+1;
    end
    MatrixKappa(class_estimated(BMU_class(i)),classes_test(i))= MatrixKappa(class_estimated(BMU_class(i)),classes_test(i))+1;
    MKappa(class_estimated(BMU_class(i)))= MKappa(class_estimated(BMU_class(i)))+1;
    NKappa(classes_test(i))= NKappa(classes_test(i))+1;
end

for i=1:Nclas
    ProbClases(i)=NKappa(i)/Total;
    ProbClasif(i)=MKappa(i)/Total;
    Po = Po + MatrixKappa(i,i)/Total;
    Pc = Pc + ProbClases(i)*ProbClasif(i);
end
Kappa= (Po-Pc)/(1-Pc)
% End of Kappa coefficient computation    

% for i=1:size(tcp_con_feat_n,1)
for i=1:size(data_test_n,1)
        if (class_estimated(BMU_class(i))==classes_test(i))
            for j=1:Nclas
                if (j==classes_test(i)) 
                    TP(j)=TP(j)+1;
                else
                    TN(j)=TN(j)+1;
                end
            end
        else
            FP(class_estimated(BMU_class(i)))=FP(class_estimated(BMU_class(i)))+1;
            for j=1:Nclas
                if (j==classes_test(i))
                    FN(j)=FN(j)+1;
                end
            end
        end
end       

for i=1:Nclas
    ACCUR(i)=(TP(i)+TN(i))/(TP(i)+TN(i)+FP(i)+FN(i));
    SENS(i)=TP(i)/(TP(i)+FN(i));
    SPEC(i)=TN(i)/(TN(i)+FP(i));
    JACCARD(i)=TP(i)/(TP(i)+TN(i)+FP(i));
end

truefalse = [TP; TN; FP; FN];
JACC=[ACCUR; SENS; SPEC; JACCARD];

% New index calculation (map units)
[hits] = som_hits(sMap,tcp_con_feat_n);
indices=find(hits);
MAX=0;
dist=zeros(size(sMap.codebook,1),size(sMap.codebook,1));

for ii = 1 : size(indices)
    for jj = 1 : size(indices)
            if (indices(ii)~=indices(jj))
                dist(indices(ii),indices(jj))=norm(sMap.codebook(indices(ii),:)-sMap.codebook(indices(jj),:));
                if (dist(indices(ii),indices(jj))>MAX) 
                    MAX=dist(indices(ii),indices(jj));
                end
            end   
    end
end
MAX=MAX+1;
for ii = 1 : size(sMap.codebook,1)
    dist(ii,ii)=MAX;
    if (hits(ii)==0)
        for jj = 1 : size(sMap.codebook,1)
            dist(ii,jj)=MAX;
            dist(jj,ii)=MAX;
        end
    end
end


DD(1,1)=indices(1);
DD(1,2)=0;
dist(:,DD(1,1))=MAX;

for ii=2 : size(indices)
    MIN = min(dist(DD(ii-1,1),:));
    test = 0;
    for jj = 1 : size(indices)
        if (test==0)
            if (dist(DD(ii-1,1),indices(jj))==MIN)
                DD(ii,1)=indices(jj);
                DD(ii,2)=MIN;
                test=1;
            end
        end
    end
    for kk = 1 : size(sMap.codebook,1)
%        dist(DD(ii,1),kk)=MAX;
        dist(kk,DD(ii,1))=MAX;
    end
end
order_dist=sort(DD(:,2));
AVG=mean(DD(:,2))
C1=0;
N1=0;
C2=0;
N2=0;

for ii = 1 : size(indices)
    if (DD(ii,2)<AVG)
        C1=C1+DD(ii,2);
        N1=N1+1;
    else
        C2=C2+(DD(ii,2));
        N2=N2+1;
    end
end
C1=C1/N1
C2=-C2/N2
cut=N1
% End of new index calculation (map units)     
       

% New index calculation (data)
MAX_d=0;
N_d=size(tcp_con_feat_n,1);
for ii = 1 : N_d
    for jj = 1 : N_d
        if (ii == jj) 
            dist_d(ii,jj)=0;
        else
            dist_d(ii,jj)=norm(tcp_con_feat_n(ii,:)-tcp_con_feat_n(jj,:));
            if (dist_d(ii,jj)>MAX_d) 
                MAX_d=dist_d(ii,jj);
            end
        end
    end
end
MAX_d=MAX_d+1;
for ii = 1 : N_d
    dist_d(ii,ii)=MAX_d;
end
DD_d(1,1)=1;
DD_d(1,2)=0;

for ii=2 : N_d
    MIN_d = min(dist_d(DD_d(ii-1,1),:));
    test = 0;
    for jj = 1 : N_d
        if (test==0)
            if (dist_d(DD_d(ii-1,1),jj)==MIN_d)
                DD_d(ii,1)=jj;
                DD_d(ii,2)=MIN_d;
                test=1;
            end
        end
    end
    for kk = 1 : N_d
        dist_d(kk,DD_d(ii,1))=MAX_d;
    end
end
length_d=sort(DD_d(:,2));
AVG_d=mean(DD_d(:,2))
C1_d=0;
N1_d=0;
C2_d=0;
N2_d=0;

for ii = 1 : N_d
    if (DD_d(ii,2)<AVG_d)
        C1_d=C1_d+DD_d(ii,2);
        N1_d=N1_d+1;
    else
        C2_d=C2_d+(DD_d(ii,2));
        N2_d=N2_d+1;
    end
end
C1_d=C1_d/N1_d
C2_d=-C2_d/N2_d
cut_d=N1_d
% End of new index calculation (data)


%hold on;
%plot(sMap.codebook(:,1),sMap.codebook(:,2),'d');
%plot(tcp_con_feat_n(:,1),tcp_con_feat_n(:,2),'+r');
%hold off
%title('Data and map after training');

% Clustering by k-means (kmeans_clusters)
[c, p, err, ind] = kmeans_clusters(sMap,20); % find clusterings (explores up to 20 clusters)
[dummy,i] = min(ind); % select the one with smallest index
%som_show(sMap,'color',{p{i},sprintf('%d clusters',i)}); % visualize it
%colormap(jet(i)), som_recolorbar % change colormap

selected_features=find(fvect);
resultados=[selected_features quality topology IDB C1 C2 C1_d C2_d correct incorrect Kappa]
quality_measures=[quality topology IDB C1 C2 C1_d C2_d correct incorrect Kappa]
save measure.txt resultados -ASCII -append;
save paths.txt order_dist -ASCII -append;
save paths.txt length_d -ASCII -append;
save qm.txt quality_measures -ASCII -append;
save bmu.txt uni_class -ASCII;
save bmu.mat uni_class;
save busca_class.mat busca_class;
save class_weight.mat class_weight;
save class_estimated.mat class_estimated;
save TPTNFPFN.mat truefalse;
save JACCARD.mat JACC;

