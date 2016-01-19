% function f = evaluate_objective_hpm(x, M, V,tcp_con_feat,labels)
function f = evaluate_objective_hpm(x, M, V);
global tcp_con_feat;
global labels;
global selec;
% Function to evaluate the objective functions for the given input vector
% x. x is an array of decision variables and f(1), f(2), etc are the
% objective functions. The algorithm always minimizes the objective
% function hence if you would like to maximize the function then multiply
% the function by negative one. M is the numebr of objective functions and
% V is the number of decision variables.
% Training patterns in "datos.mat"
%
% This functions is basically written by the user who defines his/her own
% objective function. Make sure that the M and V matches your initial user
% input. Make sure that the 
%
% An example objective function is given below. It has two six decision
% variables are two objective functions.

% f = [];
% %% Objective function one
% % Decision variables are used to form the objective function.
% f(1) = 1 - exp(-4*x(1))*(sin(6*pi*x(1)))^6;
% sum = 0;
% for i = 2 : 6
%     sum = sum + x(i)/4;
% end
% %% Intermediate function
% g_x = 1 + 9*(sum)^(0.25);
% 
% %% Objective function two
% f(2) = g_x*(1 - ((f(1))/(g_x))^2);
% Adapted for Multiobjective feature selection with SOM by Julio Ortega and
% Andrés Ortiz


tcp_con_feat=importdata('datos.mat'); % Training data are obtained from datos.mat

K=M+V;

selec=zeros(1,K);
%fprintf('Size of selec: %d...\n\n',size(selec,2));
fvect=zeros(1,size(tcp_con_feat,2));  % the component of fvect are initialized to 0
for i = 1 : size(tcp_con_feat,2)
     if x(i)==1
         fvect(i)=1;
     end
  end


tcp_con_feat2(:,find(find(fvect)))=tcp_con_feat(:,find(fvect));
%fprintf('Elementos seleccionados: %d...\n\n',size(tcp_con_feat2,2));
%
% Number of features selected for clustering
features_selected=length(find(fvect)); % One of the objectives could be the number of selected features
features=find(fvect);

% Training of the SOM by using the selected features
tcp_con_feat_n=som_normalize(tcp_con_feat2,'var');    % Normalization of vectors

sMap = som_randinit(tcp_con_feat_n,'msize',[10 10]); % 10x10 units (can be changed)

% sMap = som_make(tcp_con_feat_n,'msize',[10 10],'ep','training','long');
%Alternative training

%sMap = som_seqtrain(sMap,tcp_con_feat_n,'radius',[5 1],'trainlen',60);
sMap = som_batchtrain(sMap,tcp_con_feat_n,'radius',[5 1],'trainlen',60);

% Below several possible objective functions are computed (uncomment those
% used as objectives in a give run)

% Evaluation of the trained SOM map by som_quality
[quality,topology]=som_quality(sMap,tcp_con_feat_n);

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
        dist(kk,DD(ii,1))=MAX;
    end
end
order_dist=sort(DD(:,2));
AVG=mean(DD(:,2));
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
C1=C1/N1;
C2=-C2/N2;
cut=N1;
% End of new index calculation (map units)   

 f = [];

% Selection of two objectives
%f(1)=-features_selected; % Number of selected features (minimization)....
f(1)=C1;
f(2)=C2;

for i=3:M % Two set the rest of objectives if it would indicated more than two objectives
    f(i)=10*rand(1);
    fprintf('f(%d)=%d\n',i,f(i));
end

