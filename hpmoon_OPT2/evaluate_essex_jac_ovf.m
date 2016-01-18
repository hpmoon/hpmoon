function f = evaluate_essex_jac_ovf(x, class, M);
%global tcp_con_feat;  %removed in essex version
%global labels;  %removed in essex version
%global selec;  %removed in essex version
% Function to evaluate the objective functions for the given input vector
% x. x is an array of training patterns (with more information)
% and f(1), f(2),.. are the objective functions. 
% The algorithm always minimizes the objective
% function hence if you would like to maximize the function then multiply
% the function by negative one. M is the numebr of objective functions and
% V is the number of decision variables.
%
% This functions is basically written by the user who defines his/her own
% objective function. Make sure that the M and V matches your initial user
% input.
%

Nelect=15; % Number of electrodes
% tcp_con_feat=importdata('datos.mat'); % Training data are obtained from datos.mat

% testchannels=channels
% 
V=round(size(x,2)/Nelect);
% V=size(x,2);
% Velect(1:Nelect)=0;
K=M+V;
Nclas=max(class);

% for i = 1 : V
%    if (channels(i) > 0)
%        Velect(channels(i))=1;
%    end
% end
% C1=sum(Velect)
%%%%%%%%C1=V; (substitute for Koverf)

%%%%% Average of Kappa index


CV0 = cvpartition(class,'k',10);
% L(1:CV0.NumTestSets)=0;
for s=1 : CV0.NumTestSets
    
    trainlabels=class(CV0.training(s));
    testlabels=class(CV0.test(s));
      
    votes=[];  
    for i = 1 : V  
      
      trainset=x(CV0.training(s),(i-1)*Nelect+1:i*Nelect);
      testset=x(CV0.test(s),(i-1)*Nelect+1:i*Nelect); 
 
      linclass=ClassificationDiscriminant.fit(trainset,trainlabels);
      class_predicted=predict(linclass,testset);
      
      votes=[votes class_predicted];
    end
    vote=double(votes);
    class_estimated=mode(vote,2);

    MatrixKappa(1:Nclas,1:Nclas)=0;
    NKappa(1:Nclas)=0;
    MKappa(1:Nclas)=0;
    ProbClases(1:Nclas)=0;
    ProbClasif(1:Nclas)=0;
    Total=size(testset,1);
    Po=0;
    Pc=0;

     % Kappa coefficient computation

    for i = 1 : Total
        MatrixKappa(class_estimated(i),testlabels(i))= MatrixKappa(class_estimated(i),testlabels(i))+1;
        MKappa(class_estimated(i))= MKappa(class_estimated(i))+1;
        NKappa(testlabels(i))= NKappa(testlabels(i))+1;
    end

     for i=1:Nclas
        ProbClases(i)=NKappa(i)/Total;
        ProbClasif(i)=MKappa(i)/Total;
        Po = Po + MatrixKappa(i,i)/Total;
        Pc = Pc + ProbClases(i)*ProbClasif(i);
     end
     Koverf(s)= (Po-Pc)/(1-Pc);
end

C1=1-mean(Koverf);


%%%%%% Kappa indez as one of the cost function
votes=[];

for i = 1 : V
    LDA = x(:,(i-1)*Nelect+1:i*Nelect);
    linclass=ClassificationDiscriminant.fit(LDA,class);
    % ConfusionMatrix=confusionmat(linclass.Y,resubPredict(linclass));
    class_predicted=resubPredict(linclass);
    votes=[votes class_predicted];
end

vote=double(votes);

class_estimated=mode(vote,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CV0 = cvpartition(class,'k',10);
% L(1:CV0.NumTestSets)=0;
% for s=1 : CV0.NumTestSets   
%     trainlabels=class(CV0.training(s));
%     testlabels=class(CV0.test(s));
%     trainset=x(CV0.training(s),:);
%     testset=x(CV0.test(s),:); 
%     linclass=ClassificationDiscriminant.fit(trainset,trainlabels);
%     cvmodelo=crossval(linclass);
%     L(s)=kfoldLoss(cvmodelo);
      
      
% Continue here with the training process
% save clases.mat trainlabels;
% save clases_test.mat testlabels;
% save datos.mat trainset;
% save datos_test.mat testset;

% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------ Classification performance (labeled patterns) 
correct=0;
incorrect=0;
% TN(1:Nclas)=0;
% TP(1:Nclas)=0;
% FN(1:Nclas)=0;
% FP(1:Nclas)=0;
% ACCUR(1:Nclas)=0;
% SENS(1:Nclas)=0;
% SPEC(1:Nclas)=0;
% JACCARD(1:Nclas)=0;
MatrixKappa(1:Nclas,1:Nclas)=0;
NKappa(1:Nclas)=0;
MKappa(1:Nclas)=0;
ProbClases(1:Nclas)=0;
ProbClasif(1:Nclas)=0;
Total=size(x,1);
Po=0;
Pc=0;

% Kappa coefficient computation

for i=1:Total
    if (class_estimated(i)==class(i))
        correct=correct+1;
    else
        incorrect=incorrect+1;
    end
    MatrixKappa(class_estimated(i),class(i))= MatrixKappa(class_estimated(i),class(i))+1;
    MKappa(class_estimated(i))= MKappa(class_estimated(i))+1;
    NKappa(class(i))= NKappa(class(i))+1;
end

for i=1:Nclas
    ProbClases(i)=NKappa(i)/Total;
    ProbClasif(i)=MKappa(i)/Total;
    Po = Po + MatrixKappa(i,i)/Total;
    Pc = Pc + ProbClases(i)*ProbClasif(i);
end
Kappa= (Po-Pc)/(1-Pc);
% End of Kappa coefficient computation    

% C2=mean(L);
C2=1-Kappa;

f = [];

% Selection of two objectives
%f(1)=-features_selected; % Number of selected features (minimization)....
f(1)=C1;
f(2)=C2;
%f(2)=quality; % Quality obtained with som_quality
%f(2)=idb_k; % IDB obtained with kmeans_clusters
%f(2)=err_clust; % Error of clustering obtained with kmeans_clusters
%f(2)=IDB; % IDB obtained from sMap (without kmeans_clusters)
for i=3:M % Two set the rest of objectives if it would indicated more than two objectives
    f(i)=10*rand(1);
    fprintf('f(%d)=%d\n',i,f(i));
end

