function f = evaluate_essex_lda(x, class, channels, M);

% Supervised classification of BCI datasets of Univ. of Essex by Julio
% Ortega

Nelect=15; % Number of electrodes
% tcp_con_feat=importdata('datos.mat'); % Training data are obtained from datos.mat

% testchannels=channels
V=size(x,2);
Velect(1:Nelect)=0;
K=M+V;
Nclas=max(class);

for i = 1 : V
    if (channels(i) > 0)
        Velect(channels(i))=1;
    end
end

CV0 = cvpartition(class,'k',10);
L(1:CV0.NumTestSets)=0;
for s=1 : CV0.NumTestSets
    
      trainlabels=class(CV0.training(s));
      testlabels=class(CV0.test(s));
      trainset=x(CV0.training(s),:);
      testset=x(CV0.test(s),:); 
 
      linclass=ClassificationDiscriminant.fit(trainset,trainlabels);
      cvmodelo=crossval(linclass);
      L(s)=kfoldLoss(cvmodelo);
end

C2=mean(L);

linclass=ClassificationDiscriminant.fit(x,class);
class_estimated=resubPredict(linclass);

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

C1=1-Kappa;


f = [];

% Selection of two objectives
%f(1)=-features_selected; % Number of selected features (minimization)....
f(1)=C1;
f(2)=C2;

end

