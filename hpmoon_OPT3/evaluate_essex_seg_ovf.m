function f = evaluate_essex_seg_ovf(x, class, NFeat, Nelect, M);

% Defined for BCI datasets of Univ. of Essex by Julio Ortega

Nsegments=size(NFeat); % Number of Segments
 
% V=round(size(x,2));

Nclas=max(class);

% C1=sum(Nelect); (substituted by overfitting control cost function)

%%%%% Start overfitting control cost

CV0 = cvpartition(class,'k',10);

for s=1 : CV0.NumTestSets
    
    trainlabels=class(CV0.training(s));
    testlabels=class(CV0.test(s));

    votes=[];
    p=1;
    for i = 1 : Nsegments

      trainset= x(CV0.training(s),p:p+(NFeat(i)-1));
      testset= x(CV0.test(s),p:p+(NFeat(i)-1));
      p=p+NFeat(i);
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
%%%%%% End control of overfitting cost function


votes=[];
p=1;
for i = 1 : Nsegments
    LDA = x(:,p:p+(NFeat(i)-1));
    p=p+NFeat(i);
    linclass=ClassificationDiscriminant.fit(LDA,class);
    % class_predicted=resubPredict(linclass);
    class_predicted=predict(linclass,LDA);
    votes=[votes class_predicted];
end

vote=double(votes);
class_estimated=mode(vote,2);

% ------------------ Classification performance (labeled patterns) 
correct=0;
incorrect=0;
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

f(1)=C1;
f(2)=C2;

for i=3:M % Two set the rest of objectives if it would indicated more than two objectives
    f(i)=10*rand(1);
    fprintf('f(%d)=%d\n',i,f(i));
end

end
