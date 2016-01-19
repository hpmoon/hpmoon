function parsave (chromosome,rep,V, mI)

%save (['solution',num2str(rep),'.txt'], 'chromosome', '-ASCII');
save (['chromosome_P-',num2str(rep),'Ite-',num2str(mI),'.mat'], 'chromosome');
%save (['chromosomeFinal_P-',num2str(rep),'Ite-',num2str(mI),'.mat'], 'pCh');
%save ('chromosomefull.mat', 'chromosome', '-append');
%save (['distance',num2str(rep),'.mat'], 'avgdis');

%result=chromosome(:,(V+1):size(chromosome,2));
%save (['4-pareto_P-',num2str(rep),'Ite-',num2str(mI),'.mat'], 'result');

%for i = 1 : size(chromosome,1)
%    xyz=find(chromosome(i,1:V));
%    save (['5-features_P-',num2str(rep),'Ite-',num2str(mI),'.txt'], 'xyz', '-append', '-ASCII');
%end