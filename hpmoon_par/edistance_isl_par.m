function [chromosome] = edistance_isl_par (chromosome, pop, V)
chromosome = sortrows(chromosome,[V+1,V+2]);
distance = 0;
for reuc = 1 : pop-1
    dist1 = [chromosome(reuc,V+1),chromosome(reuc,V+2);chromosome(reuc+1,V+1),chromosome(reuc+1,V+2)];
    distance = distance + pdist(dist1);      
end 
avgdist = distance/pop;

for reuc = 2 : pop-1
    dist1 = [chromosome(reuc-1,V+1),chromosome(reuc-1,V+2);chromosome(reuc,V+1),chromosome(reuc,V+2)];
    dist2 = [chromosome(reuc,V+1),chromosome(reuc,V+2);chromosome(reuc+1,V+1),chromosome(reuc+1,V+2)];
    if (pdist(dist1)<avgdist) && (pdist(dist2)<avgdist)
        chromosome (reuc,V+5) = 1;
    else 
        chromosome (reuc,V+5) = 0;
    end
end 
