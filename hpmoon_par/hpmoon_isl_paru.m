function hpmoon_isl_paru(pop,gen,iterat,islands)

% based on function hpmoon(pop,gen)
% is a multi-objective optimization function where the input arguments are 
% pop - Population size
% gen - Total number of generations
% 
% This functions is based on evolutionary algorithm for finding the optimal
% solution for multiple objective i.e. Pareto front for the objectives. 
% Initially enter only the population size and the stoping criteria or
% the total number of generations after which the algorithm will
% automatically stopped. 
%
% You will be asked to enter the number of objective functions, and the number
% of decision variables.
% Also you will have to define your own objective funciton by editing the
% evaluate_objective() function. A sample objective function is described
% in evaluate_objective.m. 
%
% Original algorithm NSGA-II was developed by researchers in Kanpur Genetic
% Algorithm Labarotary and kindly visit their website for more information
% http://www.iitk.ac.in/kangal/

%  Copyright (c) 2009, Aravind Seshadri
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%      
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%  POSSIBILITY OF SUCH DAMAGE.
%
% This software has been adapted to the problem of parallel multiobjective feature
% selection by Andres Ortiz, Julio Ortega and Dragi Kimovski

% Number of Arguments
% Check for the number of arguments. The four input arguments are necessary
% to run this function.
if nargin < 4
    error('Please enter the population size, number of generations, number of main iterations and number of islands');
end
% Both the input arguments need to of integer data type
if isnumeric(pop) == 0 || isnumeric(gen) == 0 || isnumeric(islands) == 0 || isnumeric(iterat) == 0
    error('The three input arguments should be integer datatype');
end
% Minimum population size has to be 20 individuals and minimum number of
% generations 2
if pop < 20
    error('Minimum population for running this function is 20');
end
if gen < 2
    error('Minimum number of generations is 2');
end
% Make sure pop, gen, islands are integers
pop = round(pop);
gen = round(gen);
islands = round(islands);
iterat = round (iterat);

% The objective function description contains information about the
% objective functions. M is the dimension of the objective space, V is the
% dimension of decision variable space.
% The individuals in the population are vectors of V binary components, 
% so min_range=0 and max_range=1
% User has to define the objective functions using the decision variables. 
% Make sure to edit the function 'evaluate_objective_hpm' to suit your needs.

g = sprintf('Number of objectives: ');
% Obtain the number of objective function
M = input(g); % Number of Objectives
if M < 2
    error('This is a multi-objective optimization function hence the minimum number of objectives is two');
end
g = sprintf('\nNumber of decision variables (features): ');
% Obtain the number of decision variables
V = input(g); % Number of componets of the individuals in the population
%g = sprintf('\nNumber of modules or islands (it divides %i): ',V);
% Obtain the number of features considered by an island
%islands = input(g);
%chunk = floor(islands*rand(1));
MAXfeat=round(V/islands);

min_range=round(0); % Individuals with binary components (min_range=0)
max_range=round(1); % Individuals with bynary components (max_range=0)

K=M+V; % Components of the individuals (possible features) plus objetives 


% MATRIX Initialization
for j = 1 : V
       kkk=round(rand(1));
        for i = 1 : pop
            fn1(i,j)=kkk;
        end
end
ef=fn1;
rChromosome = [];

%% START OF THE MAIN ITERATIONS
ti=clock;
for mI=1:iterat   
Mchromosome = [];
Mchromosome = rChromosome;
if mI~=1
    ef=rChromosome(:,1:V);
end    
%START OF THE PARALLEL CODE
% This paralelisation is perfomed by using the Island parallel paradigm 
parfor chunk = 0 : islands-1
%% Initialize the population
% Population is initialized with random values which are within the
% specified range. Each chromosome consists of the decision variables. Also
% the value of the objective functions, rank and crowding distance
% information is also added to the chromosome vector but only the elements
% of the vector which has the decision variables are operated upon to
% perform the genetic operations like corssover and mutation.

if mI == 1
[chromosome] = initialize_variables_hpm_isl_par(pop, M, V, min_range, max_range, MAXfeat, chunk, ef, mI); 
%% Sort the initialized population
% Sort the population using non-domination-sort. This returns two columns
% for each individual which are the rank and the crowding distance
% corresponding to their position in the front they belong. At this stage
% the rank and the crowding distance for each chromosome is added to the
% chromosome vector for easy of computation.
chromosome = non_domination_sort_mod(chromosome, M, V);
else
chromosome = Mchromosome;
end;

%% Start the evolution process
% The following are performed in each generation
% * Select the parents which are fit for reproduction
% * Perfrom crossover and Mutation operator on the selected parents
% * Perform Selection from the parents and the offsprings
% * Replace the unfit individuals with the fit individuals to maintain a
%   constant population size.

for i = 1 : gen
    
    fprintf('GENERACION: %d \n',i);
    % Select the parents
    % Parents are selected for reproduction to generate offspring. The
    % original NSGA-II uses a binary tournament selection based on the
    % crowded-comparision operator. The arguments are 
    % pool - size of the mating pool. It is common to have this to be half the
    %        population size.
    % tour - Tournament size. Original NSGA-II uses a binary tournament
    %        selection, but to see the effect of tournament size this is kept
    %        arbitary, to be choosen by the user.
    pool = round(pop/2);
    tour = 2;
    % Selection process
    % A binary tournament selection is employed in NSGA-II. In a binary
    % tournament selection process two individuals are selected at random
    % and their fitness is compared. The individual with better fitness is
    % selcted as a parent. Tournament selection is carried out until the
    % pool size is filled. Basically a pool size is the number of parents
    % to be selected. The input arguments to the function
    % tournament_selection are chromosome, pool, tour. The function uses
    % only the information from last two elements in the chromosome vector.
    % The last element has the crowding distance information while the
    % penultimate element has the rank information. Selection is based on
    % rank and if individuals with same rank are encountered, crowding
    % distance is compared. A lower rank and higher crowding distance is
    % the selection criteria.
    
    % fprintf('Tournament %d \n',i);
    parent_chromosome = tournament_selection(chromosome, pool, tour);
    %parsaveTEST (parent_chromosome, chunk);
    % Perfrom crossover and Mutation operator
    % The original NSGA-II algorithm uses Simulated Binary Crossover (SBX) and
    % Polynomial  mutation. Crossover probability pc = 0.9 and mutation
    % probability is pm = 1/n, where n is the number of decision variables.
    % Both real-coded GA and binary-coded GA are implemented in the original
    % algorithm, while in this program only the real-coded GA is considered.
    % The distribution indeices for crossover and mutation operators as mu = 20
    % and mum = 20 respectively.
    mu = 20;
    mum = 20;
   
    % fprintf('Genetic Operator %d \n',i);
    offspring_chromosome = ...
        genetic_operator_hpm_isl_par(parent_chromosome, ...
        M, V, mu, mum, min_range, max_range, MAXfeat, chunk, ef, mI);
    % Intermediate population
    % Intermediate population is the combined population of parents and
    % offsprings of the current generation. The population size is two
    % times the initial population.
    
    % fprintf('Size %d \n',i);
    [main_pop,temp] = size(chromosome);
    [offspring_pop,temp] = size(offspring_chromosome);
    % temp is a dummy variable.
    %clear temp
    % intermediate_chromosome is a concatenation of current population and
    % the offspring population.
    
    % fprintf('Intermediate %d \n',i);
    intermediate_chromosome(1:main_pop,:) = chromosome;
    intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,1 : M+V) = ...
        offspring_chromosome;

    % Non-domination-sort of intermediate population
    % The intermediate population is sorted again based on non-domination sort
    % before the replacement operator is performed on the intermediate
    % population.
    
    % fprintf('Non domination sort %d \n',i);
    intermediate_chromosome = ...
        non_domination_sort_mod(intermediate_chromosome, M, V);
    % Perform Selection
    % Once the intermediate population is sorted only the best solution is
    % selected based on it rank and crowding distance. Each front is filled in
    % ascending order until the addition of population size is reached. The
    % last front is included in the population based on the individuals with
    % least crowding distance
    
    % fprintf('Replace chromosome %d \n',i);
    chromosome = replace_chromosome(intermediate_chromosome, M, V, pop);
    if ~mod(i,100)
        clc
        fprintf('%d generations completed \n',i);
    end
end




%% Result
% Save the result in ASCII and MAT format.
outputch {chunk+1} = chromosome; 
%parsave(chromosome,chunk,V,mI);
end

FTM = outputch; 
%% PREPARE FOR THE NEXT ITERATION

%% Concatenate the chromosomes from all processors and plot the joint matrix
clear chromosome;
clear thv;
FMT = [];FM = []; maxE = 0; maxT = 0; maxR = 0; msize = 0;
paretoCT = [];
plotChromosome = [];
rChromosome = [];
alllevels = [];
firstlevel = [];
addelements = [];
plotC = []; 
thv=etime(clock,ti);
for cnum = 1: islands
FMT = vertcat (FMT,outputch{1,cnum});
end
FM = FMT (:,1:V+2);
FM = unique(FM, 'rows');
FM = non_domination_sort_mod(FM, M, V);
msize = size (FM,1);
%--------------------------------------------------------------------------%
%% If this is the last iteration save the best sollution 
if mI == iterat
    br9 = 0;
    plotC = FM;
for i = 1 : msize
            if (plotC (i,V+3) == 1)
            plotChromosome (i-br9,:) = plotC (i,:);
            else
            br9=br9+1;    
            end
end
for w=1 : M
paretoC(:,w)= plotChromosome(:,V+w);
end
limit = [12,0];
[hv,hv1] = hypervolume (paretoC,limit);
fprintf('The hypervolume value for the last pareto front is %e /// %e \n',hv,hv1);
Volume1 (mI,1) = hv;
Volume2 (mI,1) = hv1;
Volume2 (mI,2) = thv;
save ('hypervolume1.mat','Volume1');
save ('hypervolume2.mat','Volume2');
if hv1 > archiveHV 
archive = [];
archive = plotChromosome;
archiveHV = hv1; 
end
save ('finalresult.mat','archive');
fprintf('The hypervolume value for the final pareto front is %e \n', archiveHV);
figure;
plot(archive(:,V + 1),archive(:,V + 2),'*', plotChromosome(:,V + 1),plotChromosome(:,V + 2),'r.');
title('Final chromosome');
%--------------------------------------------------------------------------%
else
%% if this is not the last iteration calculate hypervolume of the current/temporary solution
TChromosome = [];
TChromosome = FM;
br9 = 0;
msize1= size (TChromosome,1);
for i = 1 : msize1
            if (TChromosome (i,V+3) == 1)
            tempChromosome (i-br9,:) = TChromosome (i,:);
            else
            br9=br9+1;    
            end
end
paretoCT = [];
hv = 0; hv1 = 0; 
for w=1 : M
paretoCT(:,w)= tempChromosome(:,V+w);
end    
limit = [12,0];
[hv,hv1] = hypervolume (paretoCT,limit);
Volume1 (mI,1) = hv;
Volume2 (mI,1) = hv1;
Volume2 (mI,2) = thv;
%---------------------------------------% 
%% If the current solution is the best one than save it into the archive
if mI == 1
archive = [];
archive = TChromosome;
archiveHV = 0;
archiveHV = hv1;
else
if hv1 > archiveHV 
archive = [];
archive = TChromosome;
archiveHV = hv1; 
end
end
%---------------------------------------% 
%% If this is not the last iteration make a combination of all islands and return it for next iteration. 
%% Count the number of individuals that are part of the first pareto front. Additionally, count the number of individuals which were redistributed by calculating the Euclidian distance.
FM = edistance_isl_par (FM, msize, V);
Gcount = 0;
GcountD = 0;
for brj=1 : msize
    if FM (brj,V+3)==1
        Gcount= Gcount+1;
    end
    if (FM (brj,V+3)==1) && (FM (brj,V+5)==0)
        GcountD= GcountD+1;
    end
end

%% Check if the number of individuals in the first Pareto front (before and after the calculation of the Euclidian distance) is smaller than the size of the population.    
if (Gcount<=pop) && (GcountD<=pop)
    % If both values are smaller (or equal) than the population size copy all
    % individuals from the first pareto front to the output chromosome
    br1 = 0;
    for i = 1 : msize
            if FM (i,V+3) == 1
            rChromosome (i-br1,:) = FM (i,:);
            else
            br1=br1+1;    
            end
    end
    
    % If the number of individuals in the output chromosome is still smaller
    % than pop size, include elements from the second,third ... pareto
    % front.
    if size (rChromosome,1)<pop
        br3 = 0;    
        for i = 1 : msize
            if FM (i,V+3) ~= 1
            alllevels (i-br3,:) = FM (i,:);
            else
            br3=br3+1;    
            end
        end
            maxR = size (alllevels,1);
            while size(rChromosome,1)<pop
            rng shuffle;
            kkw = floor(   (   (     (maxR-1)*rand      )+1          )       ); 
            rChromosome((size(rChromosome,1)+1),:) = alllevels (kkw,:);
            end 
    end
    
  
else
    % If the Euclidian count value is smaller than the population size, copy the best       
    % individuals after the Euclidian distance calculation in the output chromosome.
    % Else, If the count of the Euclidian count value is bigger than the popolation size create
    % a pool of all first pareto front Euclidian individuals and randomply copy them
    % in the output chromosome until the size of the output chromosome is equal to the population size 
    if GcountD<=pop 
    br1 = 0;
    for i = 1 : msize
            if (FM (i,V+3) == 1) && (FM (i,V+5) == 0)
            rChromosome (i-br1,:) = FM (i,:);
            else
            br1=br1+1;    
            end
    end
    else
        br8 = 0;
         for i = 1 : msize
            if (FM (i,V+3) == 1) && (FM (i,V+5) == 0)
            firstlevel (i-br8,:) = FM (i,:);
            else
            br8=br8+1;    
            end
        end
    maxT = size (firstlevel,1);
    i = 1;
            while size(rChromosome,1)<pop
            rng shuffle;
            kkw = floor(   (   (     (maxT-1)*rand      )+1          )       ); 
            rChromosome(i,:) = firstlevel (kkw,:);
            i = i +1; 
            end      
    end
            
    % If the size of the output chromosome is still smaller than the pop size 
    % add individuals from the second/third... pareto front untill the size
    % of the output chromosome is equal to the population size. 
    if size (rChromosome,1)<pop
        br4 = 0;    
            for i = 1 : msize
            if ((FM (i,V+3) == 1) && (FM (i,V+5) == 1)) || (FM (i, V+3)~=1)
            addelements (i-br4,:) = FM (i,:);
            else
            br4=br4+1;    
            end
            end
            maxE = size (addelements,1);
            while size(rChromosome,1)<pop
            rng shuffle;
            kkw = floor(   (   (     (maxE-1)*rand      )+1          )       ); 
            rChromosome((size(rChromosome,1)+1),:) = addelements (kkw,:);
            end 
    end
end

%% Visualize the output chromosome  
    figure;
    rChromosome = non_domination_sort_mod(rChromosome, M, V);
    plot(rChromosome(:,V + 1),rChromosome(:,V + 2),'*'); 
    title(sprintf('Selected chromosome for next iteration. Iteration %d: ', mI));
     %  for i = 1 : pop
      % rChromosome(i,V + 1: K) = evaluate_objective_hpm(rChromosome(i,1:V), M, V);
       %end  
       
       %figure;
    %rChromosome = non_domination_sort_mod(rChromosome, M, V);
    %plot(rChromosome(:,V + 1),rChromosome(:,V + 2),'*'); 
    %title(sprintf('After. Iteration %d: ', mI));
fprintf('END OF THE ITERATION No: %d \n',mI);   
end

end

tf=etime(clock,ti); % End time

fprintf('Time: %d seconds\n',tf);
end    
