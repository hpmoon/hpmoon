function hpmoon_alt3_par(pop,gen,iterat,islands)


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
    error('Minimum number of generations is 5');
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
ef=zeros(pop,V);

%% START OF THE MAIN ITERATIONS
ti=clock;
for mI=1:iterat
if mI == 1
    IChromosome = [];
else
    IChromosome = [];
    IChromosome = rChromosome;
    %ef=rChromosome(:,1:V);
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
[chromosome] = initialize_variables_hpm_isl_alt_par(pop, M, V, min_range, max_range, MAXfeat, chunk, ef, mI);  
%% Sort the initialized population
% Sort the population using non-domination-sort. This returns two columns
% for each individual which are the rank and the crowding distance
% corresponding to their position in the front they belong. At this stage
% the rank and the crowding distance for each chromosome is added to the
% chromosome vector for easy of computation.

chromosome = non_domination_sort_mod(chromosome, M, V);
else
chromosome = IChromosome;  
end


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
endChromosome = [];
endChromosome = chromosome;  
chromosome = chromosome (:,(chunk*MAXfeat)+1:((chunk+1)*MAXfeat));


%% Result
% Save the result in ASCII and MAT format.
parsave_alt(chromosome,chunk,V,mI,endChromosome,iterat);
end

%% PREPARE FOR THE NEXT ITERATION

%% Concatenate the chromosomes from all processors and plot the joint matrix
clear chromosome; clear thv;
paretoC = [];
paretoC1 = [];
plotChromosome = [];
thv=etime(clock,ti);
if mI == iterat
    plotC = [];
    pCh = [];
    for cnum = 0: islands-1
    load (['chromosomeFinal_P-',num2str(cnum),'Ite-',num2str(mI),'.mat']);
    plotC = vertcat (plotC,pCh);
    end
    plotC = non_domination_sort_mod(plotC, M, V);
msize = size (plotC,1);
for i = 1 : msize
            if (plotC (i,V+3) == 1)
            plotChromosome (i,:) = plotC (i,:);    
            end
end
save ('finalresultALT.mat','plotChromosome');
figure;
plot(plotChromosome(:,V + 1),plotChromosome(:,V + 2),'*');
title('Final chromosome');
for w=1 : M
paretoC(:,w)= plotChromosome(:,V+w);
end
limit = [12,0];
[hv,hv1] = hypervolume (paretoC,limit);
fprintf('The hypervolume value for the pareto front is %e /// %e \n',hv,hv1);
Volume1 (mI,1) = hv;
Volume2 (mI,1) = hv1;
Volume2 (mI,2) = thv;
if hv1 > archiveHV 
archive = [];
archive = plotChromosome;
archiveHV = hv1; 
end
fprintf('The hypervolume value for the best pareto front is %e \n',archiveHV);
save ('hypervolume1ALT2.mat','Volume1');
save ('hypervolume2ALT2.mat','Volume2');
else
%% Calculate hypervolume value for this iteration
pCh = [];
plotT = [];
for cnum1 = 0: islands-1
    load (['chromosomeFinal_P-',num2str(cnum1),'Ite-',num2str(mI),'.mat']);
    plotT = vertcat (plotT,pCh);
end
    plotT = non_domination_sort_mod(plotT, M, V);
msize3 = size (plotT,1);
for i = 1 : msize3
            if (plotT (i,V+3) == 1)
            plotChromosome1 (i,:) = plotT (i,:);    
            end
end
for w=1 : M
paretoC1(:,w)= plotChromosome1(:,V+w);
end
limit = [12,0];
[hv,hv1] = hypervolume (paretoC1,limit);
fprintf('The hypervolume value for the pareto front is %e /// %e \n',hv,hv1);
Volume1 (mI,1) = hv;
Volume2 (mI,1) = hv1;
Volume2 (mI,2) = thv;
figure;
plot(plotChromosome1(:,V + 1),plotChromosome1(:,V + 2),'*');
title(sprintf('Temporary pareto, iteration %d: ', mI));
%% If the current solution is the best one than save it into the archive
if mI == 1
archive = [];
archive = plotT;
archiveHV = 0;
archiveHV = hv1;
else
if hv1 > archiveHV 
archive = [];
archive = plotT;
archiveHV = hv1; 
end
end
%---------------------------------------% 
%% Make vertical combination 
    
    for cnum1 = 0: islands-1
    ChooseCh = [];
    chunkCh = [];
    load (['chromosomeFinal_P-',num2str(cnum1),'Ite-',num2str(mI),'.mat']);
    chunkCh = pCh; 
            for sc = 0 : islands -1
                if sc ~=cnum1
                load (['chromosome_P-',num2str(sc),'Ite-',num2str(mI),'.mat']);
                ChooseCh = vertcat (ChooseCh,chromosome);
                end
            end
         
            for mcount = 0 : islands-2
                    for icount = pop :-1: (pop-(islands-1)) 
                    chunkCh (icount,(mcount*MAXfeat)+1:((mcount+1)*MAXfeat)) = ChooseCh ((pop*mcount)+1,:);
                    end
            end
    nextCh (:,(cnum1*MAXfeat)+1:((cnum1+1)*MAXfeat)) = chunkCh (:,(cnum1*MAXfeat)+1:((cnum1+1)*MAXfeat));  
    end 
   rChromosome = nextCh; 
    for i = 1 : pop
    rChromosome(i,V + 1: K) = evaluate_objective_hpm(rChromosome(i,:), M, V);
    end
    rChromosome = non_domination_sort_mod(rChromosome, M, V);
figure;
plot(rChromosome(:,V + 1),rChromosome(:,V + 2),'*');
title(sprintf('Send chromosome, iteration %d: ', mI));
fprintf('END OF THE ITERATION No: %d \n',mI);   
end

end

tf=etime(clock,ti); % End time

fprintf('Time: %d seconds\n',tf);
end    
