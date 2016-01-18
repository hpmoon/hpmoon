function hpmoon_essex_g(pop,gen,experiment)
global class_essex;
global data_essex;

% function hpmoon_esses_g(pop,gen)
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
% This software has been adapted to the problem of multiobjective feature
% selection by Andres Ortiz and Julio Ortega

% Number of Arguments
% Check for the number of arguments. The two input arguments are necessary
% to run this function.
if nargin < 2
    error('NSGA-II: Please enter the population size and number of generations as input arguments.');
end
% Both the input arguments need to of integer data type
if isnumeric(pop) == 0 || isnumeric(gen) == 0
    error('Both input arguments pop and gen should be integer datatype');
end
% Minimum population size has to be 10 individuals
if pop < 10
    error('Minimum population for running this function is 20');
end
if gen < 5
    error('Minimum number of generations is 5');
end
% Make sure pop and gen are integers
pop = round(pop);
gen = round(gen);

% The objective function description contains information about the
% objective functions. M is the dimension of the objective space, V is the
% dimension of decision variable space.
% The individuals in the population are vectors of V binary components, 
% so min_range=0 and max_range=1
% User has to define the objective functions using the decision variables. 
% Make sure to edit the function 'evaluate_objective_hpm' to suit your needs.

%%%%% g = sprintf('Number of objectives: '); % Uncomment to input manually
% Obtain the number of objective function
%%%%% M = input(g); % Number of Objectives (Uncomment to input manually)
M=2;
if M < 2
    error('This is a multi-objective optimization function hence the minimum number of objectives is two');
end
%%%%% g = sprintf('\nMaximum Number of decision variables (features): ');
% Obtain the number of decision variables (Uncomment for input manually)
%%%%% V = input(g); % Number of componets of the individuals in the population
V=30;
% g = sprintf('\nMaximum number of features in solution (<=%i): ',V);
% removed in the essex version
% Obtain the maximun number of features allowed to start
% MAXfeat = input(g); % removed in the essex version
% MAXfeat=V;


class_essex=importdata('class_essex_x.mat');
data_essex=importdata('data_essex_x.mat');

seg=size(data_essex,1);
lev=size(data_essex,2);
elect=size(data_essex,3);

%g = sprintf('Number of segments (1 to 20)');
%Obtain the number of segments
%seg = input(g); % Number of segments

%g = sprintf('Number of levels (1 to 7)');
% Obtain the number of levels
%lev = input(g); % Number of levels

%g = sprintf('Number of electrodes (1 to 15)');
% Obtain the number of electrodes
%elect = input(g); % Number of electrodes

min_range=round(0); % Individuals with binary components (min_range=0)
max_range=round(1); % Individuals with bynary components (max_range=0)

K=M+V; % Components of the individuals (different features) plus objetives 
selec=zeros(1,K);

%% Initialize the population
% Population is initialized with random values which are within the
% specified range. Each chromosome consists of the decision variables. Also
% the value of the objective functions, rank and crowding distance
% information is also added to the chromosome vector but only the elements
% of the vector which has the decision variables are operated upon to
% perform the genetic operations like corssover and mutation.
ti=clock; % Initial time

chromosome = initialize_essex(pop, M, V, seg, lev, elect); % change in essex version

%% Sort the initialized population
% Sort the population using non-domination-sort. This returns two columns
% for each individual which are the rank and the crowding distance
% corresponding to their position in the front they belong. At this stage
% the rank and the crowding distance for each chromosome is added to the
% chromosome vector for easy of computation.

chromosome = non_domination_sort_essex(chromosome, M, V);

%% Start the evolution process
% The following are performed in each generation
% * Select the parents which are fit for reproduction
% * Perfrom crossover and Mutation operator on the selected parents
% * Perform Selection from the parents and the offsprings
% * Replace the unfit individuals with the fit individuals to maintain a
%   constant population size.

for i = 1 : gen
    
    fprintf('Experiment: %d GENERATION: %d \n',experiment,i);
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
    parent_chromosome = tournament_selection_essex(chromosome, pool, tour);

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
% To evaluate different options for the genetic operator (uncomment the following three lines)   
  offspring_chromosome = ...
       genetic_operator_essex_exp(parent_chromosome, ...
       M, V, mu, mum, seg, lev, elect,i); % changed in essex version
% To evaluate different options for the genetic operator (comment the following three lines)  
% offspring_chromosome = ...
%         genetic_operator_essex(parent_chromosome, ...
%         M, V, mu, mum, seg, lev, elect); % changed in essex version
    % Intermediate population
    % Intermediate population is the combined population of parents and
    % offsprings of the current generation. The population size is two
    % times the initial population.
    
    % fprintf('Size %d \n',i);
    [main_pop,temp,gencode] = size(chromosome);
    [offspring_pop,temp,gencode] = size(offspring_chromosome);
    % temp is a dummy variable.
    %clear temp
    % intermediate_chromosome is a concatenation of current population and
    % the offspring population.
    
    % fprintf('Intermediate %d \n',i);
    intermediate_chromosome(1:main_pop,:,:) = chromosome; % changed in essex version
    %%% intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,1 : M+V,:) = ...
    %%%    offspring_chromosome; % changed in essex version
    intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,:,:) = ...
           offspring_chromosome; % changed in essex version

    % Non-domination-sort of intermediate population
    % The intermediate population is sorted again based on non-domination sort
    % before the replacement operator is performed on the intermediate
    % population.
    
    % fprintf('Non domination sort %d \n',i);
    intermediate_chromosome = ...
        non_domination_sort_essex(intermediate_chromosome, M, V); % changed in essex version
    % Perform Selection
    % Once the intermediate population is sorted only the best solution is
    % selected based on it rank and crowding distance. Each front is filled in
    % ascending order until the addition of population size is reached. The
    % last front is included in the population based on the individuals with
    % least crowding distance
    
    % fprintf('Replace chromosome %d \n',i);
    chromosome = replace_chromosome_essex(intermediate_chromosome, M, V, pop); % changed in essex version
    if ~mod(i,100)
        clc
        fprintf('%d generations completed \n',i);
    end
end

tf=etime(clock,ti); % End time

fprintf('Time: %d seconds\n',tf);
save time.txt tf -append -ASCII;

%% Result
% Save the result in ASCII text format.
% save solution.txt chromosome -ASCII; % Not useful in essex version
%%%%% save solutionXXX.mat chromosome XXX=number of datafile: 101, 102,...;
save (['solXXX',mat2str(experiment),'_',mat2str(gen),'.mat'], 'chromosome');

%%%%% result=chromosome(:,(V+1):size(chromosome,2),1); % changed in essex version
%%%%% save pareto.txt result -ASCII;

% Not useful in essex version
% for i = 1 : size(chromosome,1)
%    xyz=find(chromosome(i,1:V));
%    save features.txt xyz -append -ASCII;
% end

%% Visualize
% The following is used to visualize the result if objective space
% dimension is visualizable.
if M == 2
    plot(chromosome(:,V + 1,1),chromosome(:,V + 2,1),'*'); % changed in essex version
elseif M ==3
    plot3(chromosome(:,V + 1,1),chromosome(:,V + 2,1),chromosome(:,V + 3,1),'*'); % changed in essex version
end
    
