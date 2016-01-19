function f  = genetic_operator_hpm(parent_chromosome, M, V, mu, mum, l_limit, u_limit);
% genetic_operator_hpm_par
global blvect
global cent
global labels
global image
global image_seg
%% function f  = genetic_operator(parent_chromosome, M, V, mu, mum, l_limit, u_limit)
% 
% This function is utilized to produce offsprings from parent chromosomes.
% The genetic operators corssover and mutation which are carried out with
% slight modifications from the original design. For more information read
% the document enclosed. 
%
% parent_chromosome - the set of selected chromosomes.
% M - number of objective functions
% V - number of decision varaiables
% mu - distribution index for crossover (read the enlcosed pdf file)
% mum - distribution index for mutation (read the enclosed pdf file)
% l_limit - the lower limit for the corresponding decsion variables
% u_limit - the upper limit for the corresponding decsion variables
%
% The genetic operation is performed only on the decision variables, that
% is the first V elements in the chromosome vector. 

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

[N,m] = size(parent_chromosome)
save parent_chromos.txt parent_chromosome -append -ASCII;

clear m;
p = 1;
% Flags used to set if crossover and mutation were actually performed. 
was_crossover = 0;
was_mutation = 0;


for i = 1 : N
    % With 90 % probability perform crossover
    if rand(1) < 0.9
        % Initialize the children to be null vector.
        child_1 = [];
        child_2 = [];
        % Select the first parent
        parent_1 = round(N*rand(1));
        if parent_1 < 1
            parent_1 = 1;
        end
        % Select the second parent
        parent_2 = round(N*rand(1));
        if parent_2 < 1
            parent_2 = 1;
        end
        % Make sure both the parents are not the same. 
        while isequal(parent_chromosome(parent_1,:),parent_chromosome(parent_2,:))
            parent_2 = round(N*rand(1));
            if parent_2 < 1
                parent_2 = 1;
            end
        end
        % Get the chromosome information for each randomnly selected
        % parents
        parent_1 = parent_chromosome(parent_1,:);
        parent_2 = parent_chromosome(parent_2,:);
        % Perform corssover for each decision variable in the chromosome.
        cont1_1=0;
        cont1_2=0;
        for j = 1 : V
            % SBX (Simulated Binary Crossover).
            % For more information about SBX refer the enclosed pdf file.
            % Generate a random number
            u(j) = rand(1);
            if u(j) <= 0.5
                bq(j) = (2*u(j))^(1/(mu+1));
            else
                bq(j) = (1/(2*(1 - u(j))))^(1/(mu+1));
            end
            % Generate the jth element of first child
            child_1(j) = ...
                round(0.5*(((1 + bq(j))*parent_1(j)) + (1 - bq(j))*parent_2(j)));
            % Generate the jth element of second child
            child_2(j) = ...
                round(0.5*(((1 - bq(j))*parent_1(j)) + (1 + bq(j))*parent_2(j)));
            % Make sure that the generated element is within the specified
            % decision space else set it to the appropriate extrema.
            if child_1(j) > u_limit
                child_1(j) = u_limit;
            elseif child_1(j) < l_limit
                child_1(j) = l_limit;
            end
            if child_2(j) > u_limit
                child_2(j) = u_limit;
            elseif child_2(j) < l_limit
                child_2(j) = l_limit;
            end
            if (child_1(j)==1)
                cont1_1=1;
            end
            if (child_2(j)==1)
                cont1_2=1;
            end
        end
        if (cont1_1==0)
            k=round(V*rand(1));
            if (k<1)
                k=1;
            end
            child_1(k)=1;
        end   
        if (cont1_2==0)
            k=round(V*rand(1));
            if (k<1)
                k=1;
            end
            child_2(k)=1;
        end   
        % Evaluate the objective function for the offsprings and as before
        % concatenate the offspring chromosome with objective value.
        for k = 1 : V
            if (child_1(k)<=0.5)
                child_1(k)=round(0);
            else
                child_1(k)=round(1);
            end
        end
        for k = 1 : V
            if (child_2(k)<=0.5)
                child_2(k)=round(0);
            else
                child_2(k)=round(1);
            end
        end        
        %child_1(:,V + 1: M + V) = evaluate_objective_hpm(child_1, M, V);
        %child_2(:,V + 1: M + V) = evaluate_objective_hpm(child_2, M, V);
        child_1(:,V + 1: M + V) = round(0);
        child_2(:,V + 1: M + V) = round(0);
        % Set the crossover flag. When crossover is performed two children
        % are generate, while when mutation is performed only only child is
        % generated.
        was_crossover = 1;
        was_mutation = 0;
    % With 10 % probability perform mutation. Mutation is based on
    % polynomial mutation. 
    else
        % Select at random the parent.
        parent_3 = round(N*rand(1));
        if parent_3 < 1
            parent_3 = 1;
        end
        % Get the chromosome information for the randomnly selected parent.
        child_3 = parent_chromosome(parent_3,:);
        cont1_3=0;
        % Perform mutation on eact element of the selected parent.
        for j = 1 : V
           r(j) = rand(1);
           if r(j) < 0.5
               delta(j) = (2*r(j))^(1/(mum+1)) - 1;
           else
               delta(j) = 1 - (2*(1 - r(j)))^(1/(mum+1));
           end
           % Generate the corresponding child element.
           child_3(j) = child_3(j) + delta(j);
           % Make sure that the generated element is within the decision
           % space.
           if child_3(j) > u_limit
               child_3(j) = u_limit;
           elseif child_3(j) < l_limit   % here it was 'else'
               child_3(j) = l_limit;
           end
           if (child_3(j)==1)
               cont1_3=1;
           end
        end
        if (cont1_3==0)
            k=round(V*rand(1));
            if (k<1)
                k=1;
            end
            child_3(k)=1;
        end
        
        % Evaluate the objective function for the offspring and as before
        % concatenate the offspring chromosome with objective value.  
        for k = 1 : V
            if (child_3(k)<=0.5)
                child_3(k)=round(0);
            else
                child_3(k)=round(1);
            end
        end        
        %child_3(:,V + 1: M + V) = evaluate_objective_hpm(child_3, M, V);
        child_3(:,V + 1: M + V) = round(0);
        % Set the mutation flag
        was_mutation = 1;
        was_crossover = 0;
    end
    % Keep proper count and appropriately fill the child variable with all
    % the generated children for the particular generation.
    if was_crossover
        child(p,:) = child_1;
        child(p+1,:) = child_2;
        was_cossover = 0;
        p = p + 2;
    elseif was_mutation
        child(p,:) = child_3(1,1 : M + V);
        was_mutation = 0;
        p = p + 1;
    end
end


RR=size(child,1);

for ii = 1 : RR
    child(ii,V + 1 : M + V)=evaluate_objective_hpm(child(ii,1 : V), M, V);
end

RR=size(child,2)
save child_genetic.txt child -append -ASCII;

f = child;
