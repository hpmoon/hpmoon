function f  = replace_chromosome(intermediate_chromosome, M, V,pop)

% FUNCTION REPLACE_CHROMOSOME replaces the chromosomes based on rank and 
% crowding distance. Initially until the population size is reached each 
% front is added one by one until addition of a complete front which 
% results in exceeding the population size. At this point the chromosomes
% in that front are added subsequently to the population based on crowding
% distance.
%
% - INPUT PARAMETERS:
%   |_ 'intermediate_chromosome': Intermediate chromosome containing
%       original chromosome and offspring chromosome.
%   |_ 'M': Number of objective functions.
%   |_ 'V': Number of decision variables.
%   |_ 'pop': Size of population.
%
% - OUTPUT PARAMETERS:
%   |_ 'f': Replaced chromosome.
%
% -------------------------------------------------------------------------
% |||||||||||||||||||||| COPYRIGHT AND AUTHORSHIP |||||||||||||||||||||||||
% -------------------------------------------------------------------------
% Original algorithm NSGA-II was developed by researchers in Kanpur Genetic
% Algorithm Labarotary and kindly visit their website for more information
% http://www.iitk.ac.in/kangal/

%  Copyright (c) 2009, Aravind Seshadri. All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the 
%       distribution.
%      
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
%  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
%  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
%  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
%  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
%  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
%  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
%  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
%  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
%  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
%  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% - Adapted for multiobjective feature selection by:
%   |_ Dr. Prof. Andres Ortiz: 
%       * Entity:  Communications Engineering Department. University of 
%                  M�laga, Spain.
%       * Contact: aortiz@ic.uma.es
%   |_ Dr. Prof. Julio Ortega: 
%       * Entity:  Department of Computer Architecture and Computer
%                  Technology, University of Granada, Spain.
%       * Contact: jortega@ugr.es
%
% - Code restructuring and optimization:
%   |_ Dr. Alberto Olivares:
%       * Entity:   Department of Signal Theory, Telematics and
%                   Communications, University of Granada, Spain.
%       * Contact:  aolivares@ugr.es
%
% LAST MODIFICATION: 
%   |_ Date: 12/15/2014 (mm/dd/yy). 
%   |_ Location: Research Centre for Information and Communications 
%                Technologies of the University of Granada, Spain. 
%                (CITIC-UGR).
%
% -------------------------------------------------------------------------
% |||||||||||||||||||||| CODE STARTS RIGHT BELOW ||||||||||||||||||||||||||
% -------------------------------------------------------------------------

% Get size of intermediate chromosome.
[N, m] = size(intermediate_chromosome);

% Get the index for the population sort based on the rank.
[temp, index] = sort(intermediate_chromosome(:, M + V + 1));

clear temp m

% Now sort the individuals based on the index.
for i = 1 : N
    sorted_chromosome(i, :) = intermediate_chromosome(index(i), :);
end

% Find the maximum rank in the current population
max_rank = max(intermediate_chromosome(:, M + V + 1));

% Start adding each front based on rank and crowing distance until the
% whole population is filled.
previous_index = 0;
for i = 1 : max_rank
    
    % Get the index for current rank i.e the last element in the
    % 'sorted_chromosome' with rank i. 
    current_index = max(find(sorted_chromosome(:, M + V + 1) == i));
    
    % Check to see if the population is filled if all the individuals with
    % rank i are added to the population. 
    if current_index > pop
        % If so then find the number of individuals with current rank i.
        remaining = pop - previous_index;
        
        % Get information about the individuals in the current rank i.
        temp_pop = ...
            sorted_chromosome(previous_index + 1 : current_index, :);
        
        % Sort the individuals with rank i in descending order based on
        % the crowding distance.
        [temp_sort, temp_sort_index] = sort(temp_pop(:, M + V + 2), ...
            'descend');
        
        % Start filling individuals into the population in descending order
        % until the population is filled.
        for j = 1 : remaining
            f(previous_index + j, :) = temp_pop(temp_sort_index(j), :);
        end
        return;
    elseif current_index < pop
        % Add all the individuals with rank i into the population.
        f(previous_index + 1 : current_index, :) = ...
            sorted_chromosome(previous_index + 1 : current_index, :);
    else
        % Add all the individuals with rank i into the population.
        f(previous_index + 1 : current_index, :) = ...
            sorted_chromosome(previous_index + 1 : current_index, :);
        return;
    end
    % Get the index for the last added individual.
    previous_index = current_index;
end
