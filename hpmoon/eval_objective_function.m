function f = eval_objective_function(x, M, V)
% -------------------------------------------------------------------------
% ||||||||||||||||||||||||| GENERAL INFORMATION |||||||||||||||||||||||||||
% -------------------------------------------------------------------------
% FUNCTION EVALUATE_OBJECTIVE_HPM evaluates the objective functions for the
% given input vector 'x'. 'x' is an array of decision variables and f(1),
% f(2), etc. are the objective functions. The algorithm always minimizes 
% the objective function, hence if you would like to maximize the function
% multiply the function by -1. 'M' is the number of objective functions and
% 'V' is the number of decision variables. Training patterns can be found 
% in "datos.mat"
%
% This function is basically written by the user in order to define his/her
% own objective function. Make sure that 'M' and 'V' match your initial 
% user input. 
%
% ******************************* EXAMPLE *********************************
% An example objective function is given below. It has two six decision
% variables and two objective functions.
%
% f = [];
% % - Objective function one
% % Decision variables are used to form the objective function.
% f(1) = 1 - exp(-4*x(1))*(sin(6*pi*x(1)))^6;
% sum = 0;
% for i = 2 : 6
%     sum = sum + x(i)/4;
% end
% % - Intermediate function
% g_x = 1 + 9*(sum)^(0.25);
% 
% % - Objective function two
% f(2) = g_x*(1 - ((f(1))/(g_x))^2);
% *************************************************************************
%
% - INPUT VARIABLES:
%   |_'x': Array of decision variables.
%   |_'M': Number of objective functions.
%   |_'V': Number of decision variables.
%
% - OUTPUT VARIABLES:
%   |_'f': Chromosome + value of the M objective functions once evaluated.
%
% -------------------------------------------------------------------------
% |||||||||||||||||||||| COPYRIGHT AND AUTHORSHIP |||||||||||||||||||||||||
% -------------------------------------------------------------------------
% - Adapted for multiobjective feature selection by:
%   |_ Dr. Prof. Andres Ortiz: 
%       * Entity:  Communications Engineering Department. University of 
%                  Málaga, Spain.
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
%   |_ Date: 12/05/2014 (mm/dd/yy). 
%   |_ Location: Research Centre for Information and Communications 
%                Technologies of the University of Granada, Spain. 
%                (CITIC-UGR).
%
% -------------------------------------------------------------------------
% |||||||||||||||||||||| CODE STARTS RIGHT BELOW ||||||||||||||||||||||||||
% -------------------------------------------------------------------------

% 1) LOAD DATA AND SELECT ACTIVE FEATURES.
% -------------------------------------------------------------------------
% Training data are obtained from datos.mat.
features = importdata('data.mat'); 

% Select active features.
active_features = features(:, x == 1);

% 2) TRAINING OF THE SOM BY USING THE SELECTED FEATURES.
% -------------------------------------------------------------------------
% Normalization of vectors; alternative: var.
features_n = som_normalize(active_features, 'logistic');    

% Initialization of the SOM.
sMap = som_lininit(features_n, 'hexa', 'msize', [20 20]);

% Training of SOM.
sMap = som_batchtrain(sMap, features_n, 'radius_ini', 11, ...
    'radius_fin', 1, 'trainlen', 100);

% 3) DETERMINE NUMBER OF TIMES A UNIT OF THE MAP HAS BEEN ACTIVATED.
% -------------------------------------------------------------------------
% Find hits in the map.
[hits] = som_hits(sMap, features_n);

% Determine the position of the active units.
indices = find(hits);

% 4) COMPUTE MAXIMUM DISTANCE BETWEEN ACTIVE UNITS.
% -------------------------------------------------------------------------
% Initialize value of maximum distance and matrix containing distances.
MAX = 0;
dist = zeros(size(sMap.codebook, 1), size(sMap.codebook, 1));

% Compute distance between active units (one to one).
for ii = 1 : size(indices)
    for jj = 1 : size(indices)
        if (indices(ii) ~= indices(jj))
            dist(indices(ii), indices(jj)) = norm(sMap.codebook(...
                indices(ii), :) - sMap.codebook(indices(jj), :));
            if (dist(indices(ii), indices(jj)) > MAX) 
                MAX = dist(indices(ii), indices(jj));
            end
        end   
    end
end

% Populate values of the diagonal and rows and columns of units which were
% not active with a distance higher than the maximum distance. 
MAX = MAX + 1;
for ii = 1 : size(sMap.codebook, 1)
    dist(ii, ii) = MAX;
    if (hits(ii) == 0)
        for jj = 1 : size(sMap.codebook, 1)
            dist(ii, jj) = MAX;
            dist(jj, ii) = MAX;
        end
    end
end

% 5) COMPUTE OBJECTIVE FUNCTIONS.
% -------------------------------------------------------------------------
% Build a histogram calculating the shortest distances between points (once
% a point has been considered it is not reconsidered).
DD(1, 1) = indices(1);
DD(1, 2) = 0;
dist(:, DD(1, 1)) = MAX;

for ii = 2 : size(indices)
    MIN = min(dist(DD(ii - 1, 1), :));
    test = 0;
    for jj = 1 : size(indices)
        if (test == 0)
            if (dist(DD(ii - 1, 1), indices(jj)) == MIN)
                DD(ii, 1) = indices(jj);
                DD(ii, 2) = MIN;
                test = 1;
            end
        end
    end
    for kk = 1 : size(sMap.codebook, 1)
        dist(kk, DD(ii, 1)) = MAX;
    end
end

% Compute the average distance of all distances.
AVG = mean(DD(:, 2));

% Initialize value of the sum of distances which are below the average 
% distance.
C1 = 0; 
% Initializae number of distances below the average distance.
N1 = 0;
% Initialize value of the sum of distances which are over the average 
% distance.
C2 = 0;
% Initialize number of distances over the average distance.
N2 = 0;

% Compute sum and number of distances both below and over the average 
% distance and number .
for ii = 1 : size(indices)
    if (DD(ii, 2) < AVG)
        C1 = C1 + DD(ii, 2);
        N1 = N1 + 1;
    else
        C2 = C2 + (DD(ii, 2));
        N2 = N2 + 1;
    end
end

% Normalize sums of distances. We aim to minimize C1 (sum of the distances
% below the average minimum distance) and to maximize C2 (sum of the
% distances over the average minimum distance), this is why C2 has a minus
% sign.
C1 = C1 / N1;
C2 = -C2 / N2;

% Initialize vector containing the value of the objective functions.
f = zeros(1, M);

% Populate vector of objective functions.
f(1) = C1;
f(2) = C2;

% In case M was set higher than 2, add random values to the vector of
% objective functions.
for i = 3 : M 
    f(i) = 10 * rand(1);
    fprintf('f(%d)=%d\n', i, f(i));
end

% -------------------------------------------------------------------------
% ||||||||||||||| END OF EVALUATE_OBJECTIVE_HPM FUNCTION ||||||||||||||||||
% -------------------------------------------------------------------------