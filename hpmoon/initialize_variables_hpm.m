function f = initialize_variables_hpm(N, M, V, MAXfeat, showChromosomes)
% -------------------------------------------------------------------------
% ||||||||||||||||||||||||| GENERAL INFORMATION |||||||||||||||||||||||||||
% -------------------------------------------------------------------------
% FUNCTION INITIALIZE_VARIABLES initializes the chromosomes. Each 
% chromosome has the following at this stage:
%       - Set of decision variables.
%       - Objective function values.
% 
% - INPUT VARIABLES:
%   |_ 'N': Population size (scalar).
%   |_ 'M': Number of objective functions (scalar).
%   |_ 'V': Number of decision variables (scalar).
%   |_ 'MAXfeat':   Maximum number of features (scalar).
%   |_ 'showChromosomes': Flag which indicates whether the initial random
%                         chromosome composition is to be plotted ('yes').
%
% - OUTPUT VARIABLES:
%   |_ 'f': Initialized chromosome (?).
% 
% -------------------------------------------------------------------------
% |||||||||||||||||||||| COPYRIGHT AND AUTHORSHIP |||||||||||||||||||||||||
% -------------------------------------------------------------------------
%  - Copyright (c) 2009, Aravind Seshadri. All rights reserved.
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

% 0) DEFINE AUXILIARY VARIABLES.
% -------------------------------------------------------------------------
% Define 'K'. K is the total number of array elements. For ease of 
% computation decision variables and objective functions are concatenated
% to form a single array. For crossover and mutation only the decision 
% variables are used while for selection, only the objective variable are 
% utilized.
K = M + V;

% Initialize the matrix containing the chromosomes (solutions) in each one
% of its rows. The number of rows equals the size of the population, and
% the number of columns equals the number of features of each chromosome
% (M) plus the number of cost functions (V).
f = zeros(N, K);

% 1) INITIALIZE EACH CHROMOSOME.
% -------------------------------------------------------------------------
% The chromosomes are given the values 0 and 1 randomly. The value 0 means
% that the feature located at the same index is not selected, while 1 means
% it will be selected.
for i = 1 : N
  for j = 1 : V
        f(i, j) = round(0);
  end
  for jj = 1 : MAXfeat
        kk = 1 + round((V - 1) * rand(1));
        f(i, kk) = 1;
  end
    % For ease of computation and data handling the chromosome also has the
    % value of the objective function concatenated at the end. The elements
    % V + 1 to K contain the evaluated objective function. 
    % The function 'evaluate_objective_hpm' takes one chromosome at a time.
    % Only the decision variables are passed to the function along with 
    % information about the number of objective function. The function 
    % returns the value for the objective functions. These
    % values are now stored at the end of the chromosome itself.
  fprintf('Generating chromosome %d out of %d', i, N);  
  f(i, V + 1: K) = eval_objective_function(f(i, :), M, V);
end

% Finally, the generated chromosomes can be plot in a figure if the
% 'showChromosomes' variable is set to 'yes'.
if strcmpi(showChromosomes,'yes')
    % Get the matrix size.
    [r,c] = size(f);
    % Plot the image.
    imagesc((1:c)+0.5,(1:r)+0.5,f(:,1:end-2));
    % Use a gray colormap.
    colormap(gray);
    % Make axes grid sizes equal.
    axis equal
    % Change some axes properties.
    set(gca,'XTick',[1 c-2],'YTick',[1 r], 'XLim',[1 c-2],'YLim', ...
        [1 r], 'GridLineStyle','-','XGrid','on','YGrid','on');
    xlabel('Genes (Features)')
    ylabel('Chromosome number')
end
% -------------------------------------------------------------------------
% ||||||||||||||| END OF INITIALIZE_VARIABLES_HPM FUNCTION ||||||||||||||||
% -------------------------------------------------------------------------