function [f] = initialize_essex_seg(pop, M, V, seg, lev, elect);
 
% This function initializes the chromosomes. Each chromosome has the
% following at this stage
%       * set of decision variables
%       * objective function values
% 
% where,
% N - Population size
% M - Number of objective functions
% V - Maximum Number of decision variables per segment
% seg - Maximum value for Segment (1 to 20)
% lev - Maximum value for Level (1 to 7)
% elect - Maximum value for Electrode (1 to 15)
% min_range - Vector of minimum values for each decision variable.
% max_range - Vector of maximum values for each decision variable.

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
% Adapted for multiobjective feature selection by Andres Ortiz and Julio
% Ortega

%% Initialize each chromosome
% For each chromosome perform the following (N is the population size)
K=V+M; % Up to V features by segment
W=round(V/seg);
%% Uncompressed (comment the following line if applied)
%% cod=5;
%% Compressed (comment the following line if applied)
cod=4;
f(1:pop,1:K,1:cod)=0;

for i = 1 : pop   % For each individual in the population 
    % Initialize the decision variables based on the minimum and maximum
    % possible values. V is the maximum number of selected features. A random
    % number is picked between the minimum and maximum possible values for
    % the each decision variable.
    
    Features=[];
    for jj = 1 : seg
      for j = 1 : W        
         f(i,(jj-1)*W+j,1) = jj; % segment LDA
         if (randi(100)<50)
             f(i,(jj-1)*W+j,2)=1+round((lev-1)*rand(1));
             f(i,(jj-1)*W+j,3)=1+round((elect-1)*rand(1));
             f(i,(jj-1)*W+j,4)=round(rand(1));             

             % Supervised compressed SEG
             [Features Classes]=gentofen_class_compress_g(f(i,(jj-1)*W+j,1),f(i,(jj-1)*W+j,2),f(i,(jj-1)*W+j,3),f(i,(jj-1)*W+j,4),Features);     
         end
      end   
    end
    
    FeatSeg(1:seg)=0;
    Nelect(1:elect)=0;
    for j = 1 : V
        if (f(i,j,2) > 0)
            FeatSeg(f(i,j,1))=FeatSeg(f(i,j,1))+1;
            Nelect(f(i,j,3))=1;
        end 
    end    
        
    % Unsupervised
    % f(i,V + 1: K,1) = evaluate_essex(Features, M);
    % Supervised
    % f(i,V + 1: K,1) = evaluate_essex_seg(Features, Classes, FeatSeg, Nelect, M);    
    % Supervised (overfitting cost function)
    f(i,V + 1: K,1) = evaluate_essex_seg_ovf(Features, Classes, FeatSeg, Nelect, M);    

end 
%%% save f.mat f;
 
%for i = 1 : N
%  f(i,V + 1: K) = evaluate_objective_hpm_essex(f(i,:), M, V);
%end    
    

end