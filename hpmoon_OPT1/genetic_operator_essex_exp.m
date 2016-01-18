function f  = genetic_operator_essex_exp(parent_chromosome, M, V, mu, mum, seg, lev, elect,generation) % changed in essex
% genetic_operator_hpm_par
global blvec
global cent
global labels
global image
global image_seg
%% Based on function f  = genetic_operator(parent_chromosome, M, V, mu, mum, l_limit, u_limit)
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
% l_limit - the lower limit for the corresponding decision variables
% u_limit - the upper limit for the corresponding decision variables
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
%
% Moddified for multiobjective feature selection with BCI datasets of
% University if Essex by Julio Ortega

[N,m,cod] = size(parent_chromosome); % changed in essex version (cod)
%---- save parent_chromos.txt parent_chromosome -append -ASCII;

clear m;
p = 1;
% Flags used to set if crossover and mutation were actually performed. 
was_crossover = 0;
was_mutation = 0;

for i = 1 : N
    % With 90% probability perform crossover (change this in the experiments) 
    % if rand(1) < 0.9
    if rand(1) < 0.5   % 50% probability
        % Selection of parents
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
        while isequal(parent_chromosome(parent_1,:,:),parent_chromosome(parent_2,:,:)) % changed in essex version
            parent_2 = round(N*rand(1));
            if parent_2 < 1
                parent_2 = 1;
            end
        end
        % Crossover implementation in the essex version
        child_1(:,:)= parent_chromosome(parent_1,:,:); % temporal child 1 in essex version
        child_2(:,:)= parent_chromosome(parent_2,:,:); % temporal child 2 in essex version
        
        % Changes in the individuals (essex version)
        cont1_1=round(0.5*V*rand(1))+1;  % Maximum number of changes from parent_1 to parent_2 and viceversa
        for changes = 1 : cont1_1
            cont1_2=round((V-1)*rand(1))+1;
        %   temporal_chromosome(:)=child_1(cont1_2,:);
            for jjj= 1 : cod 
              temporal_chromosome(jjj)=child_1(cont1_2,jjj);
              child1(cont1_2,jjj)= child_2(cont1_2,jjj);
              child2(cont1_2,jjj)=temporal_chromosome(jjj);
            end  
        end    
        child_1(V + 1: M + V,:) = round(0);
        child_2(V + 1: M + V,:) = round(0);
        % Set the crossover flag. When crossover is performed two children
        % are generate, while when mutation is performed only only child is
        % generated.
        was_crossover = 1;
        was_mutation = 0;
    % Perform mutation.
    else
        % Select at random the parent.
        parent_3 = round(N*rand(1));
        if parent_3 < 1
            parent_3 = 1;
        end
        % Get the chromosome information for the randomnly selected parent.
        child_3(:,:) = parent_chromosome(parent_3,:,:); % Changed in essex version
        if (generation < 6) 
            cont1_3=round(V+0.2*V*randn(1));
        else
            cont1_3=round(0.7*V+0.2*V*randn(1));  % Parameters can be changed
        end
        if (cont1_3 < 1)
            cont1_3=1;
        end
        if (cont1_3 > V)
            cont1_3=V;
        end 
        indmut=randperm(V,cont1_3);
        for j = 1 : cont1_3
            for jj = 1 : cod
                   if jj == 1 
                       child_3(indmut(j),jj)= round(0.5*seg+0.2*seg*randn(1));
                       if (child_3(indmut(j),jj) < 1)
                           child_3(indmut(j),jj)=1;
                       end    
                       if (child_3(indmut(j),jj) > seg)
                           child_3(indmut(j),jj)=seg;
                       end    
                   end
                   if jj == 2 
                       child_3(indmut(j),jj)= round((lev-1)*rand(1))+1;
                   end    
                   if jj == 3 
                       child_3(indmut(j),jj)= round((elect-1)*rand(1))+1;
                   end
                   if jj == 4 
                       child_3(indmut(j),jj)= round(rand(1));
                   end    
                   if jj == 5 
                       component=2^(9-lev);
                       child_3(indmut(j),jj)= round((component-1)*rand(1))+1;
                   end                       
            end

        end    
        %%% First implementation
        % Perform mutation on each element of the selected parent.
        % for j = 1 : V
        %   r = rand(1);
        %   if r < 0.80
        %       for jj = 1 : cod
        %           if jj == 1 
        %               child_3(j,jj)= round((seg-1)*rand(1))+1;
        %           end
        %           if jj == 2 
        %               child_3(j,jj)= round((lev-1)*rand(1))+1;
        %           end    
        %           if jj == 3 
        %               child_3(j,jj)= round((elect-1)*rand(1))+1;
        %           end
        %           if jj == 4 
        %               child_3(j,jj)= round(rand(1));
        %           end    
        %           if jj == 5 
        %               component=2^(9-lev);
        %               child_3(j,jj)= round((component-1)*rand(1))+1;
        %           end                       
        %       end
        %   end
        % end   
        %%% Alternative implementation
        child_3(V + 1: M + V,1) = round(0);  % changed in essex code
        % Set the mutation flag
        was_mutation = 1;
        was_crossover = 0;
    end
    %%%%% 
    %%%%%
    % Keep proper count and appropriately fill the child variable with all
    % the generated children for the particular generation.
    if was_crossover
        child(p,:,:) = child_1(:,:); % changed in essex code
        child(p+1,:,:) = child_2(:,:); % changed in essex code 
        was_cossover = 0;
        p = p + 2;
    elseif was_mutation
        child(p,:,:) = child_3(:,:); % changed in essex code
        was_mutation = 0;
        p = p + 1;
    end
end

% ---------------------- from here is has been changed (essex code) 
RR=size(child,1);

for ii = 1 : RR % changed in essex version
    Features=[];
    for jj = 1 : V
         if (child(ii,jj,1)>0) 
             % Unsupervised uncompressed
             % Features=gentofen(child(ii,jj,1),child(ii,jj,2),child(ii,jj,3),child(ii,jj,4),child(ii,jj,5),Features);
             % Supervised uncompressed
             % [Features Classes]=gentofen_class(child(ii,jj,1),child(ii,jj,2),child(ii,jj,3),child(ii,jj,4),child(ii,jj,5),Features);
             % Supervised compressed ***************
             % [Features Classes]=gentofen_class_compress(child(ii,jj,1),child(ii,jj,2),child(ii,jj,3),child(ii,jj,4),Features);
             % Supervised compressed (with global data and classes files)
             [Features Classes]=gentofen_class_compress_g(child(ii,jj,1),child(ii,jj,2),child(ii,jj,3),child(ii,jj,4),Features);

         end
    end
    % data_patterns=Features(1:(size(Features,1)-2),1:size(Features,2));
   
    % Unsupervised
    % child(ii,V + 1: M+V,1) = evaluate_essex(Features, M);
    % Supervised
    channels=child(ii,:,3);
    child(ii,V + 1: M+V,1) = evaluate_essex_lda(Features, Classes, channels, M);
end

f = child;
