function [f, err, ind] = classes_perf_hpm(V, Nclas)
global features;
global labels;
global selec;
% FUNCTION CLASSES_PERF_HPM evaluates the classification performance for 
% the vector of selected features. 
%
% - INPUT PARAMETERS:
%   |_ 'M': number of objective functions.
%   |_ 'V': number of decision variables. 
%
% - OUTPUT PARAMETERS:
%   |_ 'f':
%   |_ 'err':
%   |_ 'ind':
%
% -------------------------------------------------------------------------
% |||||||||||||||||||||| COPYRIGHT AND AUTHORSHIP |||||||||||||||||||||||||
% -------------------------------------------------------------------------
%  - Copyright (c) 2014, Julio Ortega, Andrés Ortiz and Alberto Olivares.
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
% - Original Authors:
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
%   |_ Date: 12/16/2014 (mm/dd/yy). 
%   |_ Location: Research Centre for Information and Communications 
%                Technologies of the University of Granada, Spain. 
%                (CITIC-UGR).
%
% -------------------------------------------------------------------------
% |||||||||||||||||||||| CODE STARTS RIGHT BELOW ||||||||||||||||||||||||||
% -------------------------------------------------------------------------

% 1) LOAD SOLUTIONS TO BE EVALUATED.
% -------------------------------------------------------------------------
solution = importdata('solution.mat');

% 2) LOAD CHROMOSOMES ONE BY ONE.
% -------------------------------------------------------------------------
for ijk = 1 : size(solution, 1) 
    
    x = solution(ijk, 1 : V);

    % 2.1) LOAD FEATURES, CLASS LABELS, TEST DATA AND TEST DATA LABELS.
    % ---------------------------------------------------------------------
    features = importdata('data.mat');
    classes = importdata('classes.mat'); 
    data_test = importdata('data_test.mat');
    classes_test = importdata('classes_test.mat');

    % Get features according to active genes and normalize them. 
    active_features = features(:, x == 1);  
    features_n = som_normalize(active_features, 'var');  
    
    % Get test data according to active genes and normalize them.
    active_data_test = data_test(:, x == 1);  
    data_test_n = som_normalize(active_data_test, 'var');

    % 2.2) SOM COMPUTATION.
    % ---------------------------------------------------------------------
    % Initialize SOM.
    sMap = som_lininit(features_n, 'hexa', 'msize', [10 10]);

    % Train SOM.
    sMap = som_batchtrain(sMap, features_n, 'radius_ini', 5, ...
        'radius_fin', 1, 'trainlen', 60, 'tracking', 3);

    % Get SOM quality and topology characteristics.
    [quality, topology] = som_quality(sMap, features_n);

    % 2.3) EVALUATION OF IDB.
    % --------------------------------------------------------------------
    % Evaluation of IDB (Index of Davies-Bouldin) from the sMap (without
    % k-means clustering) from Voronoi sets.
    [BMU, QERR] = som_bmus(sMap, features_n);
    hits = som_hits(sMap, features_n);

    % Initialize variables. 
    S(1: size(hits, 1)) = 0;
    center(1 : size(hits, 1), 1 : size(features_n, 2)) = 0;
    B(1: size(hits, 1), 1 : size(hits, 1)) = 0;
    IDB = 0;
    kk = 0;

    % Assign classes to units.
    uni_class(1 : max(BMU), 1 : Nclas) = 0;
    
    for i = 1 : size(features_n, 1)
        uni_class(BMU(i), classes(i)) = uni_class(BMU(i), classes(i)) + 1;
    end

    [class_weight, class_estimated] = max(uni_class, [], 2);
    for i = 1 : size(uni_class, 2)
        bi_class(:, i) = uni_class(:, size(uni_class, 2) - i + 1);
    end    
    [biclass_weight, biclass_estimated] = max(bi_class, [], 2);

    search_class(1 : size(uni_class, 1), 1) = 0;
    for i = 1 : size(uni_class, 1)
        if (class_weight(i) == 0)
            class_estimated(i) = 0;
            search_class(i) = 1;
        else
            if (class_estimated(i) ~= (size(uni_class, 2) -  ...
                    biclass_estimated(i) + 1))
               search_class(i) = 1;
            end
        end
    end

    units_assigned = find(class_estimated);
    units_toassign = find(search_class);

    for i = 1 : size(units_toassign, 1)
        for j = 1 : size(units_assigned)
            distance = norm(sMap.codebook(units_toassign(i), :) - ...
                sMap.codebook(units_assigned(j), :));
            if (j == 1)
                MIN = distance;
                class_estimated(units_toassign(i)) = ...
                    class_estimated(units_assigned(j));
            else
                if (distance < MIN)
                    MIN = distance;
                    class_estimated(units_toassign(i)) = ...
                        class_estimated(units_assigned(j));
                end
            end 
        end  
    end
    %  End of classes assignment to units. --------------------------------

    for i = 1 : size(BMU, 1)
        for j = 1 : size(features_n, 2)
            center(BMU(i), j) = center(BMU(i), j) + ...
                (features_n(i, j) / hits(BMU(i)));
        end
    end

    % Compute centroids of clusters. 
    figure
    subplot(1, 1, 1)
    hold on;
    for i = 1 : size(features_n, 1)
        if (classes(i) == 1) 
            plot(features_n(i, 1), features_n(i, 2), 'ok');
        end    
        if (classes(i) == 2) 
            plot(features_n(i, 1), features_n(i, 2), 'og');
        end    
        if (classes(i) == 3) 
            plot(features_n(i, 1), features_n(i, 2), 'ob');
        end    
        if (classes(i) == 4) 
            plot(features_n(i, 1), features_n(i, 2), '*k');
        end   
    end
    plot(center(:, 1), center(:, 2), '+r');

    for i = 1 : size(sMap.codebook, 1) % units of the map
        for ii = 1 : size(features_n, 1) % input patterns 
            if (BMU(ii) == i)
                for jj = 1 : size(center, 2)
                    vector(jj) = features_n(ii, jj) - center(i, jj);
                end    
                S(i) = S(i) + (norm(vector) / hits(i));
            end
        end
    end

    for i = 1 : size(sMap.codebook, 1)
        for ii = 1 : size(sMap.codebook, 1)  
            if ((hits(i) ~= 0) && (hits(ii) ~= 0))  
                B(i, ii) = norm(center(i) - center(ii));
            end
        end
    end

    for i = 1 : size(sMap.codebook, 1)
        R(i) = 0;
        for ii = 1 : size(sMap.codebook, 1) 
            if (i ~= ii)
                if (B(i, ii) ~= 0)
                    SS = (S(i) + S(ii)) / B(i, ii);
                    if (R(i) < SS)
                        R(i) = SS;
                    end
                end    
            end   
        end
    end

    for i = 1 : size(sMap.codebook, 1)
        IDB = IDB + R(i);
        if (R(i) ~= 0)
            kk = kk + 1;
        end
    end
    IDB = IDB / kk;

    % 2.4) CLASSIFICATION PERFORMANCE (LABELED PATTERNS).
    % ---------------------------------------------------------------------
    correct = 0;
    incorrect = 0;
    TN(1 : Nclas) = 0;
    TP(1 : Nclas) = 0;
    FN(1 : Nclas) = 0;
    FP(1 : Nclas) = 0;
    ACCUR(1 : Nclas) = 0;
    SENS(1 : Nclas) = 0;
    SPEC(1 : Nclas) = 0;
    JACCARD(1 : Nclas) = 0;
    MatrixKappa(1 : Nclas, 1 : Nclas) = 0;
    NKappa(1 : Nclas) = 0;
    MKappa(1 : Nclas) = 0;
    ProbClases(1 : Nclas) = 0;
    ProbClasif(1 : Nclas) = 0;
    Total = size(data_test_n, 1);
    Po = 0;
    Pc = 0;

    [BMU_class, QERR_class] = som_bmus(sMap, data_test_n);
    
    % Kappa coefficient computation
    for i = 1 : size(data_test_n, 1)
        if (class_estimated(BMU_class(i)) == classes_test(i))
            correct = correct + 1;
        else
            incorrect = incorrect + 1;
        end
        MatrixKappa(class_estimated(BMU_class(i)), classes_test(i)) = ...
           MatrixKappa(class_estimated(BMU_class(i)), classes_test(i)) + 1;
        MKappa(class_estimated(BMU_class(i))) =  ...
            MKappa(class_estimated(BMU_class(i))) + 1;
        NKappa(classes_test(i)) = NKappa(classes_test(i)) + 1;
    end

    for i = 1 : Nclas
        ProbClases(i) = NKappa(i) / Total;
        ProbClasif(i) = MKappa(i) / Total;
        Po = Po + MatrixKappa(i, i) / Total;
        Pc = Pc + ProbClases(i) * ProbClasif(i);
    end
    Kappa = (Po - Pc)/(1 - Pc);

    for i = 1 : size(data_test_n, 1)
        if (class_estimated(BMU_class(i)) == classes_test(i))
            for j = 1 : Nclas
                if (j == classes_test(i)) 
                    TP(j) = TP(j) + 1;
                else
                    TN(j) = TN(j) + 1;
                end
            end
        else
            FP(class_estimated(BMU_class(i))) = ...
                FP(class_estimated(BMU_class(i))) + 1;
            for j = 1 : Nclas
                if (j == classes_test(i))
                    FN(j) = FN(j) + 1;
                end
            end
        end
    end       

    for i = 1 : Nclas
        ACCUR(i) = (TP(i) + TN(i)) / (TP(i) + TN(i) + FP(i) + FN(i));
        SENS(i) = TP(i) / (TP(i) + FN(i));
        SPEC(i) = TN(i) / (TN(i) + FP(i));
        JACCARD(i) = TP(i) / (TP(i) + TN(i) + FP(i));
    end

    truefalse = [TP; TN; FP; FN];
    JACC = [ACCUR; SENS; SPEC; JACCARD];

    % 2.5) NEX INDEX CALCULATION (MAP UNITS).
    % ---------------------------------------------------------------------
    hits = som_hits(sMap, features_n);
    indices = find(hits);
    MAX = 0;
    dist = zeros(size(sMap.codebook, 1), size(sMap.codebook, 1));

    for ii = 1 : size(indices)
        for jj = 1 : size(indices)
            if indices(ii) ~= indices(jj)
                dist(indices(ii), indices(jj)) = ...
                    norm(sMap.codebook(indices(ii), :) -  ...
                    sMap.codebook(indices(jj), :));
                if dist(indices(ii), indices(jj)) > MAX
                    MAX = dist(indices(ii), indices(jj));
                end
            end   
        end
    end
    
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
                    DD(ii ,2) = MIN;
                    test = 1;
                end
            end
        end
        
        for kk = 1 : size(sMap.codebook, 1)
            dist(kk, DD(ii, 1)) = MAX;
        end
    end
    
    order_dist = sort(DD(:, 2));
    AVG = mean(DD(:, 2));
    C1 = 0;
    N1 = 0;
    C2 = 0;
    N2 = 0;

    for ii = 1 : size(indices)
        if (DD(ii, 2) < AVG)
            C1 = C1 + DD(ii, 2);
            N1 = N1 + 1;
        else
            C2 = C2 + (DD(ii, 2));
            N2 = N2 + 1;
        end
    end
    C1 = C1 / N1;
    C2 = -C2 / N2;
    cut = N1;
 
    % 2.6) NEW INDEX CALCULATION (DATA).
    % ---------------------------------------------------------------------
    MAX_d = 0;
    N_d = size(features_n, 1);
    for ii = 1 : N_d
        for jj = 1 : N_d
            if (ii == jj) 
                dist_d(ii, jj) = 0;
            else
                dist_d(ii, jj) = norm(features_n(ii, :) - ...
                    features_n(jj, :));
                if (dist_d(ii, jj) > MAX_d) 
                    MAX_d = dist_d(ii, jj);
                end
            end
        end
    end
    
    MAX_d = MAX_d + 1;
    for ii = 1 : N_d
        dist_d(ii, ii) = MAX_d;
    end
    DD_d(1, 1) = 1;
    DD_d(1, 2) = 0;

    for ii = 2 : N_d
        MIN_d = min(dist_d(DD_d(ii - 1, 1), :));
        test = 0;
        for jj = 1 : N_d
            if (test == 0)
                if (dist_d(DD_d(ii - 1, 1), jj) == MIN_d)
                    DD_d(ii, 1) = jj;
                    DD_d(ii, 2) = MIN_d;
                    test = 1;
                end
            end
        end
        for kk = 1 : N_d
            dist_d(kk, DD_d(ii, 1)) = MAX_d;
        end
    end
    
    length_d = sort(DD_d(:,2));
    AVG_d = mean(DD_d(:,2));
    C1_d = 0;
    N1_d = 0;
    C2_d = 0;
    N2_d = 0;

    for ii = 1 : N_d
        if (DD_d(ii, 2) < AVG_d)
            C1_d = C1_d + DD_d(ii, 2);
            N1_d = N1_d + 1;
        else
            C2_d = C2_d + (DD_d(ii, 2));
            N2_d = N2_d + 1;
        end
    end
    
    C1_d = C1_d / N1_d;
    C2_d = -C2_d / N2_d;
    cut_d = N1_d;

    %hold on;
    %plot(sMap.codebook(:,1),sMap.codebook(:,2),'d');
    %plot(features_n(:,1),features_n(:,2),'+r');
    %hold off
    %title('Data and map after training');

    % 2.7) CLUSTERING BY K-MEANS.
    % ---------------------------------------------------------------------
    % Find clusterings (explores up to 20 clusters).
    [c, p, err, ind] = kmeans_clusters(sMap, 20); 
    
    % Select the one with smallest index.
    [dummy, i] = min(ind); 
    % visualize it
    som_show(sMap,'color',{p{i},sprintf('%d clusters',i)});
    colormap(jet(i)), som_recolorbar % change colormap

    % 2.8) STORE RESULTS.
    % ---------------------------------------------------------------------
    selected_features = features(:, x == 1); 
    resultados = [selected_features quality topology IDB C1 C2 C1_d ...
        C2_d correct incorrect Kappa];
    quality_measures = [quality topology IDB C1 C2 C1_d C2_d correct ...
        incorrect Kappa];
    
    save measure.txt resultados -ASCII -append;
    save qm.txt quality_measures -ASCII -append;
    save TPTNFPFN.txt truefalse -ASCII -append;
    save JACCARD.txt JACC -ASCII -append;

end