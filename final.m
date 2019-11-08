function [geneFeatures, geneFeaturesP, featuresTableMeta, operations] = final(fdMatrix, fStructInfo, subsetTitles, nanTolerance, regions, operations, isPca)
    % Usage: [gf, gfp, features, operations] = final(TS_DataMat, joinedStructInfo, [2, 4, 6, 8, 9, 12, 14], 0.15, {'Isocortex'}, Operations, false);
    %
    % Correlates gene expression data with time-series features, and creates a violin plot of the results
    % 
    % ------------------ 
    % INPUTS
    % ------------------ 
    % - Matrix of features
    % - Table of structure info (brain areas)
    % - 1xn matrix of which gene subsets to include
    % - Null threshold (e.g. 0.15 = 15%)
    % - Cell of brain regions to be selected
    % - Table of operations info
    % - Integer specifying how many PCA features are used. (= 0 if not using PCA features)
    %
    % ------------------ 
    % OUTPUTS
    % ------------------ 
    % - A feature x subset matrix of Spearman correlation coefficients
    % - A feature x subset matrix of p-values (un-adjusted)
    % - A table containing the outputs for each gene subset
    % - A table containing the outputs for each feature
    %
    % ------------------ 
    % REQUIREMENTS
    % ------------------ 
    % - LoadGeneExpressionData
    % - gd_tsd_merge
    % - filter_regions
    % - remove_nans
    % - create_acronyms
    % - BF_JitteredParallelScatter (https://github.com/benfulcher/hctsa/tree/master/PeripheryFunctions)
    % - distinguishable_colors (Timothy E. Holy, 2011) 
    % https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/38369/versions/1/previews/distinguishable_colors.m/index.html


    % Load gene data
    [gd, geneInfo, geneStructInfo] = LoadGeneExpressionData();

    if isPca
        operations = operations(1:isPca, 1);
    end

    % Merge feature with genes
    [sortedG, sortedFMatrix, structInfo] = gd_tsd_merge(gd, ...
        geneStructInfo, fdMatrix, fStructInfo);
    clear fStructInfo fdMatrix gd geneStructInfo;

    % Only retain brain areas in the `region` region
    [sortedG, sortedFMatrix, structInfo] = filter_regions(sortedG, sortedFMatrix, structInfo, regions);

    % Remove nans with > nanTolerance in both columns and rows
    nans = isnan(sortedG);
    [sortedG, sortedFMatrix, geneInfo, structInfo] = remove_nans(sortedG, ...
        sortedFMatrix, nans, nanTolerance, geneInfo, structInfo, false);
    clear nans;

    % Load in acronyms as matrix
    titles2 = {
    "Serotonin", % 1
    "Serotonin",
    "Dopamine", % 3
    "Dopamine",
    "Glutamate", % 5
    "Glutamate",
    "Acetylcholine", % 7
    "Acetylcholine",
    "GABA", % 9
    "GABA",
    "Histamine", % 11
    "Histamine",
    "Noradrenaline", % 13
    "Neurotransmitter",
    };
    [acronymsList, acronymsTitles] = create_acronyms(subsetTitles);

    % For each set of genes in acronyms, calculate PCA and add to matrix region x gene subset
    numSubsets = size(acronymsList, 2);
    subsetMatrix = zeros(size(sortedG, 1), numSubsets);

    for i = 1:numSubsets
        [~, ~, geneMask] = intersect(acronymsList{i}, geneInfo.acronym);
        subG = sortedG(:, geneMask);

        % Calculate zscore excluding NaNs
        zscore_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x, 'omitnan')), std(x, 'omitnan'));
        [coeff, score] = pca(zscore_xnan(subG)', 'Rows', 'complete');
        subsetMatrix(:, i) = coeff(:, 1);

        if false % Plots the PCA explanation figure
            figure()
            subplot(1, 2, 1)
            imagesc(zscore_xnan(subG))
            subplot(1, 2, 2)
            imagesc(coeff(:, 1))
            colorbar()
        end

    end

    % Correlate each feature column with each gene subset column to get a subset x feature matrix
    numFeatures = size(sortedFMatrix, 2);
    geneFeatures = zeros(numFeatures, numSubsets);
    geneFeaturesP = geneFeatures;

    xlabels = {
    "", % 1
    "A gene-expression proxy for serotonin receptor density",
    "", % 3
    "A gene-expression proxy for dopamine receptor density",
    "", % 5
    "A gene-expression proxy for glutamate receptor density",
    "", % 7
    "A gene-expression proxy for acetylcholine receptor density",
    "A gene-expression proxy for GABA receptor density", % 9
    "",
    "A gene-expression proxy for histamine receptor density", % 11
    "",
    "Noradrenaline (metabolism)", % 13
    "A gene-expression proxy for neurotransmitter receptor density",
    };

    for i = 1:numFeatures
        featureData = sortedFMatrix(:, i);

        for j = 1:numSubsets
            geneData = subsetMatrix(:, j);
            % Calculate correlation (Spearman)
            [c, p] = corr(featureData, geneData, 'Rows', 'pairwise', 'Type', 'Spearman');
            geneFeatures(i, j) = c;
            geneFeaturesP(i, j) = p;

            if (i == 1) & isPca & false % Plots scatters of PCA 1 vs Gene Expression
                p = figure();
                plot(geneData, featureData, '.', 'MarkerSize', 12)
                lsline
                ylabel(sprintf("PCA %i", i))
                xlabel(xlabels{subsetTitles(j)})
                title(sprintf("Timeseries PCA %i vs %s Gene Expressions (r=%.4f)", i, titles2{subsetTitles(j)}, c))
                saveas(p, sprintf('final_data/plots/pca1_%s.png', titles2{subsetTitles(j)}))
            end

        end

    end

    % Get top features for each subset by p value in a table/cell array

    topFeatures = {};

    for i = 1:numSubsets
        topFeatures{i} = operations;
        topFeatures{i}.Corr = geneFeatures(:, i);
        topFeatures{i}.Abs_Corr = abs(topFeatures{i}.Corr);
        topFeatures{i}.P_Value = geneFeaturesP(:, i);
        x = topFeatures{i};
        topFeatures{i} = x(~isnan(x.Corr), :);
        topFeatures{i}.Q_Value = mafdr(topFeatures{i}.P_Value);
        topFeatures{i} = sortrows(topFeatures{i}, {'Abs_Corr'}, 'descend');
    end

    featuresTableMeta = table(acronymsTitles, topFeatures');
    featuresTableMeta.Properties.VariableNames = {'Gene_Subset', 'Features'};

    % Plot histogram

    colors = distinguishable_colors(numSubsets, 'w');
    colors2 = {};

    for i = 1:size(colors, 1)
        colors2{i} = colors(i, :);
    end

    % Violin Plot

    figure()

    cellGF = {};

    for i = 1:numSubsets
        cellGF{i} = abs(geneFeatures(:, i));
    end

    BF_JitteredParallelScatter(cellGF, true, true, true, struct('theColors', {colors2}));
    title(sprintf("Correlations of Time Series Features with Gene Expression", numFeatures))
    ylabel("Magnitude of Spearman correlation coefficient, |\rho|")
    xticklabels(horzcat({''}, titles2{subsetTitles}));
    xlabel("Neuroreceptor Genes")

    % Plot scatter of each gene subset
    figure()
    [~, ax] = plotmatrix(subsetMatrix);

    for i = 1:numSubsets
        index = subsetTitles(i);

        if mod(i, 2) == 1 | true
            ax(i, 1).YLabel.String = titles2{index};
        end

        ax(numSubsets, i).XLabel.String = titles2{index};
    end
    title("Gene Expressions across Isocortex Areas, Faceted by Gene Subsets (PCA)")

end
