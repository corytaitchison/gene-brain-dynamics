function final_cluster(fdMatrix, fStructInfo, subsetTitles, nanTolerance, regions)
    %  Usage:
    %

    % Load gene data
    [gd, geneInfo, geneStructInfo] = LoadGeneExpressionData();

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

    % Load in acronyms as matrix [1,2,4,5,6,7,8,9,10,11,12,13,14]
    titles = {
    "Serotonin (metabolism)", % 1
    "Serotonin (receptor activity)",
    "Dopamine (metabolism)", % 3
    "Dopamine (neurotransmitter receptor activity)",
    "Glutamate (metabolism)", % 5
    "Glutamate (receptor activity)",
    "Acetylcholine (metabolism)", % 7
    "Acetylcholine (receptor activity)",
    "GABA (receptor activity)", % 9
    "GABA (receptor binding)",
    "Histamine (metabolism)", % 11
    "Histamine (receptor activity)",
    "Noradrenaline (metabolism)", % 13
    "Neurotransmitter (neurotransmitter metabolic process)",
    };
    % subsetTitles = [2, 4, 6, 8, 9, 12, 14];
    [acronymsList, acronymsTitles] = create_acronyms(subsetTitles);

    % For each set of genes in acronyms, calculate PCA and add to matrix region x gene subset
    numSubsets = size(acronymsList, 2);
    subsetMatrix = [];

    figure()
    hold on;
    zscore_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x, 'omitnan')), std(x, 'omitnan'));

    colors = distinguishable_colors(numSubsets, 'w');

    for i = 1:numSubsets
        [~, ~, geneMask] = intersect(acronymsList{i}, geneInfo.acronym);
        subG = sortedG(:, geneMask);
        % Calculate zscore and pca excluding NaNs
        % coeff = pca(zscore_xnan(subG), 'Rows', 'pairwise');
        coeff = tsne(zscore_xnan(subG)', 'Algorithm', 'barneshut', 'Distance', 'euclidean', 'NumDimensions', 2);
        % Plot
        x = coeff(:, 2);
        y = coeff(:, 1);
        plot(x, y, 'x', 'color', colors(i, :))
    end

    hold off;
    legend(acronymsTitles)
    xlabel("TSNE 1")
    ylabel("TSNE 2")
    title("PCA Clustering of Genes Expression for Each Subset")

end
