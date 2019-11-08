function final_cluster(fdMatrix, fStructInfo, subsetTitles, nanTolerance, regions)
    %  Plots a t-SNE clustering of gene expression across the regions

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

    % Load in acronyms as matrix
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
    [acronymsList, acronymsTitles] = create_acronyms(subsetTitles);

    % For each set of genes in acronyms, calculate t-SNE and add to matrix region x gene subset
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
        coeff = tsne(zscore_xnan(subG)', 'Algorithm', 'barneshut', 'Distance', 'euclidean', 'NumDimensions', 2);
        % Plot
        x = coeff(:, 2);
        y = coeff(:, 1);
        plot(x, y, 'x', 'color', colors(i, :))
    end

    hold off;
    legend(acronymsTitles)
    xlabel("t-SNE 1")
    ylabel("t-SNE 2")
    title("t-SNE Clustering of Genes Expression for Each Subset")

end
