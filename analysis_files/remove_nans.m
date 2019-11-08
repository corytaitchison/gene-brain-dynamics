function [sortedG, sortedTS, geneInfo, structInfo] = remove_nans (gd, tsd, nans, threshold, oldGeneInfo, oldStructInfo, doPrint)
    % Usage: [sortedG, sortedTS, geneInfo, mergedStructInfo] = remove_nans(sortedG, sortedTS, nans, 0.2, geneInfo, mergedStructInfo)
    % Removes nans from rows and columns above the given threshold 

    if nargin < 7 || isempty(doPrint)
        doPrint = true;
    end 
    
    [numAreas, numGenes] = size(gd);
    aggNans = zeros(1, numGenes);
    for gene = 1:numGenes
        aggNans(gene) = mean(nans(:,gene));
    end
    mask = aggNans <= threshold;
    
    aggNansA = zeros(1, numAreas);
    for area = 1:numAreas
        aggNansA(area) = mean(nans(area,:));
    end
    maskA = aggNansA <= threshold;

    if doPrint
        fprintf("%i genes contain more than %i%% NaNs\n", numGenes-sum(mask), threshold*100);
        fprintf("%i areas contain more than %i%% NaNs\n", numAreas-sum(maskA), threshold*100);
        fprintf("Areas removed: ");
        oldStructInfo(~maskA,1)
    end
    
    sortedG = gd(maskA, mask);
    sortedTS = tsd(maskA, :);
    geneInfo = oldGeneInfo(mask, :);
    structInfo = oldStructInfo(maskA, :);
end