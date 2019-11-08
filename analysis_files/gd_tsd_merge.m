function [gd, tsd, mergedStructInfo] = gd_tsd_merge(gd, structInfoG, tsd, structInfoTS)
    % Usage: [sortedG, sortedTS, mergedStructInfo] = gd_tsd_merge(gd, geneStructInfo, aggLowFreqs, timeStructInfo)
    % Merges the gene and time series data sets, ensuring rows (areas) are consistent between both 

    [~, ia, ib] = intersect(structInfoG.acronym, structInfoTS.acronym, 'stable');
    gd = gd(ia, :);
    structInfoG = structInfoG(ia, :);
    tsd = tsd(ib, :);
    structInfoTS = structInfoTS(ib, :);

    [~, ia] = sort(structInfoG.acronym);
    [~, ib] = sort(structInfoTS.acronym);
    gd = gd(ia, :);
    structInfoG = structInfoG(ia, :);
    tsd = tsd(ib, :);
    structInfoTS = structInfoTS(ib, :);
    mergedStructInfo = [structInfoG(:, [1, 3:5])];
end
