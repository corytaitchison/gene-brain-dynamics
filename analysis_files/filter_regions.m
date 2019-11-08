function [cSortedG, cSortedF, structInfo] = filter_regions(sortedG, sortedF, msi, regions)
    % Only retain rows that match the region that given by {'Region'}

    masks = zeros(size(sortedG, 1), 1);

    for i = 1:length(regions)
        mask = contains(msi.divisionLabel, regions{i}, 'Ignorecase', true);
        masks = mask + masks;
    end

    masks = logical(masks);
    cSortedG = sortedG(masks, :);
    cSortedF = sortedF(masks, :);
    structInfo = msi(masks, :);
end