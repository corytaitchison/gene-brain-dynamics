function [acronymMat, subsetTitles] = create_acronyms(indices)
    % Creates a table of acronyms for each gene subset

    titles = [
        "Serotonin (metabolism)"
        "Serotonin (receptor activity)"
        "Dopamine (metabolism)"
        "Dopamine (neurotransmitter receptor activity)"
        "Glutamate (metabolism)"
        "Glutamate (receptor activity)"
        "Acetylcholine (metabolism)"
        "Acetylcholine (receptor activity)"
        "GABA (receptor activity)"
        "GABA (receptor binding)"
        "Histamine (metabolism)"
        "Histamine (receptor activity)"
        "Noradrenaline (metabolism)"
        "Neurotransmitter (neurotransmitter metabolic process)"
        ];

    subsetTitles = titles(indices);

    acronymMat = {};

    for i = 1:length(indices)
        index = indices(i);
        acronyms = load_data(index, titles(index));
        acronymMat{i} = unique(acronyms);
    end

end

function acronyms = load_data(index, title)
    file = sprintf("final_data/%s.txt", title);
    fid = fopen(file);

    A = textscan(fid, "%s%s%s%s%s%s%s%s%s%s%*[^\n]", 'HeaderLines', 1);
    acronyms = A{2};

    fclose(fid);
    clear fid A;

end