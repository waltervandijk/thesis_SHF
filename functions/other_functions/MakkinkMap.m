function [CropFactorMap] = MakkinkMap(NHdata,daynum,factsvat_split,LGN2SIMGROtable)
    % MakkinkMap makes a map by assigning land use classes to the respective
    % crop factors provided in input table
    %
    % Walter van Dijk 2020
    %% INPUTS
    % NHdata = land use map
    % daynum = day of year
    % fact_split = table containing crop factor values for land use and day
    % LGN2SIMGROtable = table to convert LGN classes (from the land use
    % map) to SIMGRO classes (from the crop factor table)

    CropFactorMap=double(NHdata);
    LGN2SIMGROtable=table2array(LGN2SIMGROtable);

    for idx=1:43 % loop through the 43 unique categories of lgn2018nh
        CropFactorMap(CropFactorMap==LGN2SIMGROtable(idx,1))= ...
            factsvat_split.(num2str(LGN2SIMGROtable(idx,2)))(daynum);
    end
end