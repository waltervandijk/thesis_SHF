function [CropFactorMap] = MakkinkMap(NHdata,daynum,factsvat_split,LGN2SIMGROtable)
% MakkinkMap makes a map of crop factors provided the land use and
% conversion array


CropFactorMap=double(NHdata);
LGN2SIMGROtable=table2array(LGN2SIMGROtable);

for idx=1:43 % loop through the 43 unique categories of lgn2018nh
    CropFactorMap(CropFactorMap==LGN2SIMGROtable(idx,1))= ...
        factsvat_split.(num2str(LGN2SIMGROtable(idx,2)))(daynum);
end