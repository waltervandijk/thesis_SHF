function [LAIout] = constructLAI(LAIsummary,landusedata,decadenum)
%% construct LAI map for lgn2018
% Uses previously created summary of LAI values per land use category.
% The used LAI values are derived from Copernicus 300m LAI dataset.
% 
% Walter van Dijk 2020
%
% meant for lgn2018 classes
% function works for grid and vector landuse data 

time = decadenum-6;
if time > 24
    time = 24;
    disp('end of LAI data reached')
    disp('using last available LAI data')
end

if time < 1
    time = 1;
    disp('data before LAI record')
    disp('using first available LAI data')
end

LAIout=double(landusedata);
for i=1:43
    LAIout(LAIout==LAIsummary(i,5,time))=LAIsummary(i,1,time);
end

LAIout(landusedata == 16) = 0; % LAI of 0 above water
LAIout(landusedata == 17) = 0; % LAI of 0 above water
LAIout(landusedata == 23) = 0.2; % LAI of 0.2 above buildings, original value was 0.51, but is too high due to edge effects
LAIout(landusedata == 25) = 0; % LAI of 0 above roads/rails


end