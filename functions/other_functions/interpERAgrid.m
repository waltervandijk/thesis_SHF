function [DataOut] = interpERAgrid(ERA5longrid,ERA5latgrid,ERA5data,LGNlongrid,LGNlatgrid)
%% interpolation of ERA5 data
% Walter van Dijk
% 11/2019
    % 6-7 sec with griddedInterpolant
    % 4-5 sec with interpn
    % 1-2 sec in parallel
    % interpolation is using only one core, while 4 are available
    
%% gathering informaion
% xextent=1:size(LGNlatgrid,2);
% yextent=1:size(LGNlongrid,1);
% numpix=xextent(end)*yextent(end);

%% interpolation
interpdata=ERA5data(:,:,1)';
interplon=ERA5longrid';
interplat=ERA5latgrid';

% G=griddedInterpolant(fliplr(interplon),fliplr(interplat),fliplr(interpdata)); %faster than scatteredInterpolant
% DataOut=G(LGNlongrid,LGNlatgrid);

DataOut=interpn(fliplr(interplon),fliplr(interplat),fliplr(interpdata),LGNlongrid,LGNlatgrid);