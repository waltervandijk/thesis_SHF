function [DataOut] = interpERAgrid(ERA5longrid,ERA5latgrid,ERA5data,LGNlongrid,LGNlatgrid)
    %% interpolation of ERA5 data
    % 2D interpolation of gridded data.
    % This function automises the tranposing necessary to use the fast interpn
    % function.
    %
    % Walter van Dijk 2020

    %% gathering informaion
    % xextent=1:size(LGNlatgrid,2);
    % yextent=1:size(LGNlongrid,1);
    % numpix=xextent(end)*yextent(end);

    %% interpolation
    interpdata=ERA5data(:,:,1)';
    interplon=ERA5longrid';
    interplat=ERA5latgrid';

    DataOut=interpn(fliplr(interplon),fliplr(interplat),fliplr(interpdata),LGNlongrid,LGNlatgrid);
end