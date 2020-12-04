function [DataOutMatrix] = interpScattered(DataLon,DataLat,DataVal,QLon,QLat,Method)
    %% interpolation of scattered data
    % Used for matrices that are not in meshgrid format
    %
    % Walter van Dijk 2020

    %% Inputs
    % DataLon = lon of data
    % DataLat = lat of data
    % DataVal = data values
    % QLon = lon of query points
    % QLat = lat of query points
    % Method = interpolation method, can be linear or nearest

    %% gathering informaion
    xextent=1:size(QLat,2);
    yextent=1:size(QLon,1);

    %% interpolation
    interpdata=DataVal(:);
    interplon=DataLon(:);
    interplat=DataLat(:);

    F=scatteredInterpolant(interplon,interplat,interpdata,Method);
    DataOut=F(QLon(:),QLat(:));

    DataOutMatrix=zeros(size(yextent,2),size(xextent,2));
    DataOutMatrix(:)=DataOut;
end