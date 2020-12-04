function [ albedo ] = calc_albedo( NHdata,Saturation,aSoilDry,aSoilWet,aRoad,aUrban,aWater)
    % Calculates soil albedo
    % Walter van Dijk 2020
    
    %% Input
    
    % aSoilDry = dry soil albedo
    % aSoilWet = wet soil albedo
    % aRoad = road albedo
    % aUrban = urban area albedo
    
    %% Output
    
    % albedo = albedo for every pixel
    
    %% Calculate albedo
    % soil albedo is dependant on soil saturation, but can't be lower than
    % dry albedo
    albedo = Saturation.*(aSoilWet-aSoilDry)+aSoilDry;
        
    % overlay of some categories
    albedo(NHdata==16 | NHdata==17) = aWater;
    albedo(NHdata==25) = aRoad;
    albedo(NHdata==18 | NHdata==19 | NHdata==26) = aUrban;
    
end

