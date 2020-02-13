function [ FracWater, Rzone_Capacity, Rzone_Porosity ] = calc_soilmoisture( rootzone_sat, rootzone_satdef, rootzone )
    % Calculates Soil moisture.
    % Used in van der Veldt, R. (2017). Challenges of modelling soaring flight in humid landscapes
    % Edited by Rens van der Veldt, original by Casper Borgman
    
    %% [Input]
    
    % rootzone_sat = Rootzone saturation.
    % rootzone_satdef = Saturation deficiency in rootzone.
    % rootzone = Rootzone depth.
    
    %% [Output]
    
    % FracWater = Fraction of water in the soil.
    % Rzone_Capacity = Capacity of the rootzone.
    % Rzone_Porosity = Porosity of the rootzone.
    
    %% Calculate Fraction of water, rootzone capacity/-porosity.
    for i= 1:size(rootzone_sat,3)
        FracWater(:,:,i)= rootzone_sat(:,:,i) ./ rootzone;                    
    end
    
    Rzone_Capacity = rootzone_sat + rootzone_satdef;                 
    Rzone_Porosity = Rzone_Capacity(:,:,1) ./ rootzone;
end

