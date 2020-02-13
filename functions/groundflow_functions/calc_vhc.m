function [ VHC ] = calc_vhc( VHC_Mineral, VHC_Quartz, VHC_OrganicMatter,...
    VHC_Urban, VHC_Water, FracMineral, FracQuartz, FracOrganicMatter, FracUrban, FracPorosity, Saturation)
    % Calculates volumetric heat capacity.
    % Used in van der Veldt, R. (2017). Challenges of modelling soaring flight in humid landscapes
    % Edited by Rens van der Veldt, original by Casper Borgman
    
    %% [Input]
    
    % VHC_Mineral = Mineral volumetric heat capacity.
    % VHC_Quartz = Quartz volumetric heat capacity. 
    % VHC_OrganicMatter = OM volumetric heat capacity.
    % VHC_Urban = Urban volumetric heat capacity. 
    % VHC_Water = Water volumetric heat capacity. 
    % FracMineral = Mineral texture fractions. 
    % FracQuartz = Quartz texture fractions. 
    % FracOrganicMatter = OM texture fractions. 
    % FracUrban = Urban texture fractions. 
    % FracPorosity = Pore spaces texture fraction
    % Saturation = % of water in pore spaces
    
    %% [Output]
    
    % VHC = Volumetric heat capacity of the soil.       [J/m3/K]
    
    %% Calculate volumetric heat capacity.
    % Soil Volumetric Heat Capacity based on water content change           
    VHC = (VHC_Mineral * FracMineral) + (VHC_Quartz * FracQuartz) ...
    + (VHC_OrganicMatter * FracOrganicMatter) +( VHC_Urban *FracUrban)...
    + (VHC_Water * FracPorosity .* Saturation); 
    

end

