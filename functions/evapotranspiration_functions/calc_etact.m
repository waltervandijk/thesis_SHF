function [ ETact ] = calc_etact( ET0, Kcb, Ks, Ke)
    % Calculates actual evapotranspiration and evaporation.
    % Used in van der Veldt, R. (2017). Challenges of modelling soaring flight in humid landscapes
    % Edited by Rens van der Veldt, original by Bart Sweerts
    % Edited by Walter van Dijk 2020

    %% [Input]
    
    % Et0 = Reference evapotranspiration
    % Kcb = Vegetation part of crop factor
    % Ks = Water stress coefficient
    % Ke = Soil evaporation coefficient
    
    %% [Output]
    
    % ETact = actual evapotranspiration.    
        
    %% Calculate actual evapotranspiration
    
    ETact = ET0.*(Ks .* Kcb + Ke);

end

