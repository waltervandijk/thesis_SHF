function [ ETact ] = calc_etact_evap( ET0, Kcb, Ks, Ke)
    % Calculates actual evapotranspiration and evaporation.
    % Used in van der Veldt, R. (2017). Challenges of modelling soaring flight in humid landscapes
    % Edited by Rens van der Veldt, original by Bart Sweerts

    %% [Input]
    
    % Et0 = Reference evapotranspiration.            
    % KcbMid = Crop coefficients.                            
    % Ks = Water stress coefficients. 
    % Ke = Soil evaporation coefficients.
    
    %% [Output]
    
    % ETact = actual evapotranspiration.    
        
    %% Calculate actual evapotranspiration
    
    ETact = ET0.*(Ks .* Kcb + Ke);

end

