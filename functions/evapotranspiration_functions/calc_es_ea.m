function [Es,Ea] = calc_es_ea(T2m, D2m)
    % Calculates saturation and actual vapour pressure.
    % Used in van der Veldt, R. (2017). Challenges of modelling soaring flight in humid landscapes
    % Edited by Rens van der Veldt, original by Bart Sweerts
  
    %% [Input]
    
    % T2m = temperature at 2m.              [oC]
    % T2m = dewpoint temperature at 2m.     [oC]
    
    %% [Output]
    
    % Es = saturation vapour pressure.      [kPa]
    % Ea = actual vapour pressure.          [kPa]
    
    %% Calculate saturation and actual vapour pressure.
    Es = 0.6108.*(exp(17.27 .* T2m ./(T2m + 237.3)));                                   
    Ea = 0.6108.*(exp(17.27 .* D2m ./(D2m + 237.3)));
end