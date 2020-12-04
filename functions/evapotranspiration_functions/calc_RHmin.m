function RHmin = calc_RHmin(T2m, D2m)
    % Calculates relative humidity
    % Walter van Dijk 2020
  
    %% [Input]
    
    % T2m = temperature at 2m.              [oK]
    % T2m = dewpoint temperature at 2m.     [oK]
    
    %% [Output]
    
    % Es = saturation vapour pressure.      [kPa]
    % Ea = actual vapour pressure.          [kPa]
    
    T2m = T2m - 273.15; % K to C
    D2m = D2m - 273.15; % K to C
    
    %% Calculate saturation and actual vapour pressure.
    Es = 0.6108.*(exp(17.27 .* T2m ./(T2m + 237.3)));                                   
    Ea = 0.6108.*(exp(17.27 .* D2m ./(D2m + 237.3)));
    
    RHmin = (Ea./Es)*100; %[percentage]
end