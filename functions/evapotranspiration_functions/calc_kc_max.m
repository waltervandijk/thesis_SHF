function [ Kc_max ] = calc_kc_max(WindSpeed, RHmin, PlantHeight, Kcb )
    % Calculates maximum total crop factor value
    % Based on Allen et al. (1998) FAO
    %
    % Walter van Dijk 2020

    %% [Input]
    
    % WindSpeed = speed of wind             [m s-1]
    % RHmin = min relative humidity         []
    % PlantHeight = height of vegetation    [m]
    % Kcb = basal crop factor               []
    
    %% [Output]
    
    % Kc_max = maximum crop factor value based on limits of evaporation in
        % a wetted area         []
    

    %% Calculate Kc_max
    a = 0.04*(WindSpeed-2);
    b = 0.004*(RHmin-45).*(PlantHeight/3).^0.3;
    
    Kc_max = max((1.2+(a-b)), Kcb+0.05);
    
    
end

