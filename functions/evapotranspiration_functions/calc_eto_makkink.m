function [ lET0 ] = calc_eto_makkink(Rs, Ta, TimeStep)
    % Calculates Makkink reference evapotranspiration.
    % Based on de Bruin & LaBlans (1998)
    % Used in van Dijk, W. (2020)
    
    
    %% [Input]
    % variables:
    % Rs = incoming shortwave radiation                 [J m-2]
    % Ta = temperature                                  [deg C]
    
    % constants:
    C1 = 0.65; % C1 =  constante (De Bruin (1987): 0.65)           []
    ga = 0.000665 * 101.325; % ga = psychrometric constant (estimate 
                             %          avgerage P) [kPaC-1]    
    LatentHeatEvap = 2.45e+6;               % latent heat of evaporation of water [J kg-1]
	DensWater = 1000; %[kg m-3]
    
    %% [Output]
    % ET0 = reference evapotranspiration.               [kPa]
        
    %% Calculate reference evapotranspiration
    a = 6.1078; % [mbar]
    b = 17.294; % 
    c = 237.73; % [deg C]

    s = (a*b*c)./((c+Ta).^2) .* exp((b*Ta)./(c+Ta)); %slope of e at Ta
    
    ET0 = C1 *(s./(s+ga)) .* (Rs./TimeStep);      % [kg m-2 s-1] 
    lET0 = ET0 / (LatentHeatEvap.*DensWater) ; % [W m-2]
end

