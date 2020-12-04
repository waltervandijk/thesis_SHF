function [ Ks ] = calc_ks( S01, Ssd01 )
    % Calculates Soil water stress coefficients.
    % Based on Allen et al. (1998) FAO
    % Used in van der Veldt, R. (2017). Challenges of modelling soaring flight in humid landscapes
    % Edited by Rens van der Veldt, original by Bart Sweerts
    % Edited by Walter van Dijk 2020

    %% [Input]
    
    % S01 = Soil water saturation
    % Ssd01 = Soil water saturation deficit
    
    %% [Output]
    
    % Ks = Soil water stress coefficients.

    %% Calculate Ks and convert to a set of ratios.
    TAW = S01+Ssd01;
    Dr = Ssd01;
    p = 0.50; % amount that can be extracted without stress

    Ks = (TAW-Dr)./((1-p)*TAW);
    Ks(Ks>1)=1;
end

