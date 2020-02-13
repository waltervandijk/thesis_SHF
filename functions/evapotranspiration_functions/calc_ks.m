function [ Ks ] = calc_ks( S01, Ssd01 )
    % Calculates Soil water stress coefficients.
    % Used in van der Veldt, R. (2017). Challenges of modelling soaring flight in humid landscapes
    % Edited by Rens van der Veldt, original by Bart Sweerts

    %% [Input]
    
    % S01 = Soil strorage water.
    % Fpwp = Wilting point values.
    % Ssd01 = Soil storage water deficit.
    
    %% [Output]
    
    % Ks = Soil water stress coefficients.

    %% Calculate Ks and convert to a set of ratios.
    TAW = S01+Ssd01;
    Dr = Ssd01;
    p = 0.50; % amount that can be extracted without stress
    %Ks = ratio_filter((S01(:,:))./((1-Fpwp) .* (S01 + Ssd01)));

    Ks = (TAW-Dr)./((1-p)*TAW);
    Ks(Ks>1)=1;
end

