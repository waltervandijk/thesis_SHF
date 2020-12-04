function [ CoverFraction ] = calc_coverfrac( LAI, SolarAngle )
    % Calculates Cover Fraction.
    % Used in van der Veldt, R. (2017). Challenges of modelling soaring flight in humid landscapes
    % Edited by Rens van der Veldt, original by Casper Borgman
    % Edited by Walter van Dijk 2020
    
    %% Input
    
    % LAI = Leaf area index.
    % Solar_Zenith = Solar zenith.
    
    %% Output
    
    % CoverFraction = Light cover fraction.
    
    %% Calculate Cover fraction
    % Using vegetation cover fraction to calculate fraction of light
    % reaching the surface and the fraction reaching vegetation. The
    % fraction is based on the Leaf Area Index. LAI changes over time
    
    % solar zenith indicates angle of sun: range for NL = -0.8:0.8
    if SolarAngle > 0 % if is daytime
        CoverFraction = 1 - exp(-0.5 * LAI./sind(SolarAngle));  % -0.5: Choudhury et al. [1987]

    else % is night so all is 'covered'
        CoverFraction = ones(size(LAI));
    end
        
end

