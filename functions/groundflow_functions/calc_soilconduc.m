function [ Ksoil ] = calc_soilconduc( Saturation, Kdry, Porosity, K_Water )

    % Calculates Soil conductivity.
    % Used in van der Veldt, R. (2017). Challenges of modelling soaring flight in humid landscapes
    % Edited by Rens van der Veldt, original by Casper Borgman
    % Edited by Walter van Dijk 2020
    % Based on Johansen (1977)
    
    
    %% Input
    
    % PDSoil = Particle bulk density.
    % Saturation = Soil saturation %                            [0 - 1]
    % Porosity = Porosity %                                     [0 - 1]
    % K_Water = Water thermal conductivity                      [W/m/K]
    % Kdry = Dry soil thermal conductivity                      [W/m/K]
    
    %% Output
    
    % Ksoil = Soil thermal conductivity.                        [W/m/K]

    %% Calculate thermal conductivities
    
    % Calculate the Kersten Number (depends on Saturation)
    Ke = 1 +  (log10(Saturation));                                         

    % Calculate saturated conductivity
    Ksat = Kdry + K_Water * Porosity;  
    
    % Calculate Soil conductivity for each cell.
    Ksoil= ((Ksat - Kdry) .* Ke) + Kdry;
    Ksoil(Saturation<0.1)=Kdry(Saturation<0.1);
    
end

