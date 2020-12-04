function [ Flow, GroundFlow ] = calc_flows( Flow, AirDensity, AirHeatCapacity, Temp,...
    Rh, Ksoil, ThickL, ThickSurfL, NrLayer, CellX, CellY)
    % Calculates Ground flows.
    % Used in van der Veldt, R. (2017). Challenges of modelling soaring flight in humid landscapes
    % Edited by Rens van der Veldt, original by Casper Borgman
    % Edited by Walter van Dijk 2020
    
    %% Input
    
    % Flow = The initialized flow.                              [J m-2]
    % AirDensity = Density of the air.                          [kg m-3]
    % AirHeatCapacity = Heat capacity of air.                   [J kg-1 K-1]
    % Temp = Every layers temperature of the landscape.         [K]
    % Rh = Aerodynamic resistance of the landscape.             [s m-1]
    % Ksoil = Soil heat conductivity.                           [W m-1 K-1]
    % ThickL = Layer thickness.                                 [m]
    % NrLayer = Number of layers.                               []
    
    %% Output
    
    % Flow = Soil flow.                                         [W m-2]
        % positive for upward movement
    % NFlow = Net soil flow.                                    [W m-2]
	% GroundFlow = Energy transport into the soil               [W m-2]     
    
    %% Calculate Flows and Net flow.
    % Calculate veg layer - surface layer flow, needed for G
    Flow(:,:,1) = (AirDensity * AirHeatCapacity * (Temp(:,:,2)...
        - Temp(:,:,1))) ./ Rh;  %[J m-3] / [s m-1] = [J m-2 s-1] = [W m-2] 
   
    % Calculate surface - first soil layer flow.
    Flow(:,:,2) = (Temp(:,:,3) - Temp(:,:,2)) .* Ksoil ./ ...
        (0.5 .* ThickSurfL + 0.5 .* ThickL) ./ (CellX * CellY) ; % [W m-2]           
    
    % Calculate other soil layers flows.
    for k = 3 : NrLayer-1 % flow out of system at bottom not includo
        Flow(:,:,k) = (Temp(:,:,k + 1) - Temp(:,:,k)) .* Ksoil ./ ...
            ThickL ./ (CellX * CellY); % [W m-2]
    end
    
    GroundFlow = Flow(:,:,2)-Flow(:,:,1);  % [W m-2]
end

