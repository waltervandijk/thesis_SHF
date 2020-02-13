function [ Ksoil ] = calc_soilconduc( Saturation, Kdry, Porosity, K_Water )

    % Calculates Soil conductivity.
    % Used in van der Veldt, R. (2017). Challenges of modelling soaring flight in humid landscapes
    % Edited by Rens van der Veldt, original by Casper Borgman
    % Using Kersten number based on 
    
    
    %% Input
    
    % PDSoil = Particle bulk density.
    % Saturation = Soil saturation.
    % K_Quartz = Quartz thermal conductivity.
    % FracQuartz = Quartz texture fraction.
    % K_Mineral = Mineral thermal conductivity.
    % FracMineral = Mineral texture fraction.
    % FracUrban = Urban texture fraction.
    % K_Urban = Urban thermal conductivity.
    % FracOrganicMatter = OM texture fraction.
    % K_OrganicMatter = OM thermal conductivity.
    % Porosity = Porosity.
    % K_Water = Water thermal conductivity.
    
    %% Output
    
    % Kdry = Dry soil thermal conductivity.                  % [W/m/K]
    % Ksolid = Solid soil thermal conductivity.              % [W/m/K]
    % Ksat = Saturated soil thermal conductivity.            % [W/m/K]
    % Ksoil = Soil total thermal conductivity.               % [W/m/K]

    %% Calculate thermal conductivities
    
%     % Calculate the dry conductivity
%     Kdry = ((0.137 * PDSoil) + 0.0647) ./ (2.7 - (0.947 * PDSoil));
%     Kdry(Kdry == 0) = NaN;
    
    % Calculate the Kersten Number (depends on Saturation)
    Ke = 1 +  (log10(Saturation));                                         

    % Calculate saturated conductivity
    Ksat = Kdry + K_Water * Porosity;  
    
    % Calculate Soil conductivity for each cell.
%Kdry =1;
    Ksoil= ((Ksat - Kdry) .* Ke) + Kdry;
    Ksoil(Saturation<0.1)=Kdry(Saturation<0.1);
    
%     Ksoil = zeros(size(Saturation));    
%     for i = 1:size(Saturation,1)
%         for j = 1:size(Saturation,2)
%             
%             % If there is saturation under 0.1, use dry conductivity.
%             if Saturation(i,j) > 0.1
%                 Ksoil(i,j)= ((Ksat(i,j) - Ksolid(i,j)) * Ke(i,j)) + Ksolid(i,j);
%             else
%                 Ksoil(i,j) = Ksolid(i,j);
%             end
%         end
%     end
    
    % Convert conductivity to make it per hour [from W/m/K to J/m/K]
    %Ksoil(Ksoil < 0.1) = NaN;         
    %Ksoil = Ksoil * 3600; 
end

