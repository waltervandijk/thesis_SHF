function [ Kr ] = calc_kr( TEW, REW, sat, rootzone)
    % Calculates Soil evaporation reduction coefficient.
    % Used in van der Veldt, R. (2017). Challenges of modelling soaring flight in humid landscapes
    % Edited by Rens van der Veldt, original by Bart Sweerts
    % Edited by Walter van Dijk 2020

    %% [Input]
    
    % TEW = Total evaporable water. [mm]
    % REW = readily evaporable water. [mm]
    % ssd = soil saturation deficiency [%]
    % rootzone = rootzone depth [mm]
    
    %% [Output]
    
    % Kr = Soil evaporation reducing coefficients.
    % De = Soil moisture depleted
    % porvol = Pore volume;

    %% Calculate Kr and convert to proper ratios.
    De = (1-sat).*rootzone;
    Kr = (TEW-De)./(TEW-REW);
    Kr(TEW-REW<=0)=0; % REW is bigger than TEW on roads: Kr should be 0 on roads
    Kr(Kr>1)=1;
    Kr(Kr<0)=0;
    
end

