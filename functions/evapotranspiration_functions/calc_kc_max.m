function [ Kc_max ] = calc_kc_max(WindSpeed, RHmin, PlantHeight, Kcb )
    % Calculates Soil evaporation reduction coefficient.
    % Used in van der Veldt, R. (2017). Challenges of modelling soaring flight in humid landscapes
    % Edited by Rens van der Veldt, original by Bart Sweerts

    %% [Input]
    
    % WindSpeed = speed of wind [m s-1]
    % RHmin = min relative humidity []
    % PlantHeight = height of vegetation [m]
    % Kcb = basal crop factor []
    
    %% [Output]
    
    % Kc_max = maximum crop factor value based on limits of evaporation in
        % a wetted area
    

    %% Calculate Kc_max
%     [row,col] = size(Kcb);
    a = 0.04*(WindSpeed-2);
    b = 0.004*(RHmin-45).*(PlantHeight/3)^0.3;
    
    Kc_max = max((1.2+(a-b)), Kcb+0.05);
    
%     for r=1:row
%         for c=1:col
%             Kc_max = max((1.2+(a(r,c)-b(r,c))), Kcb(r,c)+0.05);
%         end
%     end
    
    
end

