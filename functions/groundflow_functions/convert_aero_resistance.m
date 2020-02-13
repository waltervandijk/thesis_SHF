function [ Rh ] = convert_aero_resistance( lgn2018, aa, bb, cc )
    % Converts lgn2018 map to aerodynamic resistance map.
    % Used in van Dijk, W. (2020)
    % Edited by Rens van der Veldt, original by Casper Borgman
    % Edit Walter van Dijk
    
    %% [Input]
    
    % lgn2018 = lgn2018 landuse map to convert to aerodynamic resistance
    % aa = [s/m] Forest Aerodynamic resistance      
    % bb = [s/m] Crop aerodynamic resistance
    % cc = [s/m] Grassland + Bare soil aerodynamic resistance
    
    %% [Output]
    
    % Rh = A map of aerodynamic resistances.
    
        
    %% Convert clusters.
    Rh = lgn2018;
       
    Rh(Rh == 1) = cc;       % grass
    Rh(Rh == 2) = bb;       % maize
    Rh(Rh == 3) = bb;       % potato
    Rh(Rh == 4) = bb;       % beets
    Rh(Rh == 5) = bb;       % cereals
    Rh(Rh == 6) = bb;       % other crops
    Rh(Rh == 8) = bb;       % greenhouses
    Rh(Rh == 9) = aa;       % boomgaarden
    Rh(Rh == 10) = bb;      % flower bulbs
    Rh(Rh == 11) = aa;      % deciduous forest
    Rh(Rh == 12) = aa;      % conifer forest
    Rh(Rh == 16) = cc;       % salt water
    Rh(Rh == 17) = cc;       % fresh water
    Rh(Rh == 18) = cc;      % buildings
    Rh(Rh == 19) = cc;      % buildings
    Rh(Rh == 20) = aa;      % forest in build up area
    Rh(Rh == 22) = aa;      % forest in build up area
    Rh(Rh == 23) = cc;      % grass in build area
    Rh(Rh == 24) = cc;      % empty ground area
    Rh(Rh == 25) = cc;      % roads
    Rh(Rh == 26) = cc;      % rural buildings
    Rh(Rh == 27) = cc;      % other rural
    Rh(Rh == 28) = cc;      % grass rural
    Rh(Rh == 30) = cc;      % kwelder
    Rh(Rh == 31) = cc;      % open sand
    Rh(Rh == 32) = cc;      % duneveg
    Rh(Rh == 33) = cc;      % duneveg
    Rh(Rh == 34) = cc;      % duneveg
    Rh(Rh == 35) = cc;      % open sand
    Rh(Rh == 36) = cc;      % heath
    Rh(Rh == 37) = cc;      % grassy heath
    Rh(Rh == 38) = cc;      % very grassy heath
    Rh(Rh == 39) = cc;      % peatland
    Rh(Rh == 40) = aa;      % forest in peat
    Rh(Rh == 41) = bb;      % other swampy veg
    Rh(Rh == 42) = bb;      % reed
    Rh(Rh == 43) = aa;      % forest in swamp
    Rh(Rh == 45) = cc;      % natural grasslands
    Rh(Rh == 46) = cc;      % coastal grasses
    Rh(Rh == 47) = cc;      % other grasses
    Rh(Rh == 61) = aa;      % tree farm
    Rh(Rh == 62) = aa;      % orchards
    Rh(Rh == 321) = bb;       % low shrub in peat
    Rh(Rh == 322) = bb;       % low shrub in marsh
    Rh(Rh == 323) = bb;       % low shrub other
    Rh(Rh == 331) = aa;       % high shrub in peat
    Rh(Rh == 332) = aa;       % high shrub in marsh
    Rh(Rh == 333) = aa;       % high shrub other
end