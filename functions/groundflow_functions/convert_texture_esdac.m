function [ FracMineral, FracQuartz, FracOrganicMatter, FracUrban, FracPorosity ] = ...
    convert_texture_esdac( Texture)
    % Convert texture to esdac values.
    % Used in van der Veldt, R. (2017). Challenges of modelling soaring flight in humid landscapes
    % Edited by Rens van der Veldt, original by Casper Borgman
    % Edited by Walter van Dijk 2020
    
    %% [Input]
    
    % Texture = Texture dataset.
    % MinFracs = ESDAC mineral fractions per class [0-5;9].
    % QuartzFracs = ESDAC quartz fractions per class [0-5;9].
    % UrbanFracs = ESDAC urban fractions per class [0-5;9].
    % OrganicMatterFracs = ESDAC organic matter fractions per class [0-5;9].
    % PorosityFracs = ESDAC porosity fractions per class [0-5;9].
    
    % ESDAC fractions per texture class (0-5;9)
        % Porosity for use in conductivity calculations 
        % (After van Wijk and de Vries (1963) & Rawls (1982)
    MinFracs =           [0,   0.12, 0.21, 0.35, 0.43, 0.50, 0.09];
    QuartzFracs =        [0,   0.45, 0.35, 0.17, 0.08, 0.02, 0.01];
    OrganicMatterFracs = [0,   0.03, 0.04, 0.03, 0.04, 0.03, 0.10];
    UrbanFracs =         [1,   0,    0,    0,    0,    0,    0   ];
    PorosityFracs =      [0,   0.4,  0.4,  0.45, 0.45, 0.45, 0.8 ];
    
    %% [Output]
    
    % FracMineral = Dataset with correct mineral fractions.
    % FracQuartz = Dataset with correct quartz fractions.
    % FracOrganicMatter = Dataset with correct OM fractions.
    % FracUrban = Dataset with correct urban fractions.
    % Porosity = Dataset with correct porosity fractions.
    
    %% Copy Textures
    FracMineral = Texture;
    FracQuartz= Texture;
    FracOrganicMatter= Texture;
    FracUrban = Texture;
    FracPorosity= Texture;
    
    %% Assign class fractions
    % Esdac texture description
    % -9999 | -128 = Waterbodies (NaN)
    % 0 = cities & impermeable surface
    % 1= sand >65% clay <18%
    % 2= medium 18-35% clay, >15% sand
    % 3= Medium fine, 35% clay, 15% sand, fine material
    % 4= Fine, clay between 35% and 60%
    % 5= very fine, >60% clay
    % 9= peat
    
    % Mineral fractions
    FracMineral (Texture == 0) = MinFracs(1);
    FracMineral (Texture == 1) = MinFracs(2);
    FracMineral (Texture == 2) = MinFracs(3);
    FracMineral (Texture == 3) = MinFracs(4);
    FracMineral (Texture == 4) = MinFracs(5);
    FracMineral (Texture == 5) = MinFracs(6);
    FracMineral (Texture == 9) = MinFracs(7);
    FracMineral (Texture == -9999) = 0;
    FracMineral (Texture == -128) = 0;
    
    % Quartz fractions
    FracQuartz (Texture == 0) = QuartzFracs(1);
    FracQuartz (Texture == 1) = QuartzFracs(2);
    FracQuartz (Texture == 2) = QuartzFracs(3);
    FracQuartz (Texture == 3) = QuartzFracs(4);
    FracQuartz (Texture == 4) = QuartzFracs(5);
    FracQuartz (Texture == 5) = QuartzFracs(6);
    FracQuartz (Texture == 9) = QuartzFracs(7);
    FracQuartz (Texture == -9999) = 0;
    FracQuartz (Texture == -128) = 0;
    
    % OM fractions
    FracOrganicMatter (Texture == 0) = OrganicMatterFracs(1);
    FracOrganicMatter (Texture == 1) = OrganicMatterFracs(2);
    FracOrganicMatter (Texture == 2) = OrganicMatterFracs(3);
    FracOrganicMatter (Texture == 3) = OrganicMatterFracs(4);
    FracOrganicMatter (Texture == 4) = OrganicMatterFracs(5);
    FracOrganicMatter (Texture == 5) = OrganicMatterFracs(6);
    FracOrganicMatter (Texture == 9) = OrganicMatterFracs(7);
    FracOrganicMatter (Texture == -9999) = 0;
    FracOrganicMatter (Texture == -128) = 0;
    
    % Urban fractions
    FracUrban (Texture == 0) = UrbanFracs(1);
    FracUrban (Texture == 1) = UrbanFracs(2);
    FracUrban (Texture == 2) = UrbanFracs(3);
    FracUrban (Texture == 3) = UrbanFracs(4);
    FracUrban (Texture == 4) = UrbanFracs(5);
    FracUrban (Texture == 5) = UrbanFracs(6);
    FracUrban (Texture == 9) = UrbanFracs(7);
    FracUrban (Texture == -9999) = 0;
    FracUrban (Texture == -128) = 0;
    
    % Porosity fractions
    FracPorosity (Texture == 0) = PorosityFracs(1);
    FracPorosity (Texture == 1) = PorosityFracs(2);
    FracPorosity (Texture == 2) = PorosityFracs(3);
    FracPorosity (Texture == 3) = PorosityFracs(4);
    FracPorosity (Texture == 4) = PorosityFracs(5);
    FracPorosity (Texture == 5) = PorosityFracs(6);
    FracPorosity (Texture == 9) = PorosityFracs(7);
    FracPorosity (Texture == -9999) = 1;   
    FracPorosity (Texture == -128) = 1;   
    
end

