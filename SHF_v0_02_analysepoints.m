%% SHF Calculation script %%
% Used in: van Dijk, W. (2020). Modelling sensible heat flux over humid
% areas. 

% This script calculates sensible heat flux (H) as a residual of the 
% energy partitioning formula: Rn = G + ET + H. This is done on a 5x5m 
% and hourly resolution. Results are aggregated to 25x25m resolution 
% while storing statistic.

% Inputs should be available in workpath, matching the names used in the
% loading section. Outputs are stored as .mat files.

%% Computing requirements
% the total area consists of 16000*11000=176 mil cells
% RAM-test for 22.5 million cells on 16GB device:
    % Idle Ram (OS+Matlab2019): 3.0 GB
    % After initialisation:     8.5 GB
    % After daily calcs:        9.4 GB
    % During hourly calcs:      12 - 14.7 GB
    

    
%% Inputs                                  Description                                          Units

% LAI                                      Cover fraction per square unit of ground.            []
% rootzone                                 Rootzone depth                                       [m]
% NHdata                                   lgn2018 land use map                                 []
    % NHlat / NHlon                            lat/long information of lgn2018
% wp                                       wilting point                                        [%]
% ERA5 data                                temp 2m, dewpoint 2m, radiation, wind
% S01                                      Rootzone saturation                                  [%]
% Ssd01                                    Rootzone saturation defficiency                      [%]
% Texture                                  Texture classes of the landscape                     []
% SIMGRO crop factors                      Kcb - crop factors to be used with the land use map

%% Outputs

% SHF_output                                 Hourly sensible heat flux                            [W/m2]

%% INITIALISATION %%
    %% Clear workspace
    close all
    clc

    %% Add folders to path
    addpath('Output',genpath('functions'),genpath('data'))
    disp('path initialised')

    %% Control Constants
    Debugging = 1;           % put on zero for more efficient memory usage
    StartDay = 0;
    StartTime = StartDay*24+1;        % Time of start  [hour]
    EndTime = StartTime+4799;          % Time of end  [hour] (endtime-starttime should be n*24-1
        %StartTime = 1 = '01-Mar-2018 00:00:00';
        %EndTime = 24 = '01-Mar-2018 23:00:00';
    TimeStep = 3600;         % [sec/hour]
    CellSizeX = 5;           % [m]
    CellSizeY = 5;           % [m]
    ThickL = 0.5;           % Layer thickness   [m]
    ThickSoilSurfL = 0.05;       % Surface layer thickness   [m]
    ThickWaterSurfL = 1;   % Thickness water surface layer [m]
    ThickSeaSurfL = 2;     % Thickness sea mixing layer [m]
    NrLayer = 10;             % Number of layers. Veg = 1, Surface = 2, Layer1-5 use 3:7
    
%bulbs, grass, cereal crop, road, deciduous forest, lake
SubRow = [4140,4560,4680,13010,10750,10750];
SubCol = [3880,5200,5260,4800,1940,3350];
    %Subset = [12001,15000,3501,6500];
    %Subset = [4140,4140,3880,3880]; 
        % Some interesting areas %
            %Subset=[1,15000,1,11000];           %full size
            %Subset=[12601,13400,4601,5100];     %A10 Amsterdam
            %Subset=[7201,8000,2301,2800];       %foresty bit bergen
            %Subset=[8601,9400,1601,2100];       %bergen aan zee
            %Subset=[10001,10800,2901,3400];     %alkmaardermeer
            %Subset=[1,2000,2701,3900];          %below Texel
            %Subset=[5001,6400,7201,8200];       %greenhouses
            %Subset = [12000,15000,3501,6500];   %amsterdam

    %% Loading and adjusting input data
    if exist ('Saturation_int', 'var') && size(Saturation_int,1) == size(SubRow,1)
        disp('using preinitialised variables')
    else
        % loading data
        disp('loading data.')
        load('LAIsummary.mat'); % LAI stats used to make LAI maps from land use
        load('SimgroCropFactors_Makkink.mat','factsvat_split','LGN2SIMGROcat') % Crop factors from SIMGRO (Makkink)
        load('ERA5_comb.mat') % ERA5 data Mar-Oct 2018
            [ERA5latgrid,ERA5longrid]=ndgrid(double(ERA5lat_comb),double(ERA5lon_comb));
            ERA5u_v_comb=sqrt(ERA5u10_comb.^2+ERA5v10_comb.^2);
            ERA5datenum = datenum(1900,1,1,double(ERA5time_comb),0,0);
            ERA5datetime = datetime(datestr(ERA5datenum),'InputFormat','dd-MM-yyyy HH:mm:ss');
            clear('ERA5u10_comb','ERA5v10_comb','ERA5datenum','ERA5time_comb')
            clear('ERA5sshf_comb','ERAstl1','ERAstl2')%,'ERAstl4')%,'ERAstl3')
        load('Hydro.mat')
            S01=S01(251:600,151:600,:);             % reduce size to NH area
            Ssd01=Ssd01(251:600,151:600,:);
            S01_lat=S01_lat(251:600,151:600,:);
            S01_lon=S01_lon(251:600,151:600,:);
            S01(isnan(S01))=1;      % nan = water
            Ssd01(isnan(Ssd01))=0;  % nan = water
        disp('loading data..')
        load('rootzone.mat')                        %[cm]
        load('texture.mat'); % Interpolated soil texture map with urban areas and water
        texture_interp_int8=texture_interp_int8(SubRow,SubCol);
        texture_interp_int8=nansum(double(texture_interp_int8).*eye(size(texture_interp_int8)));
            UrbanArea = texture_interp_int8==0;
            wp=double(texture_interp_int8); % create wilting point array based on table
                for i =1:5 % different classes of texture
                    wp(texture_interp_int8==texture_wp(i,1))=texture_wp(i,2);
                end
                wp(wp==-128)=0;
            fc=double(texture_interp_int8); % create field capacity array based on table
                for i =1:5 % different classes of texture
                    fc(texture_interp_int8==texture_fc(i,1))=texture_fc(i,2);
                end
                fc(fc==-128)=1;
        disp('loading data...')
        load('lgn2018nh') % load lgn2018; contains NHlat, NHlon & NHdata
            NHdata=NHdata(SubRow,SubCol);
                NHdata=nansum(double(NHdata).*eye(size(NHdata)));
                NHdata=double(NHdata);
                NHdata(NHdata==65535)=NaN; %or 17 to make it open water
            NHlon=NHlon(SubRow,SubCol);
                NHlon=nansum(NHlon.*eye(size(NHlon)));
            NHlat=NHlat(SubRow,SubCol);
                NHlat=nansum(NHlat.*eye(size(NHlat)));

        pause(0.1)
        fprintf('loading data successful\n');

        %% Adjust variables 
        rootzone_mm = rootzone * 10;                  % Convert rootzone from cm to mm for calculations   [mm]
        clear('rootzone') %12mb
        
        %ERA5
        ERA5str_adj = ERA5cum2add (ERA5str_comb, ERA5datetime); %[J m-2] positive downwards
        ERA5ssr_adj = ERA5cum2add (ERA5ssr_comb, ERA5datetime); %[J m-2] positive downwards
        %ERA5sshf_adj = ERA5cum2add(ERA5sshf_comb,ERA5datetime);
        ERA5net_rad = ERA5ssr_adj+ERA5str_adj;    %[J m-2]
        ERA5d2m_adj = ERA5d2m_comb; % - 273.15;        %dewpoint 2m to [deg C]
        ERA5t2m_adj = ERA5t2m_comb; % - 273.15;        %temp 2m to [deg C]
        ERA5stl4_adj = ERAstl4; % - 273.15;        %temp 2m to [deg C]
        ERA5u_v_comb_2m = ERA5u_v_comb .* 0.748;      % Convert 10m wind to 2m wind  [m s-1]             [m s-1]
        clear('ERA5d2m_comb','ERA5t2m_comb','ERA5str_comb',...
            'ERA5u_v_comb','ERA5sshf_comb') %5x 70mb

        pause(0.1)
        fprintf('adjusted data\n');

        %% Initialise system parameters
        AirDensity = 1.225;              % Air density at sea level                      [Kg/m3] 
        AirHeatCapacity = 1013;          % Specific heat cap of air at constant pressure [J/kg/K] 
        LatentHeatEvap = 2.45e+6;        % latent heat of vaporisation of water          [J kg-1]
    	DensWater = 1000; %[kg m-3]
        
        % Surface layer (thicker for water)
        ThickSurfL = zeros(size(NHdata));
        ThickSurfL(:,:) = ThickSoilSurfL;
        ThickSurfL(NHdata==16) = ThickWaterSurfL;
        ThickSurfL(NHdata==17) = ThickSeaSurfL;
        
        % Aerodynamic resistances.
        aa = 200;            % Forest Aerodynamic resistance                 [s/m] 
        bb = 100;            % Crop aerodynamic resistance                   [s/m] 
        cc = 50;             % Grassland + Bare soil aerodynamic resistance  [s/m] 

        % Convert to specific textures.
        [ FracMineral, FracQuartz, FracOrganicMatter, FracUrban, FracPorosity ] = ...
            convert_texture_esdac(double(texture_interp_int8));

        % Volumetric heat capacity                           % [J/m3/K]
        VHC_Mineral = 2.01e6; VHC_Quartz = 2.01e6; VHC_OrganicMatter= 2.51e6;
        VHC_Water= 4.18e6; VHC_Air= 1.25e3; VHC_Urban= 2.06e6;
        % VHC_Urban assumed to be that of concrete

        % Thermal conductivity (After Hillel (1982), De Vries(1963), Alter(1969)    [W/m/K] 
        K_Mineral= 2.5; K_Quartz= 8; K_OrganicMatter= 0.25;
        K_Water= 0.57; K_Air= 0.025; K_Urban= 4.6;
        EddyConducMult = 500; % K_Water multiplier for water bodies 

        % Create Aerodynamic resistance map from CLC.
        Rh = convert_aero_resistance(NHdata, aa, bb, cc); % [s m-1]

        % Interpolate rootzone map
        rootzone_mm_int = interpScattered(rootzone_lon,rootzone_lat, rootzone_mm, NHlon, NHlat, 'linear');
        clear('rootzone_mm','rootzone_lat','rootzone_lon')

        % Determine total and readily evaporable water.
        TEW = (fc-0.5*wp) * 100; % [mm] Ze = 0.1m, depth of the surface soil layer that is subject to drying by evaporation 
        TEW(TEW<0)=0;
        REW = 10; % [mm] depth of water that can be evaporated without restriction
        clear('fc','wp')

        pause(0.1)
        fprintf('initialized system parameters\n');

        %% Initialise System Variables
        % Calculate soil moisture properties
        Saturation=S01./(S01+Ssd01);            % percentage
        Saturation_int = interpScattered(S01_lon,S01_lat,Saturation(:,:,1),NHlon,NHlat,'linear');
        Saturation_int(NHdata==16 | NHdata==17) = 1; % saturation = 100% at water
        Saturation_int(UrbanArea) = 0; % no soil water on impermeable surface

        % Calculate solid (dry) conductivity
        Ksolid = (K_Quartz * FracQuartz + K_Mineral * FracMineral + FracUrban ...
            * K_Urban + FracOrganicMatter * K_OrganicMatter); 
        
        % Initial Soil Volumetric Heat Capacity                                     [J/m3/K]
        VHC_Soil_Dry = VHC_Mineral * FracMineral + VHC_Quartz * FracQuartz ...
            + VHC_OrganicMatter * FracOrganicMatter + VHC_Urban * FracUrban;
        clear('VHC_Mineral','VHC_Quartz','VHC_OrganicMatter','VHC_Urban')
        clear('FracMineral','FracOrganicMatter','FracQuartz','FracUrban')
        VHC_Soil =  VHC_Soil_Dry + VHC_Water * FracPorosity .* Saturation_int; 

        % Create a matrix with initial temperatures, flows and heatcontent
        ERA5t2m_int = interp2(ERA5longrid,ERA5latgrid,ERA5t2m_adj(:,:,StartTime)',NHlon,NHlat); % interpolate to get the right size
        Temp = zeros([size(ERA5t2m_int,[1 2]),NrLayer]);                           % Create Temperature 3d matrix                          
        Temp(:,:,1) = ERA5t2m_int;                                                 % [oC] Initial vegetation top layer temperature                                             % [oC] Initial value soil temperature assumed to be 5 oC
        ERA5stl4_adj(isnan(ERA5stl4_adj))=min(min(ERA5stl4_adj(:,:,StartTime))); % initial sea temp = 5 deg C
        for f = 2:NrLayer
            Temp(:,:,f)= interp2(ERA5longrid,ERA5latgrid,ERA5stl4_adj(:,:,StartTime)',...
                NHlon,NHlat,"linear",0);  % [oC] Soil temp from ERA5
        end
        % Initialize flows and heat content.
        Flow = zeros([size(ERA5t2m_int,[1 2]),NrLayer]);                         % [] Create the flow 3d matrix                             
        HeatCont= zeros([size(ERA5t2m_int,[1 2]),NrLayer]);                        % [] Create heatcontent matrix
        HeatCont(:,:,2) = Temp(:,:,2).*VHC_Soil.*ThickSurfL;
        HeatCont(:,:,3:NrLayer) = Temp(:,:,3:NrLayer).* VHC_Soil * ThickL; % Create initial heatcontent for soil [J/m2]

        pause(0.1)
        fprintf('initialized system variables\n');
        
    end
        %% Initialise dynamic Variables
        disp(strcat('starttime : ',string(ERA5datetime(StartTime,:))))
        disp(strcat('endtime   : ',string(ERA5datetime(EndTime,:))))

        % Set startTime
        Time = StartTime;

        % save timekey
        TimeKey = table((StartTime:EndTime)',ERA5datetime(StartTime:EndTime,:));
        save([pwd,'\Output\TimeKey'],'TimeKey');
    
        % storevars
        StoreInd = 1;
        TempStore = zeros(NrLayer,EndTime-StartTime+1);
        TempWaterStore = zeros(NrLayer,EndTime-StartTime+1);
        EvapoStore = zeros(EndTime-StartTime+1,size(texture_interp_int8,2));
        shfStore = zeros(EndTime-StartTime+1,size(texture_interp_int8,2));
        CoverStorer = zeros(EndTime-StartTime+1,size(texture_interp_int8,2));
    
%% DYNAMIC CALCULATIONS %% 

while Time <= EndTime    
    %% Daily part
    pause(0.01)
    %% Determine the time
    daynum = day(ERA5datetime(Time),'dayofyear');
    disp(['starting calculations for day ',num2str(daynum),'...'])
    decadenum = ceil(daynum/10); %ceil(day/10)=quick&dirty, not exact
    
    %% Interpolate/create daily boundary conditions
    disp('updating daily boundary conditions...')
    tic
    % shape data into the same dimensions
    Saturation_int = interpScattered(S01_lon,S01_lat,Saturation(:,:,decadenum),NHlon,NHlat,'linear');  % 8sec for 1000x1000 cells
    Saturation_int(NHdata==16 | NHdata==17) = 1; % saturation = 100% at water
    Saturation_int(UrbanArea) = 0; % no soil water on impermeable surface
    S01_int = interpScattered(S01_lon,S01_lat,S01(:,:,decadenum),NHlon,NHlat,'linear');
    S01_int(NHdata==16 | NHdata==17) = 1; % saturation = 100% at water
    S01_int(UrbanArea) = 0; % no soil water on impermeable surface
    Ssd01_int = interpScattered(S01_lon,S01_lat,Ssd01(:,:,decadenum),NHlon,NHlat,'linear');
    Ssd01_int(NHdata==16 | NHdata==17) = 0; % saturation = 100% at water
    Ssd01_int(UrbanArea) = 1; % no soil water on impermeable surface
    
    % Crop factors (Makkink)
    Kcb = MakkinkMap(NHdata,daynum,factsvat_split,LGN2SIMGROcat);
    Kcb(Kcb == 65535)=NaN;    
    
    % Make LAI grid
    LAIgrid=constructLAI(Output, NHdata, decadenum);
    % % % % % need to check LAI values for edge effects for roads&stuff
    
    % Plant Height
    PlantHeight = Rh/100; % needs real values
    
    disp('updated daily boundary conditions')
    toc
    
    %% Calculate daily variables
    disp('calculating daily variables...')
    % Calculate conductivity (Dry, solid, saturated) and total Conductivity
    Ksoil = calc_soilconduc(Saturation_int, Ksolid, FracPorosity, K_Water); %[W m-1 K-1]
    Ksoil(NHdata == 16 | NHdata == 17) = Ksoil(NHdata == 16 | NHdata == 17) .* EddyConducMult; %mixing of water
    
    % Calculate soil water stress coefficient reducing Kcb
    Ks = calc_ks(S01_int, Ssd01_int);
    
    % Calculate soil evaporation reducing coefficient
    Kr= calc_kr(TEW, REW, Saturation_int, rootzone_mm_int);
    
    % Calculate volumetric heat capacity
    VHC_Soil =  VHC_Soil_Dry + VHC_Water * FracPorosity .* Saturation_int;  % [J m-3 K-1]
        % readjust heatcontent based on temperature (change in heat
        % capacity has to be reflected in heat content)
        HeatCont(:,:,2) = Temp(:,:,2) .* (VHC_Soil .* ThickSurfL); % [J m-2]
        for k= 3:NrLayer
            HeatCont(:,:,k) = Temp(:,:,k) .* (VHC_Soil * ThickL);
        end
    
    disp('calculated daily variables')
    
    %% Hourly part
    disp('starting hourly calculations...')
    for hour = 1:24
        %% Initialise hourly data
        % Interpolate ERA5
        ERA5t2m_int = interp2(ERA5longrid,ERA5latgrid,ERA5t2m_adj(:,:,Time)',NHlon,NHlat);
        ERA5d2m_int = interp2(ERA5longrid,ERA5latgrid,ERA5d2m_adj(:,:,Time)',NHlon,NHlat);
        ERA5net_rad_int = interp2(ERA5longrid,ERA5latgrid,ERA5net_rad(:,:,Time)',NHlon,NHlat);
        ERA5ssr_adj_int = interp2(ERA5longrid,ERA5latgrid,ERA5ssr_adj(:,:,Time)',NHlon,NHlat);
        ERA5str_adj_int = interp2(ERA5longrid,ERA5latgrid,ERA5str_adj(:,:,Time)',NHlon,NHlat);
        ERA5u_v_2m_int = interp2(ERA5longrid,ERA5latgrid,ERA5u_v_comb_2m(:,:,Time)',NHlon,NHlat);
        %ERA5sshf_int = interp2(ERA5longrid,ERA5latgrid,ERA5sshf_adj(:,:,Time)',NHlon,NHlat);
        
              
        %% Dynamic boundary conditions:
        % Calculate cover fraction.
        SolarAngle = calc_SolarAngle(daynum,hour-1);
        CoverFraction = calc_coverfrac(LAIgrid, SolarAngle); % [0-1] 1 is completely covered

        % Calculate Surface radiation.
        SurfaceRad = ERA5str_adj_int + ERA5ssr_adj_int - ERA5ssr_adj_int .* CoverFraction;    % [J m-2]
        if Debugging == 0
            clear('CoverFraction')
        end
        
        % Calculate temperatures
        % ERA5t2m as upper boundary condition for the vegetation layer.
        Temp(:,:,1) = ERA5t2m_int;                   % [K]
        
%         % ERA5 soil temp layer 3 as lower boundary condition
%         Temp(:,:,NrLayer) = interp2(ERA5longrid,ERA5latgrid,ERA5stl3_adj(:,:,StartTime)',NHlon,NHlat);  % [oC]
        
        % Calculate the surface temperature(T2)   J m-2 /    [J/m2/K ]       
        Temp(:,:,2) = Temp(:,:,2) + (SurfaceRad ./ (VHC_Soil .* ThickSurfL)); % [K]
        if Debugging == 0
            clear('SurfaceRad')
        end
        % Readjust heatcontent based on temperature (change in heat
        % has to be reflected in heat content)
        HeatCont(:,:,2) = Temp(:,:,2) .* (VHC_Soil .* ThickSurfL); % [J m-2]
        
        %% Rates: Flow and net flow calculations
        % Calculate flows.
        [Flow, GroundFlow] = calc_flows(Flow, AirDensity, AirHeatCapacity, ...
            Temp, Rh, Ksoil, ThickL, ThickSurfL, NrLayer, CellSizeX, CellSizeY);           % [W m-2]

        % Update HeatContent of soil and surface layers                       
        for k = 2 : NrLayer 
            HeatCont(:,:,k) = HeatCont(:,:,k) + Flow(:,:,k) * TimeStep - ...
                Flow(:,:,k-1) * TimeStep;       % [J m-2]
        end

        % Update the temperatures of the surface and soil layers                 % [oC]
        Temp(:,:,2) = (HeatCont(:,:,2)) ./ (VHC_Soil .* ThickSurfL);
        for k= 3:NrLayer
            Temp(:,:,k) = (HeatCont(:,:,k)) ./ (VHC_Soil * ThickL);
        end
        
        %% Evaporation
        % Calculate Kc_max, the upper limit of evaporation on a surface
        RHmin = calc_RHmin(ERA5t2m_int,ERA5d2m_int); % relative humidity [%]
        Kc_max = calc_kc_max(ERA5u_v_2m_int, RHmin, PlantHeight, Kcb);   
        if Debugging == 0
            clear('ERA5d2m_int','ERA5u_v_2m_int','RHmin')
        end
        
        % Calculate Potential ET (ET0)
        lET0 = calc_eto_makkink(ERA5ssr_adj_int,ERA5t2m_int,TimeStep);     % [m s-1]
        
        if Debugging == 0
            clear('ERA5t2m_int','ET0')
        end
        
        % Calculate soil evaporation coefficient
        Ke = Kr .* (Kc_max - Kcb);
        if Debugging == 0
            clear('Kc_max')
        end
        
        % Calculate actual Evapotranspiration
        lETact = calc_etact(lET0, Kcb, Ks, Ke);  % [W m-2]
        %lETact(UrbanArea)=0; % no evaporation above impermeable surfaces
        if Debugging == 0
            clear('Ke','lET0')
        end
        
        %% Finally determine SHF for current hour.
        % Convert units to W/m2 and calculate SHF
        lE = lETact * DensWater * LatentHeatEvap;  % [W m-2]
        R = ERA5net_rad_int ./ 3600;     % [W m-2]
        G = GroundFlow;                  % [W m-2] make it positive for downward flow
        if Debugging == 0
            clear('GroundFlow','ETact')
        end
        
        
        %% The moment we've all been waiting for!
        SHF = R - lE - G;               % [W m-2]
                
        
        %% visualise
%         figure(1)
% %         subplot(1,2,1)
%         imagesc([NHlon(end,1),NHlon(end,end)],[NHlat(end,1),NHlat(end,end)], ...
%             SHF,[-1000,2000])
%         colormap('jet')
%         colorbar
%         title( ['Sensible heat', string(ERA5datetime(Time)) ] )
%         drawnow
        
%         subplot(1,2,2)
%         imagesc([NHlon(end,1),NHlon(end,end)],[NHlat(end,1),NHlat(end,end)],ERA5sshf_int./3600.*-1,[-1000,2000])
%         colormap('jet')
%         colorbar
%         title( ['ERA5 sensible heat', string(ERA5datetime(Time)) ] )
%         drawnow
        
        
        
        %% store calculated sensible heat flux
%         SHF = int16(SHF);
%         mkdir([pwd,'\Output\MatFiles\day',num2str(daynum)])
%         outputloc=[pwd,'\Output\MatFiles\day',num2str(daynum),'\SHF_hour',num2str(Time),'.mat'];
%         save(outputloc,'SHF')
        % still need to get max SHF and aggregate to 25m resolution
        
        
        TempStore(:,StoreInd) = Temp(1,1,:);
        TempWaterStore(:,StoreInd) = Temp(1,6,:);
        EvapoStore(StoreInd,:) = lE;
        shfStore(StoreInd,:) = SHF;
        CoverStorer(StoreInd,:) = CoverFraction;

        StoreInd = StoreInd +1;
        
        save('TempStoreFlowerBulb.mat','shfStore')
        fprintf('calculated hour %i\n', hour);    
        Time = Time + 1;
        
    end % end hourly calculations
    
    fprintf('calculated day %i\n', daynum);
end % end daily calculations

%% Wrapping up
% Open the profiler and give a notification that the model has finished
%profile viewer
%beep
disp('script completed')
 disp('run SHF_output to create a video')
 %SHF_output

figure(1)
plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),shfStore(:,1),'m')
hold on
plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),shfStore(:,2),'g')
plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),shfStore(:,3),'y')
plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),shfStore(:,4),'r')
plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),shfStore(:,5),'c')
plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),shfStore(:,6),'b')

title('SHF various land-use classes')
legend('bulbs', 'grass', 'cereal crop', 'road', 'deciduous forest', 'lake')
xlabel('W m-2')
hold off

figure(2)
plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),TempStore'-273.15)
title('temp')

figure(4)
plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),TempWaterStore'-273.15)
title('temp water')

figure(3)
plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),EvapoStore(:,1),'m')
hold on
plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),EvapoStore(:,2),'g')
plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),EvapoStore(:,3),'y')
plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),EvapoStore(:,4),'r')
plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),EvapoStore(:,5),'c')
plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),EvapoStore(:,6),'b')

title('Evaporations various land-use classes')
legend('bulbs', 'grass', 'cereal crop', 'road', 'deciduous forest', 'lake')
xlabel('W m-2')
hold off

% figure(4)
% plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),CoverStorer(:,1),'m')
% hold on
% plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),CoverStorer(:,2),'g')
% plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),CoverStorer(:,3),'y')
% plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),CoverStorer(:,4),'r')
% plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),CoverStorer(:,5),'c')
% plot(ERA5datetime(StartTime:StartTime+size(TempStore,2)-1),CoverStorer(:,6),'b--')
% 
% title('Cover Fractions')
% legend('bulbs', 'grass', 'cereal crop', 'road', 'deciduous forest', 'lake')
% xlabel('deg')
% hold off
