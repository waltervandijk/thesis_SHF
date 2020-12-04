% script to calculate sensible heat with SEBAL
% uses modis imagery

%% initialisation
clear
close all

% setting up environment
addpath(genpath('Data'), genpath('functions'))

% loading data
load ('MODIS060618_dayLST.mat')
    MODIS060618_data = MODIS060618_data(240:340,160:260);
    MODIS060618_lat = MODIS060618_lat(240:340,160:260);
    MODIS060618_lon = MODIS060618_lon(240:340,160:260);
load ('ERA5_comb.mat', 'ERA5t2m_comb', 'ERA5time_comb', 'ERA5lat_comb', ...
    'ERA5lon_comb', 'ERA5ssr_comb', 'ERA5str_comb')
    ERA5lat = repmat(ERA5lat_comb',46,1);
    ERA5lon = repmat(ERA5lon_comb,1,36);
    
load ('LAIstats.mat')
load ('lgn2018nh.mat')
    NHdata=double(NHdata); NHdata(NHdata==65535)=NaN;
wind_v100 = ncread('ERA5wind100m.nc','v100');
wind_u100 = ncread('ERA5wind100m.nc','u100');
    wind_100m = sqrt(wind_v100.^2+wind_u100.^2);
wind_lon_vec = ncread('ERA5wind100m.nc','longitude');
wind_lat_vec = ncread('ERA5wind100m.nc','latitude');
wind_lon = repmat(wind_lon_vec,1,15);
wind_lat = repmat(wind_lat_vec',19,1);
disp('data loaded')

% initialise time parameters
ERA5datenum = datenum(1900,1,1,double(ERA5time_comb),0,0);
ERA5datetime = datetime(datestr(ERA5datenum));%,'InputFormat','dd-MM-yyyy HH:mm:ss');
datainfo = hdfinfo('Data\MODIS\MOD11A1.A2018157.h18v03.006.2018158191853.hdf');
LST_time = hdfread(datainfo.Vgroup.Vgroup(1).SDS(3)); 
disp('Time = 10.6')
disp(ERA5datetime(2341))
disp(ERA5datetime(2342))

% adjust data
ERA5str_adj = ERA5cum2add (ERA5str_comb, ERA5datetime); %[J m-2]
ERA5ssr_adj = ERA5cum2add (ERA5ssr_comb, ERA5datetime); %[J m-2]
ERA5Rn = ERA5ssr_adj+ERA5str_adj;

% interpolate parameters in time
Rn_interp_A = 0.6*ERA5Rn(:,:,2342)+0.4*ERA5Rn(:,:,2341); % account for time
T2m_interp_A = 0.6*ERA5t2m_comb(:,:,2342)+0.4*ERA5t2m_comb(:,:,2341); % account for time
wind100m_interp_A = 0.6*wind_100m(:,:,2342)+0.4*wind_100m(:,:,2341); % account for time

% interpolate parameters in space
Rn_interp = interpn(ERA5lat',ERA5lon',Rn_interp_A',MODIS060618_lat,MODIS060618_lon);
T2m_interp = interpn(ERA5lat',ERA5lon',T2m_interp_A',MODIS060618_lat,MODIS060618_lon);
wind100m_interp = interpn(wind_lat',wind_lon',wind100m_interp_A',MODIS060618_lat,MODIS060618_lon);

disp('data initialised')

%% adjust data
LST_day = MODIS060618_data;
LST_day = 0.02*LST_day; % scale factor
LST_day(LST_day == 0) = NaN;
LST_day = LST_day - 273.15; % K to C
LST_day(LST_day < 0) = NaN;

%% control constants
decadenum = 16;

%% system constants
Zx = 100;           % m (height of wind measurements)
z2 = 2;             % m (height of air temp)
z1 = 0.1;           % m
Cp = 1004;          % J Kg-1 K-1
Karman = 0.41;      % Karman constant

%% calculations

disp('starting calculations')

% calculate Zom
LAI = constructLAI(LAIdata,double(NHdata),decadenum);
Zom_moment = make_roughness(LAI, double(NHdata));           %lgn2018 resolution

% conversion to modis resolution
modis_pix = ones(size(MODIS060618_data));
modis_pix(:) = find(modis_pix);                     % number the modis pixels
lgn_modis_pix = interp2(MODIS060618_lat,MODIS060618_lon,modis_pix,NHlat,NHlon,'nearest');   % assign lgn pixels to modis pixels with nearest neighbour

Zom_moment_modis = modis_pix;   % initialise new matix with index values
Zom_moment_modis(modis_pix<min(lgn_modis_pix(:)) | modis_pix>max(lgn_modis_pix(:)))=NaN; % remove indices outside lgn area
disp('converting lgn resolution to modis resolution')
for i = min(lgn_modis_pix(:)):max(lgn_modis_pix(:))
    Zom_moment_modis(i) = mean(Zom_moment(lgn_modis_pix==i));
    if mod(i,50)<1
        disp([num2str(((i-min(lgn_modis_pix(:)))/(max(lgn_modis_pix(:))-min(lgn_modis_pix(:))))*100),' %'])
    end
end

% calculate friction velocity
FricVel = (Karman.*wind100m_interp) ./ log(Zx./Zom_moment_modis);  %modis resolution

% calculate Rah
Rah = (log(z2/z1))./(FricVel.*Karman); %modis resolution [s m-1]

% calculate dT with hot/cold pixel
% cold: 31,71:  lake Alkamardermeer         (H = 0)
% hot:  50,94:  city centre of Amsterdam
coldH = 0;
hotH = 0.5*Rn_interp(50,94)/3600;   %0.2*Rn (Rigo&Parlow, 2006)

dT_cold = 0;
%dT_cold = T2m_interp(31,71)-LST_day(31,71);
AirDensHot = 101325./(287.05.*(T2m_interp(50,94)));
dT_hot = hotH*Rah(50,94)/(AirDensHot .* Cp);
%dT_hot = LST_day(50,94)-T2m_interp(50,94);

dTmdl = fitlm([LST_day(31,71), LST_day(50,94)],[dT_cold,dT_hot]);

dT = dTmdl.Coefficients.Estimate(1)+dTmdl.Coefficients.Estimate(2).*LST_day;
% first estimate of H
%dT = LST_day - T2m_interp; %modis resolution
AirDens = 101325./(287.05.*(273.15+LST_day+dT)); % kg m-3
SensHeat = (AirDens .* Cp .* dT) ./ Rah; % J s-1 m-2

disp('finished calculations')
disp('making figure')
figure(1)
    imagesc(SensHeat')
    colorbar
    colormap jet
    title('first estimate')

%% iterative process

% Monin-Obukhov Length (L)
L = -1 * (AirDens.*Cp.*(FricVel.^3).*(LST_day+273.15))./(Karman.*9.81.*SensHeat);

for iterate = 1:10
    MomentCorr100m = zeros(size(L));
    HeatTransCorr2m = zeros(size(L));
    HeatTransCorr01m = zeros(size(L));
    for i = 1:size(L,1)
        for j = 1:size(L,2)
            if L(i,j)<0
                x100m = (1-16*(100./L(i,j))).^0.25;
                MomentCorr100m(i,j) = 2*log((1+x100m)/2)+log((1+x100m.^2)/2)-2.*atan(x100m+0.5*pi);

                x2m = (1-16*(2./L(i,j))).^0.25;
                HeatTransCorr2m(i,j) = 2*log((1+x2m.^2)/2);

                x01m = (1-16*(0.1./L(i,j))).^0.25;
                HeatTransCorr01m(i,j) = 2*log((1+x01m.^2)/2);

            elseif L(i,j)>0
                MomentCorr100m(i,j) = -5*(2/L(i,j));
                HeatTransCorr2m(i,j) = -5*(2/L(i,j));
                HeatTransCorr01m(i,j) = -5*(0.1/L(i,j));

            elseif L(i,j)==0
                MomentCorr100m(i,j) = 0;
                HeatTransCorr2m(i,j) = 0;
                HeatTransCorr01m(i,j) = 0;
            end
        end
    end

    % correct friction velocity
    FricVelCorr = (Karman.*wind100m_interp) ./ (log(Zx./Zom_moment_modis)-MomentCorr100m);

    % correct Rah
    RahPrev = Rah;
    Rah = (log(z2./z1)-HeatTransCorr2m+HeatTransCorr01m)./(FricVelCorr.*Karman);
    % correct Sensible heat
    SensHeat = (AirDens .* Cp .* dT) ./ Rah;
    figure(4)
        imagesc(SensHeat',[-100 1500])
        colorbar
        title('SHF SEBAL MODIS')
        drawnow
        pause(1)
    % for re-iteration
    L = -1 * (AirDens.*Cp.*(FricVelCorr.^3).*(LST_day+273.15))./(Karman.*9.81.*SensHeat);
end

% remove outliers
SensHeat(SensHeat>2500)=NaN;
    figure(4)
        imagesc(SensHeat',[-100 1500])
        colorbar
        title('SHF SEBAL MODIS')
        drawnow
        pause(2)
        
%% comparison
load('D:\Documenten\AcademicStuff\EarthScienceFPES\ThesisThermal_LisaOut\OutputMerged\MatFiles\day157\SHF_hour2341.mat')
SHF_comb_alt=0.4*SHF_comb;
load('D:\Documenten\AcademicStuff\EarthScienceFPES\ThesisThermal_LisaOut\OutputMerged\MatFiles\day157\SHF_hour2342.mat')
SHF_comb_alt=SHF_comb_alt+0.6*SHF_comb;
SHF_comb_alt(SHF_comb_alt>2000)=0;

    figure(3)
        imagesc(SHF_comb_alt,[-100 1500])
        colorbar
        drawnow
        pause(2)

% aggregate mdoel output to modis resolution
SHF_comb_modis = modis_pix;   % initialise new matix with index values
SHF_comb_modis (modis_pix<min(lgn_modis_pix(:)) | modis_pix>max(lgn_modis_pix(:)))=NaN; % remove indices outside lgn area
disp('converting SHF model to modis resolution')
for i = min(lgn_modis_pix(:)):max(lgn_modis_pix(:))
    SHF_comb_modis (i) = mean(SHF_comb_alt(lgn_modis_pix==i));
    if mod(i,50)<1
        disp([num2str(((i-min(lgn_modis_pix(:)))/(max(lgn_modis_pix(:))-min(lgn_modis_pix(:))))*100),' %'])
    end
end

    figure(3)
        imagesc(SHF_comb_modis',[-100 1500])
        colorbar
        title('SHF model')
        drawnow
        pause(2)

disp('completion')

%% saving results

SHF_from_model = SHF_comb_modis';
SHF_from_SEBAL = SensHeat';
model_cells=SHF_from_model(~(isnan(SHF_from_SEBAL)));
SEBAL_cells=SHF_from_SEBAL(~(isnan(SHF_from_SEBAL)));
save('compare_model_SEBAL_MODIS_1', 'SHF_from_model', 'SHF_from_SEBAL', 'model_cells','SEBAL_cells')


%% visualisation
figure(11)
    subplot(1,2,1)
        imagesc(SHF_from_model, [-100 1000])
        title('SHF from model')
        colorbar
    subplot(1,2,2)
        imagesc(SHF_from_SEBAL, [-100 1000])
        title('SHF from SEBAL')
        colorbar
        
    figure(12)
        boxplot([model_cells,SEBAL_cells])
        xticklabels([{'Model'} , {'SEBAL'}])
        ylabel('W m^{-2}')
        title('Box plot of the model and SEBAL n=4057')

    figure(13)
        plot(model_cells,SEBAL_cells,'.k', 'MarkerSize', 1)
        hold on
        plot([-400:800],[-400:800],'b')
        xlabel('SHF from model [W m^{-2}]')
        ylabel('SHF from SEBAL [W m^{-2}]')
        title('Correlation between the model and SEBAL')
        
%% figure for in report
figure(100)
    set(gcf,'position',[10,10,1080,1080])
    subplot(2,2,1)
        imagesc(SHF_from_model, [-100 1000])
        title('SHF from model')
        c=colorbar;
        title(c, 'W m^{-2}')
    subplot(2,2,2)
        imagesc(SHF_from_SEBAL, [-100 1000])
        title('SHF from SEBAL')
        c=colorbar;
        title(c, 'W m^{-2}')
    subplot(2,2,3)
        boxplot([model_cells,SEBAL_cells])
        xticklabels([{'Model'} , {'SEBAL'}])
        ylabel('W m^{-2}')
        title('Box plot of the model and SEBAL n=4057')       
    subplot(2,2,4)
        plot(model_cells,SEBAL_cells,'.k', 'MarkerSize', 1)
        hold on
        plot([-400:800],[-400:800],'b')
        xlabel('SHF from model [W m^{-2}]')
        ylabel('SHF from SEBAL [W m^{-2}]')
        title('Correlation between the model and SEBAL')        
saveas(gcf,'figure_run1.fig')
saveas(gcf,'figure_run1.png')