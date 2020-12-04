function [Solar_Angle] = calc_SolarAngle(day,hour)
    %SOLARZENITH Calculates the solar zenith for given hour and day
    % applicable to the latitude of Noord Holland, the Netherlands 
    %
    % Walter van Dijk 2020

    %% INPUTS
    % day = day of year
    % hour = hour of day
    
    %% OUTPUT
    % Solar_Angle = angle of the sun - positive values for daytime, max
        % possible is 90 deg, however this does not occur in the Netherlands
    
    %% calculations
    
    latitude = 52.5; 
    longitude = 4.8;

    % difference with greenwich time
    LSTM = 0;

    % equation of time
    B = (360/365)*(day-81);
    EoT = 9.87*sind(2*B) - 7.53*cosd(B) - 1.5*sind(B);

    % the Earth rotates 1° every 4 minutes
    % time correction
    TC = 4*(longitude-LSTM)+EoT; % [min]

    % local solar time
    LST = hour + TC/60; % [hour]

    % hour angle, noon = 0deg
    HRA = 15 * (LST-12); % [deg]

    % declination
    declination = 23.45 * sind(B); % [deg]

    % Elevation and Azimuth
    a = asind(sind(declination).*sind(latitude)+cosd(declination) .* ...
        cosd(latitude).*cosd(HRA));
    %Azimuth = acosd((sind(declination).*cosd(latitude)+cosd(declination) .* ...
     %   sind(latitude).*cosd(HRA)) ./ cosd(a));

    Solar_Angle = a;

end

