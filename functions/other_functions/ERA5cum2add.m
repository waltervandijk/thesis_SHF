function ERA5Out = ERA5cum2add(ERA5In,ERA5datetime)
% ERA5 data adjustment
% Walter van Dijk
% University of Amsterdam

%% Notes
% All data are daily cumulative values.
% The first one '01-May-2018 00:00:00' is the total for '30-April-2018'
% To get hourly data the previous value has to be subtracted from the
% current value. This has to be done for all values exept the one of which
% the hour is '01:00:00'.
% Also, data is transformed

%% Inputs
% ERA5In = Input variable with cumulative values

%% Outputs
% ERA5Out = Output variable with hourly data
ERA5Out(:,:,:)=zeros(size(ERA5In,1,2,3));
steps = size(ERA5In,3);

for i=3:steps
    datechar = char(ERA5datetime(i));
    if sum(datechar(end-7:end)~=char("01:00:00"))==0
        ERA5Out(:,:,i)=ERA5In(:,:,i);
    else
        ERA5Out(:,:,i)=ERA5In(:,:,i)-ERA5In(:,:,i-1);
    end
end
    
    


end

