function [d] = calc_d(T2m)
    % Calculates slope vapour pressure.
    % Used in van der Veldt, R. (2017). Challenges of modelling soaring flight in humid landscapes
    % Edited by Rens van der Veldt, original by Bart Sweerts

    %% [Input]
    
    % T2m = temperature at 2m.              [oC]
    
    %% [Output]
    
    % d = slope vapour pressure             [kPaC-1]
    
    %% Calculate slope vapour pressure
    d = (4098 .* (0.6108 * exp(17.27 .* T2m ./ (T2m + 237.3)))) ./ ((T2m + 237.3) .^ 2);
end