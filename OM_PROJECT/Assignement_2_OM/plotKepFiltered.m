function plotKepFiltered(t, kep, mu, filterMode)
%PLOTKEPRAWFILTERED  Plot keplerian element evolution (raw vs filtered).
%
% INPUT
%   t          time vector [s]
%   kep        Nx6 matrix [a e i Omega omega f] (a km, angles rad)
%   mu         gravitational parameter [km^3/s^2]
%   filterMode string specifying what to emphasize:
%              "secular"  -> highlights slow drift (J2 / long-term evolution)
%              "events"   -> highlights jumps / rapid changes (TLE discontinuities, maneuvers)
%              "mixed"    -> compromise (good default for long datasets)
%
% OUTPUT
%   (none)     Generates figures for each keplerian element

    a0 = kep(1,1);
    T  = 2*pi*sqrt(a0^3/mu);
    tau = t / T;
    N = length(tau);

    dt = median(diff(t)); % seconds between samples

    % -------------------- FILTERING CHOICE (justified) --------------------
    % "secular": use 1 orbit window -> removes periodic oscillations, keeps drift (J2, drag trend)
    % "events" : use fixed short window (6h) -> reduces noise but preserves sharp jumps in TLE
    % "mixed"  : use 1 day window -> smooth enough for trends, but still responsive to changes
    %
    % Why not same window always?
    % - 1 orbit window is best for physical secular behavior but can hide sudden TLE jumps
    % - 6h window is best for detecting abrupt changes but still noisy in angular elements
    % - 1 day is a robust compromise for multi-year TLE datasets

    switch lower(filterMode)
        case "secular"
            w = max(5, round(T/dt));         % ~1 orbit
            filterLabel = sprintf("Filtered (%.1f orbit)", w*dt/T);

        case "events"
            w = max(5, round((6*3600)/dt));  % 6 hours
            filterLabel = sprintf("Filtered (%.1f h)", w*dt/3600);

        case "mixed"
            w = max(5, round((24*3600)/dt)); % 1 day
            filterLabel = sprintf("Filtered (%.1f d)", w*dt/86400);

        otherwise
            error('filterMode must be "secular", "events" or "mixed"');
    end

    % -------------------- PREPARE ARRAYS FOR PLOT --------------------
    kep_plot = zeros(N,6);

    kep_plot(:,1) = kep(:,1);
    kep_plot(:,2) = kep(:,2);

    for j = 3:5
        kep_plot(:,j) = rad2deg(unwrap(kep(:,j)));   % i, Omega, omega normali
    end
    
    ang = kep(:,6);                 % f in rad
    d   = diff(ang);
    d   = mod(d + 2*pi, 2*pi);      % forza incrementi sempre positivi (0..2pi)
    ang_cont = [ang(1); ang(1)+cumsum(d)];
    kep_plot(:,6) = rad2deg(ang_cont);


    kep_filt = zeros(N,6);
    for j=1:6
        kep_filt(:,j) = movmean(kep_plot(:,j), w, 'Endpoints','shrink');
    end

    names = {'a','e','i','Raan','omega','f'};
    ylab  = {'a [km]','e [-]','i [deg]','Raan [deg]','omega [deg]','f [deg]'};

    % -------------------- PLOTS --------------------
for j = 1:6

    figure('Name', ['Keplerian evolution - ' names{j}]);
    hold on;

    plot(tau, kep_plot(:,j), 'LineWidth', 0.9);
    plot(tau, kep_filt(:,j), 'LineWidth', 1.4);

    % ---- secular trend line (linear fit) ----
    p = polyfit(tau, kep_filt(:,j), 1);
    trend = polyval(p, tau);
    plot(tau, trend, 'LineWidth', 2.2);

    grid on;
    axis tight;
    xlabel('time [T]');
    ylabel(ylab{j});
    title(['Keplerian element: ' names{j} '  ' filterMode '']);
    legend('Not Filtered', filterLabel, 'Secular linear trend', 'Location','south');

end

end
