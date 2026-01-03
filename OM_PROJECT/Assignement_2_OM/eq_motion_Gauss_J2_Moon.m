function plotKepFiltered(t, kep, mu, filterMode)
%PLOTKEPFILTERED Plot Keplerian element evolution (raw vs filtered) and secular trend.
%
%   plotKepFiltered(t, kep, mu, filterMode)
%
% PROTOTYPE
%   plotKepFiltered(t, kep, mu, filterMode)
%
% DESCRIPTION
%   This function visualises the time history of the six classical Keplerian
%   elements provided as a time series. The elements are plotted both in
%   their processed (continuous-angle) form and after low-pass filtering via
%   a moving average (movmean), in order to separate short-period
%   oscillations from long-period and/or secular behaviour.
%   A linear fit is also computed on the filtered signal to highlight the
%   mean secular trend over the selected time span. Time is normalised by
%   the orbital period computed from the initial semi-major axis.
%   The filtering window is selected according to the user-defined
%   filterMode to emphasise either secular evolution, abrupt events, or a
%   compromise between the two.
%
% INPUT
%   t          [Nx1] Time vector                                     [s]
%   kep        [Nx6] Keplerian element history:
%                   kep(:,1) = a     Semi-major axis                 [km]
%                   kep(:,2) = e     Eccentricity                    [-]
%                   kep(:,3) = i     Inclination                     [rad]
%                   kep(:,4) = Omega RAAN                            [rad]
%                   kep(:,5) = omega Argument of pericenter          [rad]
%                   kep(:,6) = f     True anomaly                    [rad]
%
%   mu         [1x1] Gravitational parameter of the central body     [km^3/s^2]
%
%   filterMode [char/string] Filtering mode selector:
%                   "secular" -> ~1-orbit window, suppresses periodic terms
%                   "events"  -> ~6-hour window, preserves fast changes
%                   "mixed"   -> ~1-day window, robust compromise
%
% OUTPUT
%   (none)     Generates one figure for each Keplerian element, showing:
%             (i) processed (raw) signal, (ii) filtered signal, and
%             (iii) linear secular trend estimated from the filtered signal.
%
% ASSUMPTIONS
%   - The input Keplerian angles are expressed in radians.
%   - Angles (i, Omega, omega) are unwrapped to avoid 2*pi discontinuities.
%   - The true anomaly is made continuous by enforcing positive 0..2*pi
%     increments (monotone evolution), which is suitable for visualisation
%     over multiple revolutions.
%   - The time step is approximately uniform (dt estimated by median(diff(t))).
%
% CONTRIBUTORS
%   Luca Deli
%
% VERSION
%   2026-01-03

    a0 = kep(1,1);
    T  = 2*pi*sqrt(a0^3/mu);
    tau = t / T;
    N = length(tau);

    dt = median(diff(t)); % seconds between samples

    % -------------------- FILTERING CHOICE (justified) --------------------
    % "secular": use 1 orbit window -> removes periodic oscillations, keeps drift (J2 / long-term evolution)
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
        kep_plot(:,j) = rad2deg(unwrap(kep(:,j)));   % i, Omega, omega
    end
    
    ang = kep(:,6);                 % f in rad
    d   = diff(ang);
    d   = mod(d + 2*pi, 2*pi);      % force positive increments (0..2pi)
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
