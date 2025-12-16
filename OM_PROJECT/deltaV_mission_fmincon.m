function [MISSION, x_opt, fval, exitflag, output] = deltaV_mission_fmincon(MISSION, x0)
% DELTAVMISSIONFMINCON Continuous optimization of the minimum total DeltaV between two time windows
% using fmincon, with constraints on t2 > t1 and ||deltaV1|| <= V_launcher.
%
% INPUT:  MISSION : Structure that must contain:
%           .TimeWindow.Departure_v (2x6 datetime or [y m d] vectors)
%           .TimeWindow.Arrival_v   (2x6 datetime or [y m d] vectors)
%           .Data.dep_planet, .Data.target_planet
%           .Data.mu_sun
%           .Data.v_launch
%         x0  : initial guess [t1_0, t2_0] in MJD2000.
%             If not provided, uses the midpoint of the windows.
% OUTPUT: MISSION : updated structure with MostEfficient fields
%         x_opt   : optimal solution [t1*, t2*] in MJD2000
%         fval    : DeltaV_tot(t1*, t2*)
%         exitflag, output : diagnostics from fmincon
%
% For other details, see DELTAV_MISSION function.

% Extract time windows and convert to MJD2000 ----
t1 = [date2mjd2000(MISSION.TimeWindow.Departure_v(1,:)), ...
      date2mjd2000(MISSION.TimeWindow.Departure_v(2,:))];
t2 = [date2mjd2000(MISSION.TimeWindow.Arrival_v(1,:)), ...
      date2mjd2000(MISSION.TimeWindow.Arrival_v(2,:))];

lb = [t1(1), t2(1)];            % lower bounds [t1_min, t2_min]
ub = [t1(2), t2(2)];            % upper bounds [t1_max, t2_max]

% Initial guess (midpoint of windows) if not provided
if nargin < 2 || isempty(x0)
    x0 = [mean(t1), mean(t2)];
end

%Setup fmincon
fun     = @(x) objectiveDeltaV(x, MISSION);
nonlcon = @(x) nonlinearConstraints(x, MISSION);

options = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'Display','iter', ...
    'MaxIterations', 500, ...
    'MaxFunctionEvaluations', 5000, ...
    'OptimalityTolerance', 1e-12, ...
    'ConstraintTolerance', 1e-13, ...
    'StepTolerance', 1e-14);

%Call fmincon 
[x_opt, fval, exitflag, output] = fmincon(fun, x0, [], [], [], [], lb, ub, nonlcon, options);

% Populate MISSION.MostEfficient with results 
[DV_tot, DV1_norm, DV2_norm, TOF_s, v1_t, v2_t, v1_p, v2_p] = evaluateTransfer(x_opt(1), x_opt(2), MISSION);

MISSION.MostEfficient.deltaV_TOT = DV_tot;
MISSION.MostEfficient.deltaV1    = DV1_norm;
MISSION.MostEfficient.deltaV2    = DV2_norm;
MISSION.MostEfficient.TOF        = TOF_s;                    % [s]
MISSION.MostEfficient.t1         = x_opt(1);                 % MJD2000
MISSION.MostEfficient.t2         = x_opt(2);                 % MJD2000
MISSION.MostEfficient.Departure  = datetime( mjd20002date(x_opt(1)) );
MISSION.MostEfficient.Arrival    = datetime( mjd20002date(x_opt(2)) );

[MISSION.MostEfficient.Eph_dep_planet,~] = uplanet(MISSION.MostEfficient.t1, MISSION.Data.dep_planet);
[MISSION.MostEfficient.Eph_target_planet,~] = uplanet(MISSION.MostEfficient.t2, MISSION.Data.target_planet);

MISSION.MostEfficient.v_t1      = v1_p;
MISSION.MostEfficient.v_t2      = v2_p;
MISSION.MostEfficient.v_t1_tr   = v1_t;
MISSION.MostEfficient.v_t2_tr   = v2_t;

% Objective function: total DeltaV 
function DV = objectiveDeltaV(x, M)
    t1_loc = x(1);
    t2_loc = x(2);
    if t2_loc <= t1_loc
        DV = 1e9; % heavy penalty if arrival before departure
        return;
    end
    [DV_tot, ~, ~, ~, ~, ~, ~, ~] = evaluateTransfer(t1_loc, t2_loc, M);
    if ~isfinite(DV_tot)
        DV = 1e9;
    else
        DV = DV_tot;
    end
end

% Nonlinear constraints
function [c, ceq] = nonlinearConstraints(x, M)
    t1_loc = x(1);
    t2_loc = x(2);

    % Constraint 1: t2 >= t1
    eps_days = 1e-9;
    c1 = (t1_loc - t2_loc) + eps_days;

    % If v_launch = inf, ignore deltaV1 constraint
    if isinf(M.Data.v_launch)
        c2 = -1e6; % Always satisfied
    else
        [~, DV1_norm, ~, ~, ~, ~, ~, ~] = evaluateTransfer(t1_loc, t2_loc, M);
        if ~isfinite(DV1_norm)
            c2 = 1e6; % infeasible if Lambert fails
        else
            c2 = DV1_norm - M.Data.v_launch;
        end
    end

    c   = [c1; c2];
    ceq = [];
end

% Transfer evaluation 
function [DV_tot, DV1_norm, DV2_norm, TOF_s, v_t1_tr, v_t2_tr, v_t1_pl, v_t2_pl] = evaluateTransfer(t1_eval, t2_eval, M)
    [eph1,~] = uplanet(t1_eval, M.Data.dep_planet);
    eph1 = [eph1(1) eph1(2) rad2deg(eph1(3)) rad2deg(eph1(4)) rad2deg(eph1(5)) rad2deg(eph1(6))];
    [s1,~] = kep2car(eph1(1),eph1(2),eph1(3),eph1(4),eph1(5),eph1(6),M.Data.mu_sun);
    r_t1   = s1(1:3);
    v_t1_pl = s1(4:6);

    [eph2,~] = uplanet(t2_eval, M.Data.target_planet);
    eph2 = [eph2(1) eph2(2) rad2deg(eph2(3)) rad2deg(eph2(4)) rad2deg(eph2(5)) rad2deg(eph2(6))];
    [s2,~] = kep2car(eph2(1),eph2(2),eph2(3),eph2(4),eph2(5),eph2(6),M.Data.mu_sun);
    r_t2   = s2(1:3);
    v_t2_pl = s2(4:6);

    TOF_s = seconds(days(t2_eval - t1_eval));
    if TOF_s <= 0
        DV_tot    = NaN; DV1_norm = NaN; DV2_norm = NaN;
        v_t1_tr   = [NaN NaN NaN]; v_t2_tr = [NaN NaN NaN];
        return;
    end

    try
        [~,~,~,~,v_t1_tr, v_t2_tr, ~, ~] = lambertMR(r_t1, r_t2, TOF_s, M.Data.mu_sun, 0);
        v_t1_tr = v_t1_tr'; v_t2_tr = v_t2_tr';
    catch
        DV_tot    = NaN; DV1_norm = NaN; DV2_norm = NaN;
        v_t1_tr   = [NaN NaN NaN]; v_t2_tr = [NaN NaN NaN];
        return;
    end

    deltaV1   = v_t1_tr - v_t1_pl;
    deltaV2   = v_t2_pl - v_t2_tr;
    DV1_norm  = norm(deltaV1);
    DV2_norm  = norm(deltaV2);
    DV_tot    = DV1_norm + DV2_norm;
end

end
