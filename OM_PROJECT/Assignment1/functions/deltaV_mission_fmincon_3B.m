function [MISSION, x_opt, fval, exitflag, output] = deltaV_mission_fmincon_3B(MISSION, x0)
% DELTAV_MISSION_3B_FMINCON Continuous optimization of a 3-body mission with powered flyby
% in the degrees of freedom x = [t_dep, TOF1, TOF2] in MJD2000 days
% The objective minimized is the total Delta-V: DV_tot = DV1 + DV_flyby + DV2
% - DV1  computed from a Lambert transfer between dep planet and flyby planet
% - DV2  computed from a Lambert transfer between flyby planet and asteroid
% - DV_flyby computed via powered flyby model (flyby_powered)
%
% INPUTS : MISSION : struct
%           Mission definition and parameters. Required fields:
%
%           MISSION.LEG1.TimeWindow.Departure_v : [2x6] datetime or date-vector
%                                                  Departure window bounds (start/end). 
%                                                   Converted to MJD2000 [days].
%
%           MISSION.LEG1.TimeWindow.TOF1_min    : minimum TOF1  [days]
%           MISSION.LEG1.TimeWindow.TOF1_max    : maximum TOF1  [days]
%
%           MISSION.LEG2.TimeWindow.TOF2_min    : minimum TOF2  [days]
%           MISSION.LEG2.TimeWindow.TOF2_max    : maximum TOF2  [days]
%
%           MISSION.LEG1.Data.dep_planet        : ID of departure planet (used by uplanet).
%           MISSION.FLYBY.Data.flyby_planet     : ID of flyby planet (used by uplanet).
%
%           MISSION.LEG2.Data.target_asteroid   : ID Target asteroid identifier (used by ephAsteroids).
%
%           MISSION.LEG1.Data.mu_sun            :   [km^3/s^2]
%           MISSION.LEG2.Data.mu_sun            :    [km^3/s^2]
%                                                   Gravitational parameter of the Sun
%
%           MISSION.LEG1.Data.v_launch          : Optional launcher constraint on DV1 [km/s]:
%                                                   DV1 <= v_launch
%                                                 If Inf, constraint is disabled.
%
%           MISSION.FLYBY.Data.mu_planet        :  [km^3/s^2]
%           MISSION.FLYBY.Data.R_planet         :  [km]
%           MISSION.FLYBY.Data.h_min            :  [km]
%                                                   Flyby planet constants and minimum 
%                                                   altitude constraint for
%                                                   powered flyby computation.
%
%            x0 : [1x3] double  [MJD2000 days, days, days]
%                           Initial guess for the optimizer:
%                           x0 = [t_dep0, TOF10, TOF20]
%                 If omitted or empty, a midpoint-based guess is used.
%
% OUTPUTS :  MISSION : struct
%           Input struct updated with solution fields under:
%           MISSION.MostEfficient.t_dep        : [1x1] double  [MJD2000 days]
%           MISSION.MostEfficient.TOF1         : [1x1] double  [days]
%           MISSION.MostEfficient.TOF2         : [1x1] double  [days]
%           MISSION.MostEfficient.DELTA_V1     : [1x1] double  [km/s]
%           MISSION.MostEfficient.DELTA_VFB    : [1x1] double  [km/s]
%           MISSION.MostEfficient.DELTA_V2     : [1x1] double  [km/s]
%           MISSION.MostEfficient.DELTA_V_TOT  : [1x1] double  [km/s]
%
%           x_opt       :   [1x3] double  [MJD2000 days, days, days]
%                           Optimal decision vector [t_dep*, TOF1*, TOF2*].
%   
%           fval        :   [1x1] double  [km/s]
%                           Optimal objective value (DV_tot at x_opt).
%
%            exitflag   :   [1x1] double/int
%                           fmincon termination flag.
%
%           output      :   struct
%                           fmincon diagnostic output structure.
%
% N.B.:
%   - States are computed via uplanet/ephAsteroids in Keplerian elements, then
%     converted to Cartesian (kep2car) with mu_sun.
%   - Lambert transfers are computed with lambertMR; TOFs are converted as:
%         dt = seconds(days(TOF)).
%   - Nonlinear constraints enforce:
%       TOF1 > 0, TOF2 > 0, and (optionally) DV1 <= v_launch.
%
% AUTHORS :     Amura Fabio
%

% Time windows and bounds 

t_dep_bounds = [ date2mjd2000(MISSION.LEG1.TimeWindow.Departure_v(1,:)), ...
                 date2mjd2000(MISSION.LEG1.TimeWindow.Departure_v(2,:)) ];

lb = [t_dep_bounds(1), ...
      MISSION.LEG1.TimeWindow.TOF1_min, ...
      MISSION.LEG2.TimeWindow.TOF2_min ];

ub = [t_dep_bounds(2), ...
      MISSION.LEG1.TimeWindow.TOF1_max, ...
      MISSION.LEG2.TimeWindow.TOF2_max ];

% Initial guess 
if nargin < 2 || isempty(x0)
    x0 = [ mean(t_dep_bounds), ...
           mean([lb(2), ub(2)]), ...
           mean([lb(3), ub(3)]) ];
end

% fmincon setup
fun = @(x) objective3B(x, MISSION);
nonlcon = @(x) nonlinearConstraints3B(x, MISSION);

options = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'Display','iter', ...
    'MaxIterations', 300, ...
    'MaxFunctionEvaluations', 4000, ...
    'OptimalityTolerance', 1e-12, ...
    'StepTolerance', 1e-14);

[x_opt, fval, exitflag, output] = ...
    fmincon(fun, x0, [], [], [], [], lb, ub, nonlcon, options);

% Store solution 
RESULT = evaluateMission3B(x_opt, MISSION);

MISSION.MostEfficient.t_dep  = x_opt(1);
MISSION.MostEfficient.TOF1   = x_opt(2);
MISSION.MostEfficient.TOF2   = x_opt(3);
MISSION.MostEfficient.DELTA_V1  = RESULT.DV1;
MISSION.MostEfficient.DELTA_VFB = RESULT.DVFB;
MISSION.MostEfficient.DELTA_V2  = RESULT.DV2;
MISSION.MostEfficient.DELTA_V_TOT = RESULT.DVtot;

end


function J = objective3B(x, M)

RES = evaluateMission3B(x, M);

if ~RES.feasible
    J = 1e6 + 1e3*(~isfinite(RES.DVtot));
    return
end

J = RES.DVtot;

end


function [c, ceq] = nonlinearConstraints3B(x, M)

% Retrieving Time of Flights
TOF1 = x(2);
TOF2 = x(3);

% Enforce positivity (robustness)
c1 = -TOF1;
c2 = -TOF2;

% Optional launcher constraint
RES = evaluateMission3B(x, M);

if isinf(M.LEG1.Data.v_launch)
    c3 = -1;
elseif ~RES.feasible || ~isfinite(RES.DV1)
    c3 = 1e9; 
else
    c3 = RES.DV1 - M.LEG1.Data.v_launch;
end

c = [c1; c2; c3];
ceq = [];

end


function RES = evaluateMission3B(x, M)

t_dep = x(1);
TOF1  = x(2);
TOF2  = x(3);

RES.feasible = false;
RES.DV1  = NaN;
RES.DVFB = NaN;
RES.DV2  = NaN;
RES.DVtot = NaN;

t_FB  = t_dep + TOF1;
t_arr = t_FB  + TOF2;

try
    % Planet states
    [kep_M,~] = uplanet(t_dep, M.LEG1.Data.dep_planet);
    [kep_E,~] = uplanet(t_FB,  M.FLYBY.Data.flyby_planet);
    [kep_A,~] = ephAsteroids(t_arr, M.LEG2.Data.target_asteroid);
    
    % Conversion to degrees
    kep_M(3:6) = rad2deg(kep_M(3:6));
    kep_E(3:6) = rad2deg(kep_E(3:6));
    kep_A(3:6) = rad2deg(kep_A(3:6));
    
    % Computing the state
    [sM,~] = kep2car(kep_M(1),kep_M(2),kep_M(3),kep_M(4),kep_M(5),kep_M(6),M.LEG1.Data.mu_sun);
    [sE,~] = kep2car(kep_E(1),kep_E(2),kep_E(3),kep_E(4),kep_E(5),kep_E(6),M.LEG1.Data.mu_sun);
    [sA,~] = kep2car(kep_A(1),kep_A(2),kep_A(3),kep_A(4),kep_A(5),kep_A(6),M.LEG2.Data.mu_sun);

    % Lambert legs
    dt1 = seconds(days(TOF1));
    dt2 = seconds(days(TOF2));

    [~,~,~,~,vM_t,vE_t1] = lambertMR(sM(1:3), sE(1:3), dt1, M.LEG1.Data.mu_sun, 0);
    [~,~,~,~,vE_t2,vA_t] = lambertMR(sE(1:3), sA(1:3), dt2, M.LEG2.Data.mu_sun, 0);

    vM_t = vM_t';
    vE_t1 = vE_t1';
    vE_t2 = vE_t2';
    vA_t = vA_t';

    DV1 = norm(vM_t - sM(4:6));
    DV2 = norm(vA_t - sA(4:6));
    
    % Flyby
    vinf_min  = vE_t1 - sE(4:6);
    vinf_plus = vE_t2 - sE(4:6);

    FLYBY = flyby_powered(vinf_min, vinf_plus, ...
        M.FLYBY.Data.mu_planet, ...
        M.FLYBY.Data.R_planet, ...
        M.FLYBY.Data.h_min);

    if isnan(FLYBY.rp)
        return
    end

    RES.DV1  = DV1;
    RES.DVFB = FLYBY.Deltavp;
    RES.DV2  = DV2;
    RES.DVtot = DV1 + RES.DVFB + DV2;
    RES.feasible = true;

catch
    return
end

end
