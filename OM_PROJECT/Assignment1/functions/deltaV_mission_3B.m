function [MISSION] = deltaV_mission_3B(MISSION,t_dep_step,tof_step)
%DELTAV_MISSION_3B This function exploits a 3D grid search for the mission specified. 
% The Injection Orbit from departure planet and
% Insertion Orbit to arrival planet ARE NOT CONSIDERED.
%
% INPUTS : MISSION : This struct contains all the informations for the
%                   mission, and should be fullfilled with the following
%                   data: 
%                   DATA : All info about planetary constants, ID of
%                   planets, V_launcher constraint, rp of flyby
%                   constraint, budget...
%                   TIMEWINDOW : Info about departure time window and
%                   TOF1_min, TOF1_max,TOF2_min,TOF2_max, specified in datetime format and
%                   vector: (yyyy,mm,d,h,m,s), and in MJD2000 date for TOFs
%          t_dep_step : scalar for departure time step  
%                           1 : each time step = 1 day (not recommended)
%                           2 : each time step = 7 days (not recommended)
%                           3 : each time step = 14 days
%                           4 : each time step = 30 days 
%                           5 : each time step = 60 days
%          tof_step : scalar for tof time step
%                           1 : each time step = 12 hours (not recommended)
%                           2 : each time step = 1 day (not recommended)
%                           3 : each time step = 4 days
%                           4 : each time step = 10 days 
%                           5 : each time step = 20 days
% 
% 
% OUTPUTS : MISSION : struct returned with the solutions of the 3-body mission
%             obtained via 3D grid search (t_dep, TOF1, TOF2).
%
%               The struct is organized by mission legs as follows:
%
%               ---------------- LEG 1 : Mars → Earth ----------------
%               MISSION.LEG1.Solution.DELTA_V1 :
%               3D matrix [N_dep x N_TOF1 x N_TOF2] containing the ΔV required at
%               Mars departure for each combination of (t_dep, TOF1, TOF2)
%               [km/s]
%
%               MISSION.LEG1.Solution.t_dep_v :
%               Vector of departure times from Mars in MJD2000 days
%
%               MISSION.LEG1.Solution.TOF1_v :
%               Vector of times of flight for LEG1 (Mars → Earth) in MJD2000 days
%
%
%               ---------------- FLYBY : Earth Powered Fly-by ----------------
%               MISSION.FLYBY.Solution.DELTAVFB :
%               3D matrix [N_dep x N_TOF1 x N_TOF2] containing the ΔV required at
%               Earth pericentre for the powered gravity assist
%               [km/s]
%
%
%               ---------------- LEG 2 : Earth → Asteroid ----------------
%               MISSION.LEG2.Solution.DELTAV2 :
%               3D matrix [N_dep x N_TOF1 x N_TOF2] containing the ΔV required at
%               arrival to the target asteroid
%               [km/s]
%
%               MISSION.LEG2.Solution.TOF2_v :
%               Vector of times of flight for LEG2 (Earth → Asteroid) in MJD2000 days
%
%
%               ---------------- GLOBAL MISSION COST ----------------
%               MISSION.DELTAV_TOT :
%               3D matrix [N_dep x N_TOF1 x N_TOF2] containing the total mission
%               ΔV, defined as:
%               ΔV_TOT = ΔV1 + ΔV_FB + ΔV2
%               [km/s]
%
%               MISSION.MostEfficient :
%               Struct containing the best (minimum ΔV) solution found in the grid,
%               useful for sanity checks and as an initial guess for refinement
%               methods (e.g., fmincon). Not intended as final mission design.
% 
% 
% N.B.: the dimensions of the matrices depends on the time_step_case
% specified
% 
% AUTHORS :     Amura Fabio
%




if(nargin<2)
    disp('Error: unspecified time step case');
end

% Departure time in MJD2000 format for LEG 1
t_dep = [date2mjd2000(MISSION.LEG1.TimeWindow.Departure_v(1,:)), date2mjd2000(MISSION.LEG1.TimeWindow.Departure_v(2,:))];

% Discretization of times. MJD2000 format
switch t_dep_step
    case 1
        t_dep_step = days(days(1)); % 1 day
    case 2
        t_dep_step = days(days(7)); % 7 days
    case 3 
        t_dep_step = days(days(14)); % 14 days
    case 4
        t_dep_step = days(days(30)); % 30 days
    case 5 
        t_dep_step = days(days(60)); % 60 days
    otherwise 
        disp('Error: time step case not valid');
end
switch tof_step
    case 1
        tof_step = days(hours(12)); % 12 hours in days
    case 2
        tof_step = days(days(1)); % 1 day
    case 3 
        tof_step = days(days(4)); % 4 days
    case 4
        tof_step = days(days(10)); % 10 days
    case 5 
        tof_step = days(days(20)); % 20 days
    otherwise 
        disp('Error: time step case not valid');
end

% Introducing degrees of freedom
t_dep_v = t_dep(1):t_dep_step:t_dep(2);
% TOFs in MJD2000
TOF1_v = MISSION.LEG1.TimeWindow.TOF1_min:tof_step:MISSION.LEG1.TimeWindow.TOF1_max;
TOF2_v = MISSION.LEG2.TimeWindow.TOF2_min:tof_step:MISSION.LEG2.TimeWindow.TOF2_max;

% Preallocation of the matrices
N1 = length(t_dep_v);
N2 = length(TOF1_v);
N3 = length(TOF2_v);

DELTA_V_TOT = NaN(N1, N2, N3);
DELTA_V1    = NaN(N1, N2, N3);
DELTA_VFB   = NaN(N1, N2, N3);
DELTA_V2    = NaN(N1, N2, N3); 
rp_hyp = NaN(N1,N2,N3);

% 3D GRID SEARCH 
% Cycle for departure from Mars
for ii = 1:N1
    
    t_dep = t_dep_v(ii);
    
    % Mars State @t_dep
    [kep_M,~] = uplanet(t_dep, MISSION.LEG1.Data.dep_planet);
    % Conversion to degrees for incl,raan,omega,theta
    kep_M(3:6) = rad2deg(kep_M(3:6));
    % Creating state vector 
    [s_M,~] = kep2car(kep_M(1),kep_M(2),kep_M(3), ...
                      kep_M(4),kep_M(5),kep_M(6), ...
                      MISSION.LEG1.Data.mu_sun);
    r_M = s_M(1:3);
    v_M = s_M(4:6);
    
    % Cycle for TOF1
    for jj = 1:N2
        
        TOF1 = TOF1_v(jj);
        t_FB = t_dep + TOF1;

        % Earth State @flyby
        [kep_E,~] = uplanet(t_FB, MISSION.FLYBY.Data.flyby_planet);
        % Conversion to degrees for incl,raan,omega,theta
        kep_E(3:6) = rad2deg(kep_E(3:6));
        % Creating state vector 
        [s_E,~] = kep2car(kep_E(1),kep_E(2),kep_E(3), ...
                          kep_E(4),kep_E(5),kep_E(6), ...
                          MISSION.LEG1.Data.mu_sun);
        r_E = s_E(1:3);
        v_E = s_E(4:6);
        
        % Converting TOF in seconds
        deltaT1 = seconds(days(TOF1));
        % Solve Lambert's arc (LEG1) for t_dep = t_dep_v(ii), TOF1 = TOF1_v(jj),
        % PROGRADE ORBIT, tau = 0
        % Solver takes by itself into account if deltaT < TPAR -> NO
        % SOLUTION
        try
            [~,~,~,~,v_M_t,v_E_t1,~,~] = lambertMR(r_M,r_E,deltaT1,MISSION.LEG1.Data.mu_sun,0);
        catch
            continue
        end
        v_M_t = v_M_t';
        v_E_t1 = v_E_t1';
        % Compute DeltaV1
        deltaV1 = norm(v_M_t - v_M);

        % CONSTRAINT: v_launcher
        if(deltaV1 < MISSION.LEG1.Data.v_launch)
            
            % Cycle for TOF2
            for kk = 1:N3
                TOF2 = TOF2_v(kk);
                t_arr = t_FB + TOF2;
                % Asteroid state @arrival
                [kep_A, ~, ~] = ephAsteroids(t_arr, MISSION.LEG2.Data.target_asteroid);
                % Conversion to degrees for incl,raan,omega,theta
                kep_A(3:6) = rad2deg(kep_A(3:6));
                % Creating state vector 
                [s_A,~] = kep2car(kep_A(1),kep_A(2),kep_A(3),kep_A(4),kep_A(5),kep_A(6),MISSION.LEG2.Data.mu_sun);
                % Extract position and velocity vector at the final orbit (arrival)
                r_A = s_A(1:3);     
                v_A = s_A(4:6);
                
                % Converting TOF in seconds
                deltaT2 = seconds(days(TOF2));
                % Solve Lambert's arc (LEG2) for t_dep = t_dep_v(ii), TOF1 = TOF1_v(jj),
                % TOF2 = TOF2_v(kk), PROGRADE ORBIT tau = 0
                % Solver takes by itself into account if deltaT < TPAR -> NO
                % SOLUTION
                try
                    [~,~,~,~,v_E_t2,v_A_t,~,~] = lambertMR(r_E,r_A,deltaT2,MISSION.LEG2.Data.mu_sun,0);
                catch
                    continue
                end
                v_A_t = v_A_t';
                v_E_t2 = v_E_t2';
                % Compute DeltaV2
                deltaV2 = norm(v_A_t - v_A);

                % Characterisation of Flyby Hyperbola
                % Hyperbolic excess velocities
                vinf_min = v_E_t1 - v_E;
                vinf_plus = v_E_t2 - v_E;

                % Solving the flyby
                % CONSTRAINT : rp_min : automatically discarded by find_rp
                % function into flyby_powered function
                FLYBY = flyby_powered(vinf_min, vinf_plus, ...
                                MISSION.FLYBY.Data.mu_planet, ...
                                MISSION.FLYBY.Data.R_planet, ...
                                MISSION.FLYBY.Data.h_min );
               
                if isnan(FLYBY.rp)
                    continue
                end

                rp = FLYBY.rp;
                deltaV_FB = FLYBY.Deltavp;
                DELTA_V1(ii,jj,kk)  = deltaV1;
                DELTA_VFB(ii,jj,kk)= deltaV_FB;
                DELTA_V2(ii,jj,kk) = deltaV2;
                rp_hyp(ii,jj,kk) = rp;
                DELTA_V_TOT(ii,jj,kk) = deltaV1 + deltaV_FB + deltaV2;
            end
        end
    end
end

% Checking Global Feasability
validMask = ~isnan(DELTA_V_TOT) & DELTA_V_TOT > 0;
if ~any(validMask, 'all')
    MISSION.Status.Feasible = false;
    warning(['No feasible transfer found within the provided time windows ' ...
             'given the constraints.']);
    return; % Exit
else
    MISSION.Status.Feasible = true;
end

% Computing outputs
MISSION.LEG1.Solution.DELTA_V1 = DELTA_V1;
MISSION.LEG1.Solution.t_dep_v = t_dep_v;
MISSION.LEG1.Solution.TOF1_v  = TOF1_v;

MISSION.LEG2.Solution.DELTA_V2 = DELTA_V2;
MISSION.LEG2.Solution.TOF2_v  = TOF2_v;

MISSION.FLYBY.Solution.DELTA_VFB = DELTA_VFB;
MISSION.FLYBY.Solution.rp = rp_hyp;

MISSION.DELTAV_TOT = DELTA_V_TOT;

% Computing Most Efficient ONLY FOR DEBUGGING AND CHECK 
% IMPORTANT : THIS IS NOT THE GLOBAL MINIMA
validMask = DELTA_V_TOT > 0 & ~isnan(DELTA_V_TOT);
[minDV, idx] = min(DELTA_V_TOT(validMask));
if ~isempty(idx)
    [i,j,k] = ind2sub(size(DELTA_V_TOT), find(DELTA_V_TOT == minDV,1));
    MISSION.MostEfficient.t_dep = t_dep_v(i);
    MISSION.MostEfficient.TOF1  = TOF1_v(j);
    MISSION.MostEfficient.TOF2  = TOF2_v(k);
    MISSION.MostEfficient.DELTA_V1   = DELTA_V1(i,j,k);
    MISSION.MostEfficient.DELTA_VFB  = DELTA_VFB(i,j,k);
    MISSION.MostEfficient.DELTA_V2   = DELTA_V2(i,j,k);
    MISSION.MostEfficient.DELTA_V_TOT = minDV;
end