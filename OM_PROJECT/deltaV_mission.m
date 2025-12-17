function [MISSION] = deltaV_mission(MISSION,time_step_case)
%DELTAV_MISSION This function computes the deltaV required for a certain
% specified mission. The Injection Orbit from departure planet and
% Insertion Orbit to arrival planet ARE NOT CONSIDERED.
%
% INPUTS : MISSION : This struct contains all the informations for the
%                   mission, and should be fullfilled with the following
%                   data: 
%                   DATA : All info about planetary constants, ID of
%                   planets and V_launcher constraint 
%                   TIMEWINDOW : Info about departure time window and
%                   arrival time window, specified in datetime format and
%                   vector: (yyyy,mm,d,h,m,s)
%          time_step_case : vector for time step case 
%                           1 : each time step = 1 minute (not recommended)
%                           2 : each time step = 15 minutes (not recommended)
%                           3 : each time step = 30 minutes 
%                           4 : each time step = 1 hour 
%                           5 : each time step = 1 day
%                           6 : each time step = 7 day
%                           7 : each time step = 14 days

% OUPTUTS : MISSION : The struct is returned with the added solution of the
%                     Lambert problem with the 2 degrees of freedom of
%                     times: 
%           MISSION.Solution.DELTAV_TOT : matrix with the total cost of the mission
%                                         deltaVTOT [km/s]
%           MISSION.Solution.DELTAV1 : matrix with the cost of initial
%                                      transfer manouvre
%                                         deltaV1 [km/s]
%           MISSION.Solution.DELTAV2 : matrix with the cost of final
%                                      transfer manouvre
%                                         deltaV2 [km/s]
%           MISSION.Solution.TOF : matrix with the time of flight for the
%                                  transfer arc
%                                         ToF [s]
%           MISSION.Solution.t1_v : vector of dimension Nx1 where N depends
%                                   on the time_step_case specified in
%                                   input, represents the time vector for
%                                   orbit 1 in MJD2000 days
%           MISSION.Solution.t1_v : vector of dimension Nx1 where N depends
%                                   on the time_step_case specified in
%                                   input, represents the time vector for
%                                   orbit 1 in MJD2000 days
%           MISSION.MostEfficient : This field contains the same
%                                   informations as the field solution, 
%                                   but for the most efficient
%                                   (in terms of Δv) transfer

if(nargin<2)
    disp('Error: unspecified time step case');
end

% Departure times in MJD2000 format
t1 = [date2mjd2000(MISSION.TimeWindow.Departure_v(1,:)), date2mjd2000(MISSION.TimeWindow.Departure_v(2,:))];
t2 = [date2mjd2000(MISSION.TimeWindow.Arrival_v(1,:)), date2mjd2000(MISSION.TimeWindow.Arrival_v(2,:))];

% Discretization of times. MJD2000 format
switch time_step_case
    case 1
        t_step = days(minutes(1)); % 1 minute in days
        t1_v = t1(1):t_step:t1(2);
        t2_v = t2(1):t_step:t2(2);
    case 2
        t_step = days(minutes(15)); % 15 minutes in days
        t1_v = t1(1):t_step:t1(2);
        t2_v = t2(1):t_step:t2(2);
    case 3 
        t_step = days(minutes(30)); % 30 minutes in days
        t1_v = t1(1):t_step:t1(2);
        t2_v = t2(1):t_step:t2(2);
    case 4
        t_step = days(hours(1)); % 1 hour in days
        t1_v = t1(1):t_step:t1(2);
        t2_v = t2(1):t_step:t2(2);
    case 5 
        t_step = days(days(1)); % 1 day
        t1_v = t1(1):t_step:t1(2);
        t2_v = t2(1):t_step:t2(2);
    case 6 
        t_step = days(days(7)); % 7 days
        t1_v = t1(1):t_step:t1(2);
        t2_v = t2(1):t_step:t2(2);
    case 7 
        t_step = days(days(14)); % 14 days
        t1_v = t1(1):t_step:t1(2);
        t2_v = t2(1):t_step:t2(2);
    otherwise 
        disp('Error: time step case not valid');
end

% Pre-allocate matrices
DELTAV_TOT = zeros(length(t1_v),length(t2_v));
TOF = zeros(length(t1_v),length(t2_v));
DELTAV1 = zeros(length(t1_v),length(t2_v));
DELTAV2 = zeros(length(t1_v),length(t2_v));

% Cycle for solution 
for ii = 1:length(t1_v)

    % Compute ephemerides of the first orbit (departure), not dependant on
    % jj, depending on time t1
    [eph1,~] = uplanet(t1_v(ii), MISSION.Data.dep_planet);
    % Converting radiants in degrees
    eph1 = [eph1(1) eph1(2) rad2deg(eph1(3)) rad2deg(eph1(4)) rad2deg(eph1(5)) rad2deg(eph1(6)) ];
    % Computing cartesian state of orbit 1
    [s1,~] = kep2car(eph1(1),eph1(2),eph1(3),eph1(4),eph1(5),eph1(6),MISSION.Data.mu_sun);
    % Extract position and velocity vector at the initial orbit (departure)
    r_t1 = s1(1:3);
    v_t1 = s1(4:6);
    
    for jj = 1:length(t2_v)
        if(t2_v(jj) > t1_v(ii)) % Checking that arrival date is AFTER departure date
        
        % Compute ephemerides of the second orbit (arrival), depending on
        % time t2
        [eph2,~] = uplanet(t2_v(jj), MISSION.Data.target_planet);
        % Converting radiants in degrees
        eph2 = [eph2(1) eph2(2) rad2deg(eph2(3)) rad2deg(eph2(4)) rad2deg(eph2(5)) rad2deg(eph2(6)) ];
        % Computing cartesian state of orbit 2
        [s2,~] = kep2car(eph2(1),eph2(2),eph2(3),eph2(4),eph2(5),eph2(6),MISSION.Data.mu_sun);
       % Extract position and velocity vector at the final orbit (arrival)
        r_t2 = s2(1:3);     
        v_t2 = s2(4:6);
        
        % Computing Time of Flight (ToF) 
        deltaT = seconds(days(t2_v(jj) - t1_v(ii)));
        
        % Solve Lambert's problem for t1 = t1_v(ii), t2 = t2_v(jj),
        % PROGRADE ORBIT
        % Solver takes by itself into account if deltaT < TPAR -> NO
        % SOLUTION
        [~,~,~,~,v_t1_t,v_t2_t,~,~] = lambertMR(r_t1,r_t2,deltaT,MISSION.Data.mu_sun,0);
        v_t1_t = v_t1_t';
        v_t2_t = v_t2_t';

        % Compute deltaVs
        deltaV1 = v_t1_t - v_t1; 
        deltaV2 = v_t2 - v_t2_t;
        deltaV = norm(deltaV1) + norm(deltaV2);
        
        % Update matrices
        % Check feasability in terms of ||deltaV1|| <= V_launcher
        if(norm(deltaV1) <= MISSION.Data.v_launch)
            DELTAV_TOT(ii,jj) = deltaV;
            DELTAV1(ii,jj) = norm(deltaV1);
            DELTAV2(ii,jj) = norm(deltaV2);
            TOF(ii,jj) = deltaT; 
        else
            DELTAV_TOT(ii,jj) = NaN;
            DELTAV1(ii,jj) = NaN;
            DELTAV2(ii,jj) = NaN; 
            TOF(ii,jj) = NaN; 
        end
            
        else % Case of Arrival Time BEFORE of Departure Time 
        DELTAV_TOT(ii,jj) = NaN; 
        DELTAV1(ii,jj) = NaN;
        DELTAV2(ii,jj) = NaN;
        TOF(ii,jj) = NaN;
        end

    end

end

% Checking Feasability
validMask = ~isnan(DELTAV_TOT) & DELTAV_TOT > 0;
if ~any(validMask, 'all')
    MISSION.Status.Feasible = false;
    warning(['No feasible transfer found within the provided time windows ' ...
             'given the launch constraint (||deltaV1|| <= ', num2str(MISSION.Data.v_launch), ' km/s).']);
    return; % Exit
else
    MISSION.Status.Feasible = true;
end


% Computing outputs
MISSION.Solution.DELTAV_TOT = DELTAV_TOT;
MISSION.Solution.DELTAV1 = DELTAV1;
MISSION.Solution.DELTAV2 = DELTAV2;
MISSION.Solution.TOF = TOF;
MISSION.Solution.t1_v = t1_v; 
MISSION.Solution.t2_v = t2_v;

% Finding most efficient transfer in terms of deltaV_TOT
deltaV_TOT_min = min(MISSION.Solution.DELTAV_TOT(:));
[idx_i, idx_j] = find(MISSION.Solution.DELTAV_TOT == deltaV_TOT_min);
if numel(idx_i)>1 || numel(idx_j)>1
    idx_i = idx_i(1); % Take the first in case of multiple correspondances
    idx_j = idx_j(1);
end
MISSION.MostEfficient.deltaV_TOT = MISSION.Solution.DELTAV_TOT(idx_i,idx_j);
MISSION.MostEfficient.deltaV1 = MISSION.Solution.DELTAV1(idx_i,idx_j);
MISSION.MostEfficient.deltaV2 = MISSION.Solution.DELTAV2(idx_i,idx_j);
MISSION.MostEfficient.TOF = MISSION.Solution.TOF(idx_i,idx_j);
MISSION.MostEfficient.t1 = MISSION.Solution.t1_v(idx_i); %MJD2000
MISSION.MostEfficient.t2 = MISSION.Solution.t2_v (idx_j); %MJD2000
MISSION.MostEfficient.Departure = datetime( mjd20002date(MISSION.MostEfficient.t1) ); % Date
MISSION.MostEfficient.Arrival = datetime( mjd20002date(MISSION.MostEfficient.t2) ); % Date
% Computing also ephemerides of orbit 1, orbit 2 at the MOST EFFICIENT
% TRANSFER for orbit propagation
[MISSION.MostEfficient.Eph_dep_planet,~] = uplanet(MISSION.MostEfficient.t1,...
                                                   MISSION.Data.dep_planet);
[MISSION.MostEfficient.Eph_target_planet,~] = uplanet(MISSION.MostEfficient.t2,...
                                                   MISSION.Data.target_planet);

end