%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%                                                         %%%%%%%%%
%%%%%%%%%                  ORBITAL MECHANICS PROJECT              %%%%%%%%%
%%%%%%%%%                       A.A. 2025-2026                    %%%%%%%%%
%%%%%%%%%                                                         %%%%%%%%%
%%%%%%%%%                        GROUP 2519                       %%%%%%%%%
%%%%%%%%%                       ASSIGNMENT 1                      %%%%%%%%%
%%%%%%%%%                                                         %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CONTRIBUTORS
% Amura Fabio, 10830618, fabio.amura@mail.polimi.it
% Crisanti Costanza, 10911209, costanza.crisanti@mail.polimi.it
% Deli Luca, 10907638, luca.deli@mail.polimi.it
% Tomas Fadista, 11027292, tomas.nascimentodasilva@mail.polimi.it

% REQUIREMENTS : Assignment 1
% Asteroid N.470864 
% Departure date : 00:00:00 01/01/2030 
% Arrival date: 00:00:00 01/01/2060 
% Departure planet : Mars
% Flyby planet : Earth

clear
close all
clc
format longE

%% --------------- ASSIGNMENT 1 ---------------

fprintf('\n----PRELIMINARY ANALYSIS FOR THE INTERPLANETARY TRANSFER MISSION DESIGN----\n')
% The preliminary analysis is needed in order to touch approximately the
% problem and understand what could be some reasonable values in terms of physical quantities
% like TOF, deltaVs,...

% --------- PRELIMINARY: LEG 1 --------
% Defining a struct containing all the data of the LEG1
PRE_LEG1.Data.mu_sun = astroConstants(4); %[km^3/s^2]
PRE_LEG1.Data.dep_planet = 4; % uplanet ID Mars
PRE_LEG1.Data.target_planet = 3; % uplanet ID Earth
% ARBITRARLY CHOOSE THE DATES FOR NOW
PRE_LEG1.TimeWindow.Departure = [datetime(2030,1,1,0,0,0); datetime(2060,1,1,0,0,0)];
PRE_LEG1.TimeWindow.Departure_v = [2030,1,1,0,0,0; 2060,1,1,0,0,0];
PRE_LEG1.TimeWindow.Arrival = [datetime(2030,1,1,0,0,0); datetime(2060,1,1,0,0,0)];
PRE_LEG1.TimeWindow.Arrival_v = [2030,1,1,0,0,0; 2060,1,1,0,0,0];
% Eventual constraint on v_launch : UNCONSTRAINED
PRE_LEG1.Data.v_launch = inf; %[km/s]

% Computing the DeltaV of the first leg 
time_step_case = 7; % 1 time step each 14 days
[PRE_LEG1] = deltaV_mission(PRE_LEG1,time_step_case);
% Refine solution using fmincon : initial guess is the one contained in the
% MostEfficient field of LEG1
x0 = [PRE_LEG1.MostEfficient.t1, PRE_LEG1.MostEfficient.t2];
[PRE_LEG1, ~, ~, ~, ~] = deltaV_mission_fmincon(PRE_LEG1, x0);

% Minimum deltaV transfer for LEG1
disp(PRE_LEG1.MostEfficient);

% ------- PRELIMINARY : LEG 2 ----------
% Defining a struct containing all the data of the LEG1
PRE_LEG2.Data.mu_sun = astroConstants(4); %[km^3/s^2]
PRE_LEG2.Data.dep_planet = 3; % uplanet ID Earth
PRE_LEG2.Data.target_asteroid = 470864; % uplanet ID Earth
% ARBITRARLY CHOOSE THE DATES FOR NOW
PRE_LEG2.TimeWindow.Departure = [datetime(2030,1,1,0,0,0); datetime(2060,1,1,0,0,0)];
PRE_LEG2.TimeWindow.Departure_v = [2030,1,1,0,0,0; 2060,1,1,0,0,0];
PRE_LEG2.TimeWindow.Arrival = [datetime(2030,1,1,0,0,0); datetime(2060,1,1,0,0,0)];
PRE_LEG2.TimeWindow.Arrival_v = [2030,1,1,0,0,0; 2060,1,1,0,0,0];
% Eventual constraint on v_launch : UNCONSTRAINED
PRE_LEG2.Data.v_launch = inf; %[km/s]

% Computing the DeltaV of the first leg 
time_step_case = 7; % 1 time step each 14 days
[PRE_LEG2] = deltaV_mission_asteroid(PRE_LEG2,time_step_case);
% Refine solution using fmincon : initial guess is the one contained in the
% MostEfficient field of LEG1
x0 = [PRE_LEG2.MostEfficient.t1, PRE_LEG2.MostEfficient.t2];
[PRE_LEG2, ~, ~, ~, ~] = deltaV_mission_fmincon_asteroid(PRE_LEG2, x0);

% Minimum deltaV transfer for LEG
disp(PRE_LEG2.MostEfficient);

%% -------- 2030-2040 -------------

% ---------- MISSION DESIGN - 2030-2040 ------------

fprintf('\n----INTERPLANETARY TRANSFER MISSION DESIGN 2030-2040----\n')
% The initial step is to consider a very 'rough' 3D grid search algorithm.
% Defining the struct MISSION1 containing all the data 
MISSION1.LEG1.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION1.LEG1.Data.dep_planet = 4; % uplanet ID Mars
MISSION1.LEG1.Data.target_planet = 3; % uplanet ID Earth
MISSION1.LEG2.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION1.LEG2.Data.dep_planet = 3; % uplanet ID Earth
MISSION1.LEG2.Data.target_asteroid = 470864; % ephAsteroid ID
MISSION1.FLYBY.Data.flyby_planet =  3; % uplanet ID Earth
MISSION1.FLYBY.Data.R_planet = 6378; % Equatorial radius of Earth, [km]
MISSION1.FLYBY.Data.mu_planet = astroConstants(13); % [km^3/s^2]

% Time window : ONLY FOR MARS DEPARTURE
% First 10 years, for computational efficiency, than repeat the analysis
% for the other 10-10 years
MISSION1.LEG1.TimeWindow.Departure = [datetime(2030,1,1,0,0,0); datetime(2040,1,1,0,0,0)];
MISSION1.LEG1.TimeWindow.Departure_v = [2030,1,1,0,0,0; 2040,1,1,0,0,0];
% The time window for the departure from Mars identifies the first degree
% of freedom of our mission design 
% For the 3D grid search, we have (t_dep, TOF1, TOF2) as degrees of
% freedom, and this is useful because reduces drastically the amount of
% time and calculations to do at each cycle with respet to the case of
% having (t_dep, t_fb, t_arr)
% We have to take a guess for the ranges of the TOFs, and we can use for
% example this:
MISSION1.LEG1.TimeWindow.TOF1_min = days(seconds(PRE_LEG1.MostEfficient.TOF/2));
MISSION1.LEG1.TimeWindow.TOF1_max = days(seconds(PRE_LEG1.MostEfficient.TOF*2));
% And also: 
MISSION1.LEG2.TimeWindow.TOF2_min = days(seconds(PRE_LEG2.MostEfficient.TOF/2));
MISSION1.LEG2.TimeWindow.TOF2_max = days(seconds(PRE_LEG2.MostEfficient.TOF*2));

% Setting constraints on the v_launcher, altitude of flyby, ...
MISSION1.LEG1.Data.v_launch = 5; %[km/s]
MISSION1.FLYBY.Data.h_min = 400; % Minimum altitude of flyby hyperbola's pericentre[km]  
...

% 3D Grid search for first estimation of possible location of the minima 
% t_dep_step : 30 days (4)
% tof_step : 20 days (5)
MISSION1 = deltaV_mission_3B(MISSION1,4,5)

% Display the minimum found for later considerations
fprintf('\nMinimum Δv_TOT for MISSION1 (coarse 3D Grid Search Algorithm): %.4f km/s\n',...
    MISSION1.MostEfficient.DELTA_V_TOT);
fprintf('Total TOF: \nTOF1: %d days\n',MISSION1.MostEfficient.TOF1);
fprintf('TOF2: %d days\n',...
    MISSION1.MostEfficient.TOF2);

% Now the solution for the first 10 years is in terms of
% DELTAV_TOT(t_dep,TOF1,TOF2). 
% We can use the porkchop plots tool to identify more precisely sub-optimal
% regions in which the refinement can be done in order to try to catch more
% precisely the local minima found by the coarse 3D grid search.
% In order to do that, we have to let one of the three dimensions 'collapse'
% (by using the min function).
% In other words, we are freezing one of the three variables in order to
% make a 2D representation (porkchop plots)

% Struct PORKCHOP1_1 (associated to MISSION1 for first 10 years, and leg 1
% driven)
% Porkchop LEG1 driven, x axis : t_dep, y axis : TOF1
% Physical meaning : whatever will be the arrival, what are the best
% regions for the departure and the flyby (t_FB = t_dep + TOF1)?
PORKCHOP1_1.Solution.DELTAV_TOT = min(MISSION1.DELTAV_TOT, [], 3); % [N_dep x N_TOF1]
PORKCHOP1_1.Solution.t1_v = MISSION1.LEG1.Solution.t_dep_v;
PORKCHOP1_1.Solution.t2_v = MISSION1.LEG1.Solution.TOF1_v;

porkchop_date_tof(PORKCHOP1_1)
fig = gcf;
set(fig, 'Name', 'GRID SEARCH 2030-2040, LEG1 driven', 'NumberTitle', 'off');


% Porkchop LEG2 driven , x axis : TOF1, y axis : TOF2
% Physical meaning: whatever will be the departure, what are the best
% regions for the flyby and the arrival (t_arr = t_FB + TOF2)?
PORKCHOP1_2.Solution.DELTAV_TOT = squeeze(min(MISSION1.DELTAV_TOT, [], 1)); % [N_TOF1 x N_TOF2]
PORKCHOP1_2.Solution.t1_v = MISSION1.LEG1.Solution.TOF1_v;
PORKCHOP1_2.Solution.t2_v = MISSION1.LEG2.Solution.TOF2_v;

porkchop_tof_tof(PORKCHOP1_2)
fig = gcf;
set(fig, 'Name', 'GRID SEARCH 2030-2040, LEG2 driven', 'NumberTitle', 'off');

% From the first porkchop two sub-optimal regions are worth for further
% investigations: 1) t_dep : [2036,4,1,0,0,0; 2038,1,1,0,0,0]
%                 2) t_dep : [2039,1,1,0,0,0; 2040,1,1,0,0,0]
% Instead, from the second porkchop, we can see that limiting the TOF1
% would be dangerous in the sense that it could cancel solutions. Even for
% the TOF2 we can't retrieve nice limits.
% So we can assume : 1) TOF1 : [160, 570]
%                    2) TOF2 : [200, 600]

% ---------- MISSION REFINEMENT - 2030-2040 ------------

fprintf('\n----INTERPLANETARY TRANSFER MISSION REFINEMENT 2030-2040----\n')

% In this phase, always considering only the first 10 years, we can refine
% the solution before using the optimization tool in order to have the best
% possible initial guess according to our discretization and analysis

% Defining the struct MISSION1_REF1 for the refinement 1
% MISSION1_REF1 refines the first interesting good region in terms of
% departure, MISSION1_REF2 refines instead the second one.
% The TOFs are fixed for both the refinements
MISSION1_REF1.LEG1.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION1_REF1.LEG1.Data.dep_planet = 4; % uplanet ID Mars
MISSION1_REF1.LEG1.Data.target_planet = 3; % uplanet ID Earth
MISSION1_REF1.LEG2.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION1_REF1.LEG2.Data.dep_planet = 3; % uplanet ID Earth
MISSION1_REF1.LEG2.Data.target_asteroid = 470864; % ephAsteroid ID
MISSION1_REF1.FLYBY.Data.flyby_planet =  3; % uplanet ID Earth
MISSION1_REF1.FLYBY.Data.R_planet = 6378; % Equatorial radius of Earth, [km]
MISSION1_REF1.FLYBY.Data.mu_planet = astroConstants(13); % [km^3/s^2]

% Time window : ONLY MARS DEPARTURE
MISSION1_REF1.LEG1.TimeWindow.Departure = [datetime(2036,4,1,0,0,0); datetime(2038,1,1,0,0,0)];
MISSION1_REF1.LEG1.TimeWindow.Departure_v = [2036,4,1,0,0,0; 2038,1,1,0,0,0];
MISSION1_REF1.LEG1.TimeWindow.TOF1_min = 160;
MISSION1_REF1.LEG1.TimeWindow.TOF1_max = 570;
% And also: 
MISSION1_REF1.LEG2.TimeWindow.TOF2_min = 200;
MISSION1_REF1.LEG2.TimeWindow.TOF2_max = 600;

% Setting constraints on the v_launcher, altitude of flyby, ...
MISSION1_REF1.LEG1.Data.v_launch = 5; %[km/s]
MISSION1_REF1.FLYBY.Data.h_min = 400; % Minimum altitude of flyby hyperbola's pericentre[km]  
...

% 3D Grid search for first estimation of possible location of the minima
% t_dep_step : 7 days (2)
% tof_step : 4 days (3)
MISSION1_REF1 = deltaV_mission_3B(MISSION1_REF1,2,3)

% Display the minimum found
fprintf('\nMinimum Δv_TOT for MISSION1_REF1 (refinement 1, first good region): %.4f km/s\n',...
    MISSION1_REF1.MostEfficient.DELTA_V_TOT);
fprintf(['Note : this local minimum is less promising than the grid search one, also in terms of ' ...
    ' Total TOF.\nTOF1: %d days\n'],MISSION1_REF1.MostEfficient.TOF1);
fprintf('TOF2: %d days\n',...
    MISSION1_REF1.MostEfficient.TOF2);


% In the same way of before we pass to the porkchop plot analysis:
% Struct PORKCHOP1_REF1_1 (associated to MISSION1 for first 10 years,
% refinement 1 for departure date and leg 1 driven)
% Porkchop LEG1 driven, x axis : t_dep, y axis : TOF1
% Physical meaning : whatever will be the arrival, what are the best
% regions for the departure and the flyby (t_FB = t_dep + TOF1)?
PORKCHOP1_REF1_1.Solution.DELTAV_TOT = min(MISSION1_REF1.DELTAV_TOT, [], 3); % [N_dep x N_TOF1]
PORKCHOP1_REF1_1.Solution.t1_v = MISSION1_REF1.LEG1.Solution.t_dep_v;
PORKCHOP1_REF1_1.Solution.t2_v = MISSION1_REF1.LEG1.Solution.TOF1_v;

porkchop_date_tof(PORKCHOP1_REF1_1)
fig = gcf;
set(fig, 'Name', 'REFINED 1ST REGION 2030-2040, LEG1 driven', 'NumberTitle', 'off');

% Porkchop LEG2 driven , x axis : TOF1, y axis : TOF2
% Physical meaning: whatever will be the departure, what are the best
% regions for the flyby and the arrival (t_arr = t_FB + TOF2)?
PORKCHOP1_REF1_2.Solution.DELTAV_TOT = squeeze(min(MISSION1_REF1.DELTAV_TOT, [], 1)); % [N_TOF1 x N_TOF2]
PORKCHOP1_REF1_2.Solution.t1_v = MISSION1_REF1.LEG1.Solution.TOF1_v;
PORKCHOP1_REF1_2.Solution.t2_v = MISSION1_REF1.LEG2.Solution.TOF2_v;

porkchop_tof_tof(PORKCHOP1_REF1_2)
fig = gcf;
set(fig, 'Name', 'REFINED 1ST REGION 2030-2040, LEG2 driven', 'NumberTitle', 'off');

% Defining the struct MISSION1_REF2 for the refinement 2
MISSION1_REF2.LEG1.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION1_REF2.LEG1.Data.dep_planet = 4; % uplanet ID Mars
MISSION1_REF2.LEG1.Data.target_planet = 3; % uplanet ID Earth
MISSION1_REF2.LEG2.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION1_REF2.LEG2.Data.dep_planet = 3; % uplanet ID Earth
MISSION1_REF2.LEG2.Data.target_asteroid = 470864; % ephAsteroid ID
MISSION1_REF2.FLYBY.Data.flyby_planet =  3; % uplanet ID Earth
MISSION1_REF2.FLYBY.Data.R_planet = 6378; % Equatorial radius of Earth, [km]
MISSION1_REF2.FLYBY.Data.mu_planet = astroConstants(13); % [km^3/s^2]

% Time window : ONLY MARS DEPARTURE
MISSION1_REF2.LEG1.TimeWindow.Departure = [datetime(2039,1,1,0,0,0); datetime(2040,1,1,0,0,0)];
MISSION1_REF2.LEG1.TimeWindow.Departure_v = [2039,1,1,0,0,0; 2040,1,1,0,0,0];
MISSION1_REF2.LEG1.TimeWindow.TOF1_min = 160;
MISSION1_REF2.LEG1.TimeWindow.TOF1_max = 570;
% And also: 
MISSION1_REF2.LEG2.TimeWindow.TOF2_min = 200;
MISSION1_REF2.LEG2.TimeWindow.TOF2_max = 600;

% Setting constraints on the v_launcher, altitude of flyby, ...
MISSION1_REF2.LEG1.Data.v_launch = 5; %[km/s]
MISSION1_REF2.FLYBY.Data.h_min = 400; % Minimum altitude of flyby hyperbola's pericentre[km]  
...

% 3D Grid search for first estimation of possible location of the minima 
% t_dep_step : 7 days (2)
% tof_step : 4 days (3)
MISSION1_REF2 = deltaV_mission_3B(MISSION1_REF2,2,3)

% Display the minimum found
fprintf('\nMinimum Δv_TOT for MISSION1_REF2 (refinement 2, second good region): %.4f km/s\n',...
    MISSION1_REF2.MostEfficient.DELTA_V_TOT);
fprintf(['Note : this local minimum is promising in terms of Δv_TOT. ' ...
    '\nTotal TOF:\nTOF1: %d days\n'],MISSION1_REF2.MostEfficient.TOF1);
fprintf('TOF2: %d days\n',...
    MISSION1_REF2.MostEfficient.TOF2);


% Struct PORKCHOP1_REF2_1 (associated to MISSION1 for first 10 years,
% refinement 2 for departure date and leg 1 driven)
% Porkchop LEG1 driven, x axis : t_dep, y axis : TOF1
% Physical meaning : whatever will be the arrival, what are the best
% regions for the departure and the flyby (t_FB = t_dep + TOF1)?
PORKCHOP1_REF2_1.Solution.DELTAV_TOT = min(MISSION1_REF2.DELTAV_TOT, [], 3); % [N_dep x N_TOF1]
PORKCHOP1_REF2_1.Solution.t1_v = MISSION1_REF2.LEG1.Solution.t_dep_v;
PORKCHOP1_REF2_1.Solution.t2_v = MISSION1_REF2.LEG1.Solution.TOF1_v;

porkchop_date_tof(PORKCHOP1_REF2_1)
fig = gcf;
set(fig, 'Name', 'REFINED 2ND REGION 2030-2040, LEG1 driven', 'NumberTitle', 'off');

% Porkchop LEG2 driven , x axis : TOF1, y axis : TOF2
% Physical meaning: whatever will be the departure, what are the best
% regions for the flyby and the arrival (t_arr = t_FB + TOF2)?
PORKCHOP1_REF2_2.Solution.DELTAV_TOT = squeeze(min(MISSION1_REF2.DELTAV_TOT, [], 1)); % [N_TOF1 x N_TOF2]
PORKCHOP1_REF2_2.Solution.t1_v = MISSION1_REF2.LEG1.Solution.TOF1_v;
PORKCHOP1_REF2_2.Solution.t2_v = MISSION1_REF2.LEG2.Solution.TOF2_v;

porkchop_tof_tof(PORKCHOP1_REF2_2)
fig = gcf;
set(fig, 'Name', 'REFINED 2ND REGION 2030-2040, LEG2 driven', 'NumberTitle', 'off');

% ---------- MISSION OPTIMIZATION - 2030-2040 ------------

fprintf('\n----INTERPLANETARY TRANSFER MISSION OPTIMIZATION 2030-2040----\n')

% Now, we can finally do the optimization algorithm with fmincon taking
% always into account the constraints on v_launcher and on the minimum
% altitude for the powered flyby.
% Observing the results up to now, it's evident that the first region after
% the refinement could not be a possible optimal-minima. The second one,
% instead, could be considered for the optimization. 

% Defining the struct MISSION1_OPT  (only one case for these first 10
% years, the case of the refiniment 2)
MISSION1_OPT.LEG1.Data = MISSION1_REF2.LEG1.Data;
MISSION1_OPT.LEG2.Data = MISSION1_REF2.LEG2.Data;
MISSION1_OPT.FLYBY.Data = MISSION1_REF2.FLYBY.Data;
MISSION1_OPT.LEG1.TimeWindow = MISSION1_REF2.LEG1.TimeWindow;
MISSION1_OPT.LEG2.TimeWindow = MISSION1_REF2.LEG2.TimeWindow;

% Setting up the initial guess for the optimization tool to start from,
% chosen as x0 = [t_dep, TOF1, TOF2] in MJD2000 format (3 DOFs).
% The initial guess is found in MISSION1_REF2.MostEfficient field
x0 = [MISSION1_REF2.MostEfficient.t_dep, ...
      MISSION1_REF2.MostEfficient.TOF1,...
      MISSION1_REF2.MostEfficient.TOF2];
% Call optimization function
[MISSION1_OPT, ~, ~, ~, ~] = deltaV_mission_fmincon_3B(MISSION1_OPT, x0);

% Display Optimal solution found
disp(MISSION1_OPT.MostEfficient)

% The analysis for the first 10 years is now concluded, the very last thing
% to do is to store the deltaV for the final comparison (after having
% repeated the analysis for the other cases) in the struct DV_STORE, and
% then considering the next case
DV_STORE.OPT1 = MISSION1_OPT.MostEfficient;

%% -------- 2040-2050 -------------

% ---------- MISSION DESIGN - 2040-2050 ------------

fprintf('\n----INTERPLANETARY TRANSFER MISSION DESIGN 2040-2050----\n')

% Defining the struct MISSION2 containing all the data 
MISSION2.LEG1.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION2.LEG1.Data.dep_planet = 4; % uplanet ID Mars
MISSION2.LEG1.Data.target_planet = 3; % uplanet ID Earth
MISSION2.LEG2.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION2.LEG2.Data.dep_planet = 3; % uplanet ID Earth
MISSION2.LEG2.Data.target_asteroid = 470864; % ephAsteroid ID
MISSION2.FLYBY.Data.flyby_planet =  3; % uplanet ID Earth
MISSION2.FLYBY.Data.R_planet = 6378; % Equatorial radius of Earth, [km]
MISSION2.FLYBY.Data.mu_planet = astroConstants(13); % [km^3/s^2]

% Time window : ONLY FOR MARS DEPARTURE
% Consider now the window in the middle : 2040-2050
MISSION2.LEG1.TimeWindow.Departure = [datetime(2040,1,1,0,0,0); datetime(2050,1,1,0,0,0)];
MISSION2.LEG1.TimeWindow.Departure_v = [2040,1,1,0,0,0; 2050,1,1,0,0,0];
% The time window for the departure from Mars identifies the first degree
% of freedom of our mission design 
% For the 3D grid search, we have (t_dep, TOF1, TOF2) as degrees of
% freedom, and this is useful because reduces drastically the amount of
% time and calculations to do at each cycle with respet to the case of
% having (t_dep, t_fb, t_arr)
% The guess for the ranges of the TOFs is the same of the previous case,
% dictated by the preliminary analysis
MISSION2.LEG1.TimeWindow.TOF1_min = days(seconds(PRE_LEG1.MostEfficient.TOF/2));
MISSION2.LEG1.TimeWindow.TOF1_max = days(seconds(PRE_LEG1.MostEfficient.TOF*2));
% And also: 
MISSION2.LEG2.TimeWindow.TOF2_min = days(seconds(PRE_LEG2.MostEfficient.TOF/2));
MISSION2.LEG2.TimeWindow.TOF2_max = days(seconds(PRE_LEG2.MostEfficient.TOF*2));

% Setting constraints on the v_launcher, altitude of flyby, ...
MISSION2.LEG1.Data.v_launch = 5; %[km/s]
MISSION2.FLYBY.Data.h_min = 400; % Minimum altitude of flyby hyperbola's pericentre[km]  
...

% 3D Grid search for first estimation of possible location of the minima 
% t_dep_step : 30 days (4)
% tof_step : 20 days (5)
MISSION2 = deltaV_mission_3B(MISSION2,4,5)

% Display the minimum found for later considerations
fprintf('\nMinimum Δv_TOT for MISSION2 (coarse 3D Grid Search Algorithm): %.4f km/s\n',...
    MISSION2.MostEfficient.DELTA_V_TOT);
fprintf('Total TOF: \nTOF1: %d days\n',MISSION2.MostEfficient.TOF1);
fprintf('TOF2: %d days\n',...
    MISSION2.MostEfficient.TOF2);

% The solution is always a 3D-function in the form:
% DELTAV_TOT(t_dep,TOF1,TOF2). 
% We can use the porkchop plots tool to identify more precisely sub-optimal
% regions in which the refinement can be done in order to try to catch more
% precisely the local minima found by the coarse 3D grid search.
% In order to do that, we have to let one of the three dimensions 'collapse'
% (by using the min function).
% In other words, we are freezing one of the three variables in order to
% make a 2D representation (porkchop plots)

% Struct PORKCHOP2_1 (associated to MISSION2 for second 10 years, and leg 1
% driven)
% Porkchop LEG1 driven, x axis : t_dep, y axis : TOF1
% Physical meaning : whatever will be the arrival, what are the best
% regions for the departure and the flyby (t_FB = t_dep + TOF1)?
PORKCHOP2_1.Solution.DELTAV_TOT = min(MISSION2.DELTAV_TOT, [], 3); % [N_dep x N_TOF1]
PORKCHOP2_1.Solution.t1_v = MISSION2.LEG1.Solution.t_dep_v;
PORKCHOP2_1.Solution.t2_v = MISSION2.LEG1.Solution.TOF1_v;

porkchop_date_tof(PORKCHOP2_1)
fig = gcf;
set(fig, 'Name', 'GRID SEARCH 2040-2050, LEG1 driven', 'NumberTitle', 'off');

% Porkchop LEG2 driven , x axis : TOF1, y axis : TOF2
% Physical meaning: whatever will be the departure, what are the best
% regions for the flyby and the arrival (t_arr = t_FB + TOF2)?
PORKCHOP2_2.Solution.DELTAV_TOT = squeeze(min(MISSION2.DELTAV_TOT, [], 1)); % [N_TOF1 x N_TOF2]
PORKCHOP2_2.Solution.t1_v = MISSION2.LEG1.Solution.TOF1_v;
PORKCHOP2_2.Solution.t2_v = MISSION2.LEG2.Solution.TOF2_v;

porkchop_tof_tof(PORKCHOP2_2)
fig = gcf;
set(fig, 'Name', 'GRID SEARCH 2040-2050, LEG2 driven', 'NumberTitle', 'off');

% From the first porkchop, this time we can see a significant good region
% and a very borderline region to be considered good, but it will be anyway
% counted in for the refinement and maybe discarded after (if the local
% minima won't be satisfying).
% 1) t_dep : [2041,1,1,0,0,0; 2042,7,1,0,0,0]
% 2) t_dep : [2043,1,1,0,0,0; 2044,4,1,0,0,0]
% Instead, from the second porkchop, we can get useful information
% regarding the TOFs and see that we have a 'sub-optimal rectangle' that can
% be assumed with limits : 1) TOF1 : [160, 450]
%                          2) TOF2 : [200, 400]

% ---------- MISSION REFINEMENT - 2040-2050 ------------

fprintf('\n----INTERPLANETARY TRANSFER MISSION REFINEMENT 2040-2050----\n')

% In this phase, always considering the second 10 years, we can refine
% the solution before using the optimization tool in order to have the best
% possible initial guess according to our discretization and analysis

% Defining the struct MISSION1_REF1 for the refinement 1
% MISSION1_REF1 refines the first interesting good region in terms of
% departure, MISSION_REF2 refines instead the second one.
% The TOFs are fixed for both the refinements
MISSION2_REF1.LEG1.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION2_REF1.LEG1.Data.dep_planet = 4; % uplanet ID Mars
MISSION2_REF1.LEG1.Data.target_planet = 3; % uplanet ID Earth
MISSION2_REF1.LEG2.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION2_REF1.LEG2.Data.dep_planet = 3; % uplanet ID Earth
MISSION2_REF1.LEG2.Data.target_asteroid = 470864; % ephAsteroid ID
MISSION2_REF1.FLYBY.Data.flyby_planet =  3; % uplanet ID Earth
MISSION2_REF1.FLYBY.Data.R_planet = 6378; % Equatorial radius of Earth, [km]
MISSION2_REF1.FLYBY.Data.mu_planet = astroConstants(13); % [km^3/s^2]

% Time window : ONLY MARS DEPARTURE
MISSION2_REF1.LEG1.TimeWindow.Departure = [datetime(2041,1,1,0,0,0); datetime(2042,7,1,0,0,0)];
MISSION2_REF1.LEG1.TimeWindow.Departure_v = [2041,1,1,0,0,0; 2042,7,1,0,0,0];
MISSION2_REF1.LEG1.TimeWindow.TOF1_min = 160;
MISSION2_REF1.LEG1.TimeWindow.TOF1_max = 450;
% And also: 
MISSION2_REF1.LEG2.TimeWindow.TOF2_min = 200;
MISSION2_REF1.LEG2.TimeWindow.TOF2_max = 400;

% Setting constraints on the v_launcher, altitude of flyby, ...
MISSION2_REF1.LEG1.Data.v_launch = 5; %[km/s]
MISSION2_REF1.FLYBY.Data.h_min = 400; % Minimum altitude of flyby hyperbola's pericentre[km]  
...

% 3D Grid search for first estimation of possible location of the minima
% t_dep_step : 7 days (2)
% tof_step :  days (3)
MISSION2_REF1 = deltaV_mission_3B(MISSION2_REF1,2,3)

% Display the minimum found
fprintf('\nMinimum Δv_TOT for MISSION2_REF1 (refinement 1, first good region): %.4f km/s\n',...
    MISSION2_REF1.MostEfficient.DELTA_V_TOT);
fprintf(['Note : this local minimum is very promising and is most likely to be the grid search one.' ...
    ' Total TOF:\nTOF1: %d days\n'],MISSION2_REF1.MostEfficient.TOF1);
fprintf('TOF2: %d days\n',...
    MISSION2_REF1.MostEfficient.TOF2);

% In the same way of before we pass to the porkchop plot analysis:
% Struct PORKCHOP2_REF1_1 (associated to MISSION2 for second 10 years,
% refinement 1 for departure date and leg 1 driven)
% Porkchop LEG1 driven, x axis : t_dep, y axis : TOF1
% Physical meaning : whatever will be the arrival, what are the best
% regions for the departure and the flyby (t_FB = t_dep + TOF1)?
PORKCHOP2_REF1_1.Solution.DELTAV_TOT = min(MISSION2_REF1.DELTAV_TOT, [], 3); % [N_dep x N_TOF1]
PORKCHOP2_REF1_1.Solution.t1_v = MISSION2_REF1.LEG1.Solution.t_dep_v;
PORKCHOP2_REF1_1.Solution.t2_v = MISSION2_REF1.LEG1.Solution.TOF1_v;

porkchop_date_tof(PORKCHOP2_REF1_1)
fig = gcf;
set(fig, 'Name', 'REFINED 1ST REGION 2040-2050, LEG1 driven', 'NumberTitle', 'off');

% Porkchop LEG2 driven , x axis : TOF1, y axis : TOF2
% Physical meaning: whatever will be the departure, what are the best
% regions for the flyby and the arrival (t_arr = t_FB + TOF2)?
PORKCHOP2_REF1_2.Solution.DELTAV_TOT = squeeze(min(MISSION2_REF1.DELTAV_TOT, [], 1)); % [N_TOF1 x N_TOF2]
PORKCHOP2_REF1_2.Solution.t1_v = MISSION2_REF1.LEG1.Solution.TOF1_v;
PORKCHOP2_REF1_2.Solution.t2_v = MISSION2_REF1.LEG2.Solution.TOF2_v;

porkchop_tof_tof(PORKCHOP2_REF1_2)
fig = gcf;
set(fig, 'Name', 'REFINED 1ST REGION 2040-2050, LEG2 driven', 'NumberTitle', 'off');


% Defining the struct MISSION1_REF2 for the refinement 2
MISSION2_REF2.LEG1.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION2_REF2.LEG1.Data.dep_planet = 4; % uplanet ID Mars
MISSION2_REF2.LEG1.Data.target_planet = 3; % uplanet ID Earth
MISSION2_REF2.LEG2.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION2_REF2.LEG2.Data.dep_planet = 3; % uplanet ID Earth
MISSION2_REF2.LEG2.Data.target_asteroid = 470864; % ephAsteroid ID
MISSION2_REF2.FLYBY.Data.flyby_planet =  3; % uplanet ID Earth
MISSION2_REF2.FLYBY.Data.R_planet = 6378; % Equatorial radius of Earth, [km]
MISSION2_REF2.FLYBY.Data.mu_planet = astroConstants(13); % [km^3/s^2]

% Time window : ONLY MARS DEPARTURE
MISSION2_REF2.LEG1.TimeWindow.Departure = [datetime(2043,1,1,0,0,0); datetime(2044,4,1,0,0,0)];
MISSION2_REF2.LEG1.TimeWindow.Departure_v = [2043,1,1,0,0,0; 2044,4,1,0,0,0];
MISSION2_REF2.LEG1.TimeWindow.TOF1_min = 160;
MISSION2_REF2.LEG1.TimeWindow.TOF1_max = 450;
% And also: 
MISSION2_REF2.LEG2.TimeWindow.TOF2_min = 200;
MISSION2_REF2.LEG2.TimeWindow.TOF2_max = 400;

% Setting constraints on the v_launcher, altitude of flyby, ...
MISSION2_REF2.LEG1.Data.v_launch = 5; %[km/s]
MISSION2_REF2.FLYBY.Data.h_min = 400; % Minimum altitude of flyby hyperbola's pericentre[km]  
...

% 3D Grid search for first estimation of possible location of the minima 
% t_dep_step : 7 days (2)
% tof_step : 4 days (3)
MISSION2_REF2 = deltaV_mission_3B(MISSION2_REF2,2,3)

% Display the minimum found
fprintf('\nMinimum Δv_TOT for MISSION2_REF2 (refinement 2, second good region): %.4f km/s\n',...
    MISSION2_REF2.MostEfficient.DELTA_V_TOT);
fprintf(['Note : this local minimum is not promising at all, it wont be taken into account for the optimization.' ...
    '\nTotal TOF:\nTOF1: %d days\n'],MISSION2_REF2.MostEfficient.TOF1);
fprintf('TOF2: %d days\n',...
    MISSION2_REF2.MostEfficient.TOF2);


% Struct PORKCHOP2_REF2_1 (associated to MISSION2 for second 10 years,
% refinement 2 for departure date and leg 1 driven)
% Porkchop LEG1 driven, x axis : t_dep, y axis : TOF1
% Physical meaning : whatever will be the arrival, what are the best
% regions for the departure and the flyby (t_FB = t_dep + TOF1)?
PORKCHOP2_REF2_1.Solution.DELTAV_TOT = min(MISSION2_REF2.DELTAV_TOT, [], 3); % [N_dep x N_TOF1]
PORKCHOP2_REF2_1.Solution.t1_v = MISSION2_REF2.LEG1.Solution.t_dep_v;
PORKCHOP2_REF2_1.Solution.t2_v = MISSION2_REF2.LEG1.Solution.TOF1_v;

porkchop_date_tof(PORKCHOP2_REF2_1)
fig = gcf;
set(fig, 'Name', 'REFINED 2ND REGION 2040-2050, LEG1 driven', 'NumberTitle', 'off');

% Porkchop LEG2 driven , x axis : TOF1, y axis : TOF2
% Physical meaning: whatever will be the departure, what are the best
% regions for the flyby and the arrival (t_arr = t_FB + TOF2)?
PORKCHOP2_REF2_2.Solution.DELTAV_TOT = squeeze(min(MISSION2_REF2.DELTAV_TOT, [], 1)); % [N_TOF1 x N_TOF2]
PORKCHOP2_REF2_2.Solution.t1_v = MISSION2_REF2.LEG1.Solution.TOF1_v;
PORKCHOP2_REF2_2.Solution.t2_v = MISSION2_REF2.LEG2.Solution.TOF2_v;

porkchop_tof_tof(PORKCHOP2_REF2_2)
fig = gcf;
set(fig, 'Name', 'REFINED 2ND REGION 2040-2050, LEG2 driven', 'NumberTitle', 'off');


% ---------- MISSION OPTIMIZATION - 2040-2050 ------------

fprintf('\n----INTERPLANETARY TRANSFER MISSION OPTIMIZATION 2040-2050----\n')

% The optimization algorithm with fmincon is considered taking
% always into account the constraints on v_launcher and on the minimum
% altitude for the powered flyby.
% We will take into account only the first region refined for this second case, which led to 
% the most promising minima.

% Defining the struct MISSION2_OPT
MISSION2_OPT.LEG1.Data = MISSION2_REF1.LEG1.Data;
MISSION2_OPT.LEG2.Data = MISSION2_REF1.LEG2.Data;
MISSION2_OPT.FLYBY.Data = MISSION2_REF1.FLYBY.Data;
MISSION2_OPT.LEG1.TimeWindow = MISSION2_REF1.LEG1.TimeWindow;
MISSION2_OPT.LEG2.TimeWindow = MISSION2_REF1.LEG2.TimeWindow;

% Setting up the initial guess for the optimization tool to start from,
% chosen as x0 = [t_dep, TOF1, TOF2] in MJD2000 format (3 DOFs).
% The initial guess is found in MISSION2_REF1.MostEfficient field
x0 = [MISSION2_REF1.MostEfficient.t_dep, ...
      MISSION2_REF1.MostEfficient.TOF1,...
      MISSION2_REF1.MostEfficient.TOF2];

% Call optimization function
[MISSION2_OPT, ~, ~, ~, ~] = deltaV_mission_fmincon_3B(MISSION2_OPT, x0);

% Display Optimal solution found
disp(MISSION2_OPT.MostEfficient)

% The analysis for the second 10 years is concluded
% Storing the Optimal solution found
DV_STORE.OPT2 = MISSION2_OPT.MostEfficient;


%% -------- 2050-2060 -------------

% ---------- MISSION DESIGN - 2050-2060 ------------

fprintf('\n----INTERPLANETARY TRANSFER MISSION DESIGN 2050-2060----\n')

% Defining the struct MISSION3 containing all the data 
MISSION3.LEG1.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION3.LEG1.Data.dep_planet = 4; % uplanet ID Mars
MISSION3.LEG1.Data.target_planet = 3; % uplanet ID Earth
MISSION3.LEG2.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION3.LEG2.Data.dep_planet = 3; % uplanet ID Earth
MISSION3.LEG2.Data.target_asteroid = 470864; % ephAsteroid ID
MISSION3.FLYBY.Data.flyby_planet =  3; % uplanet ID Earth
MISSION3.FLYBY.Data.R_planet = 6378; % Equatorial radius of Earth, [km]
MISSION3.FLYBY.Data.mu_planet = astroConstants(13); % [km^3/s^2]

% Time window : ONLY FOR MARS DEPARTURE
% Consider now the last time window : 2050-2060
MISSION3.LEG1.TimeWindow.Departure = [datetime(2050,1,1,0,0,0); datetime(2060,1,1,0,0,0)];
MISSION3.LEG1.TimeWindow.Departure_v = [2050,1,1,0,0,0; 2060,1,1,0,0,0];
% The guess for the ranges of the TOFs is the same of the previous cases,
% dictated by the preliminary analysis
MISSION3.LEG1.TimeWindow.TOF1_min = days(seconds(PRE_LEG1.MostEfficient.TOF/2));
MISSION3.LEG1.TimeWindow.TOF1_max = days(seconds(PRE_LEG1.MostEfficient.TOF*2));
% And also: 
MISSION3.LEG2.TimeWindow.TOF2_min = days(seconds(PRE_LEG2.MostEfficient.TOF/2));
MISSION3.LEG2.TimeWindow.TOF2_max = days(seconds(PRE_LEG2.MostEfficient.TOF*2));

% Setting constraints on the v_launcher, altitude of flyby, ...
MISSION3.LEG1.Data.v_launch = 5; %[km/s]
MISSION3.FLYBY.Data.h_min = 400; % Minimum altitude of flyby hyperbola's pericentre[km]  
...

% 3D Grid search for first estimation of possible location of the minima 
% t_dep_step : 30 days (4)
% tof_step : 20 days (5)
MISSION3 = deltaV_mission_3B(MISSION3,4,5)

% Display the minimum found for later considerations
fprintf('\nMinimum Δv_TOT for MISSION3 (coarse 3D Grid Search Algorithm): %.4f km/s\n',...
    MISSION3.MostEfficient.DELTA_V_TOT);
fprintf('Total TOF: \nTOF1: %d days\n',MISSION3.MostEfficient.TOF1);
fprintf('TOF2: %d days\n',...
    MISSION3.MostEfficient.TOF2);

% Struct PORKCHOP3_1 (associated to MISSION3 for the last 10 years, and leg 1
% driven)
% Porkchop LEG1 driven, x axis : t_dep, y axis : TOF1
% Physical meaning : whatever will be the arrival, what are the best
% regions for the departure and the flyby (t_FB = t_dep + TOF1)?
PORKCHOP3_1.Solution.DELTAV_TOT = min(MISSION3.DELTAV_TOT, [], 3); % [N_dep x N_TOF1]
PORKCHOP3_1.Solution.t1_v = MISSION3.LEG1.Solution.t_dep_v;
PORKCHOP3_1.Solution.t2_v = MISSION3.LEG1.Solution.TOF1_v;

porkchop_date_tof(PORKCHOP3_1)
fig = gcf;
set(fig, 'Name', 'GRID SEARCH 2050-2060, LEG1 driven', 'NumberTitle', 'off');

% Porkchop LEG2 driven , x axis : TOF1, y axis : TOF2
% Physical meaning: whatever will be the departure, what are the best
% regions for the flyby and the arrival (t_arr = t_FB + TOF2)?
PORKCHOP3_2.Solution.DELTAV_TOT = squeeze(min(MISSION3.DELTAV_TOT, [], 1)); % [N_TOF1 x N_TOF2]
PORKCHOP3_2.Solution.t1_v = MISSION3.LEG1.Solution.TOF1_v;
PORKCHOP3_2.Solution.t2_v = MISSION3.LEG2.Solution.TOF2_v;

porkchop_tof_tof(PORKCHOP3_2)
fig = gcf;
set(fig, 'Name', 'GRID SEARCH 2050-2060, LEG2 driven', 'NumberTitle', 'off');

% From the first porkchop plot, we can see that this last time window is
% probably the best window in terms of departure date, and this time the
% regions that are worth it for further investigations are 3
% 1) t_dep : [2054,1,1,0,0,0; 2055,1,1,0,0,0]
% 2) t_dep : [2056,1,1,0,0,0; 2057,1,1,0,0,0]
% 3) t_dep : [2058,1,1,0,0,0; 2059,3,1,0,0,0] , ATTENTION TO THE ARRIVAL
% DATE 
% Instead, from the second porkchop, we can get useful information
% regarding the TOFs and see that we have a 'sub-optimal rectangle' that can
% be assumed with limits : 1) TOF1 : [160, 470]
%                          2) TOF2 : [200, 600]


% ---------- MISSION REFINEMENT - 2050-2060 ------------

fprintf('\n----INTERPLANETARY TRANSFER MISSION REFINEMENT 2050-2060----\n')

% In this phase, considering the last 10 years, we can refine
% the solution before using the optimization tool in order to have the best
% possible initial guess according to our discretization and analysis

% Defining the struct MISSION3_REF1 for the refinement 1
% MISSION3_REF1 refines the first interesting good region in terms of
% departure, MISSION3_REF2 refines instead the second one and lastly MISSION3_REF3
% refines the last one.
% The TOFs are fixed for all the refinements
MISSION3_REF1.LEG1.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION3_REF1.LEG1.Data.dep_planet = 4; % uplanet ID Mars
MISSION3_REF1.LEG1.Data.target_planet = 3; % uplanet ID Earth
MISSION3_REF1.LEG2.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION3_REF1.LEG2.Data.dep_planet = 3; % uplanet ID Earth
MISSION3_REF1.LEG2.Data.target_asteroid = 470864; % ephAsteroid ID
MISSION3_REF1.FLYBY.Data.flyby_planet =  3; % uplanet ID Earth
MISSION3_REF1.FLYBY.Data.R_planet = 6378; % Equatorial radius of Earth, [km]
MISSION3_REF1.FLYBY.Data.mu_planet = astroConstants(13); % [km^3/s^2]

% Time window : ONLY MARS DEPARTURE
MISSION3_REF1.LEG1.TimeWindow.Departure = [datetime(2054,1,1,0,0,0); datetime(2055,1,1,0,0,0)];
MISSION3_REF1.LEG1.TimeWindow.Departure_v = [2054,1,1,0,0,0; 2055,1,1,0,0,0];
MISSION3_REF1.LEG1.TimeWindow.TOF1_min = 160;
MISSION3_REF1.LEG1.TimeWindow.TOF1_max = 470;
% And also: 
MISSION3_REF1.LEG2.TimeWindow.TOF2_min = 200;
MISSION3_REF1.LEG2.TimeWindow.TOF2_max = 600;

% Setting constraints on the v_launcher, altitude of flyby, ...
MISSION3_REF1.LEG1.Data.v_launch = 5; %[km/s]
MISSION3_REF1.FLYBY.Data.h_min = 400; % Minimum altitude of flyby hyperbola's pericentre[km]  
...

% 3D Grid search for first estimation of possible location of the minima
% t_dep_step : 7 days (2)
% tof_step :  days (3)
MISSION3_REF1 = deltaV_mission_3B(MISSION3_REF1,2,3)

% Display the minimum found
fprintf('\nMinimum Δv_TOT for MISSION3_REF1 (refinement 1, first good region): %.4f km/s\n',...
    MISSION3_REF1.MostEfficient.DELTA_V_TOT);
fprintf(['Note : this is not such a good minimum. It wont be considered for the optimization.' ...
    ' Total TOF:\nTOF1: %d days\n'],MISSION3_REF1.MostEfficient.TOF1);
fprintf('TOF2: %d days\n',...
    MISSION3_REF1.MostEfficient.TOF2);
fprintf('Furthermore, also the TOF2 is really high.\n')

% In the same way of before we pass to the porkchop plot analysis:
% Struct PORKCHOP3_REF1_1 (associated to MISSION3 for last 10 years,
% refinement 1 for departure date and leg 1 driven)
% Porkchop LEG1 driven, x axis : t_dep, y axis : TOF1
% Physical meaning : whatever will be the arrival, what are the best
% regions for the departure and the flyby (t_FB = t_dep + TOF1)?
PORKCHOP3_REF1_1.Solution.DELTAV_TOT = min(MISSION3_REF1.DELTAV_TOT, [], 3); % [N_dep x N_TOF1]
PORKCHOP3_REF1_1.Solution.t1_v = MISSION3_REF1.LEG1.Solution.t_dep_v;
PORKCHOP3_REF1_1.Solution.t2_v = MISSION3_REF1.LEG1.Solution.TOF1_v;

porkchop_date_tof(PORKCHOP3_REF1_1)
fig = gcf;
set(fig, 'Name', 'REFINED 1ST REGION 2050-2060, LEG1 driven', 'NumberTitle', 'off');

% Porkchop LEG2 driven , x axis : TOF1, y axis : TOF2
% Physical meaning: whatever will be the departure, what are the best
% regions for the flyby and the arrival (t_arr = t_FB + TOF2)?
PORKCHOP3_REF1_2.Solution.DELTAV_TOT = squeeze(min(MISSION3_REF1.DELTAV_TOT, [], 1)); % [N_TOF1 x N_TOF2]
PORKCHOP3_REF1_2.Solution.t1_v = MISSION3_REF1.LEG1.Solution.TOF1_v;
PORKCHOP3_REF1_2.Solution.t2_v = MISSION3_REF1.LEG2.Solution.TOF2_v;

porkchop_tof_tof(PORKCHOP3_REF1_2)
fig = gcf;
set(fig, 'Name', 'REFINED 1ST REGION 2050-2060, LEG2 driven', 'NumberTitle', 'off');


% Defining the struct MISSION3_REF2 for the refinement 2
MISSION3_REF2.LEG1.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION3_REF2.LEG1.Data.dep_planet = 4; % uplanet ID Mars
MISSION3_REF2.LEG1.Data.target_planet = 3; % uplanet ID Earth
MISSION3_REF2.LEG2.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION3_REF2.LEG2.Data.dep_planet = 3; % uplanet ID Earth
MISSION3_REF2.LEG2.Data.target_asteroid = 470864; % ephAsteroid ID
MISSION3_REF2.FLYBY.Data.flyby_planet =  3; % uplanet ID Earth
MISSION3_REF2.FLYBY.Data.R_planet = 6378; % Equatorial radius of Earth, [km]
MISSION3_REF2.FLYBY.Data.mu_planet = astroConstants(13); % [km^3/s^2]

% Time window : ONLY MARS DEPARTURE
MISSION3_REF2.LEG1.TimeWindow.Departure = [datetime(2056,1,1,0,0,0); datetime(2057,1,1,0,0,0)];
MISSION3_REF2.LEG1.TimeWindow.Departure_v = [2056,1,1,0,0,0; 2057,1,1,0,0,0];
MISSION3_REF2.LEG1.TimeWindow.TOF1_min = 160; 
MISSION3_REF2.LEG1.TimeWindow.TOF1_max = 470;
% And also: 
MISSION3_REF2.LEG2.TimeWindow.TOF2_min = 200;
MISSION3_REF2.LEG2.TimeWindow.TOF2_max = 600;

% Setting constraints on the v_launcher, altitude of flyby, ...
MISSION3_REF2.LEG1.Data.v_launch = 5; %[km/s]
MISSION3_REF2.FLYBY.Data.h_min = 400; % Minimum altitude of flyby hyperbola's pericentre[km]  
...

% 3D Grid search for first estimation of possible location of the minima 
% t_dep_step : 7 days (2)
% tof_step : 4 days (3)
MISSION3_REF2 = deltaV_mission_3B(MISSION3_REF2,2,3)

% Display the minimum found
fprintf('\nMinimum Δv_TOT for MISSION3_REF2 (refinement 2, second good region): %.4f km/s\n',...
    MISSION3_REF2.MostEfficient.DELTA_V_TOT);
fprintf(['Note : this local minimum could be interesting and will be considered for the' ...
    ' optimisation.\nAltough, the TOF2 is quite large.' ...
    '\nTotal TOF:\nTOF1: %d days\n'],MISSION3_REF2.MostEfficient.TOF1);
fprintf('TOF2: %d days\n',...
    MISSION3_REF2.MostEfficient.TOF2);


% Struct PORKCHOP3_REF2_1 (associated to MISSION3 for last 10 years,
% refinement 2 for departure date and leg 1 driven)
% Porkchop LEG1 driven, x axis : t_dep, y axis : TOF1
% Physical meaning : whatever will be the arrival, what are the best
% regions for the departure and the flyby (t_FB = t_dep + TOF1)?
PORKCHOP3_REF2_1.Solution.DELTAV_TOT = min(MISSION3_REF2.DELTAV_TOT, [], 3); % [N_dep x N_TOF1]
PORKCHOP3_REF2_1.Solution.t1_v = MISSION3_REF2.LEG1.Solution.t_dep_v;
PORKCHOP3_REF2_1.Solution.t2_v = MISSION3_REF2.LEG1.Solution.TOF1_v;

porkchop_date_tof(PORKCHOP3_REF2_1)
fig = gcf;
set(fig, 'Name', 'REFINED 2ND REGION 2050-2060, LEG1 driven', 'NumberTitle', 'off');

% Porkchop LEG2 driven , x axis : TOF1, y axis : TOF2
% Physical meaning: whatever will be the departure, what are the best
% regions for the flyby and the arrival (t_arr = t_FB + TOF2)?
PORKCHOP3_REF2_2.Solution.DELTAV_TOT = squeeze(min(MISSION3_REF2.DELTAV_TOT, [], 1)); % [N_TOF1 x N_TOF2]
PORKCHOP3_REF2_2.Solution.t1_v = MISSION3_REF2.LEG1.Solution.TOF1_v;
PORKCHOP3_REF2_2.Solution.t2_v = MISSION3_REF2.LEG2.Solution.TOF2_v;

porkchop_tof_tof(PORKCHOP3_REF2_2)
fig = gcf;
set(fig, 'Name', 'REFINED 2ND REGION 2050-2060, LEG2 driven', 'NumberTitle', 'off');


% Defining the struct MISSION3_REF3 for the refinement 3
MISSION3_REF3.LEG1.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION3_REF3.LEG1.Data.dep_planet = 4; % uplanet ID Mars
MISSION3_REF3.LEG1.Data.target_planet = 3; % uplanet ID Earth
MISSION3_REF3.LEG2.Data.mu_sun = astroConstants(4); %[km^3/s^2]
MISSION3_REF3.LEG2.Data.dep_planet = 3; % uplanet ID Earth
MISSION3_REF3.LEG2.Data.target_asteroid = 470864; % ephAsteroid ID
MISSION3_REF3.FLYBY.Data.flyby_planet =  3; % uplanet ID Earth
MISSION3_REF3.FLYBY.Data.R_planet = 6378; % Equatorial radius of Earth, [km]
MISSION3_REF3.FLYBY.Data.mu_planet = astroConstants(13); % [km^3/s^2]

% Time window : ONLY MARS DEPARTURE
MISSION3_REF3.LEG1.TimeWindow.Departure = [datetime(2058,1,1,0,0,0); datetime(2059,3,1,0,0,0)];
MISSION3_REF3.LEG1.TimeWindow.Departure_v = [2058,1,1,0,0,0; 2059,3,1,0,0,0];
MISSION3_REF3.LEG1.TimeWindow.TOF1_min = 160; 
MISSION3_REF3.LEG1.TimeWindow.TOF1_max = 470;
% And also: 
MISSION3_REF3.LEG2.TimeWindow.TOF2_min = 200;
MISSION3_REF3.LEG2.TimeWindow.TOF2_max = 600;

% Setting constraints on the v_launcher, altitude of flyby, ...
MISSION3_REF3.LEG1.Data.v_launch = 5; %[km/s]
MISSION3_REF3.FLYBY.Data.h_min = 400; % Minimum altitude of flyby hyperbola's pericentre[km]  
...

% 3D Grid search for first estimation of possible location of the minima 
% t_dep_step : 7 days (2)
% tof_step : 4 days (3)
MISSION3_REF3 = deltaV_mission_3B(MISSION3_REF3,2,3)

% Display the minimum found
fprintf('\nMinimum Δv_TOT for MISSION3_REF3 (refinement 3, third good region): %.4f km/s\n',...
    MISSION3_REF3.MostEfficient.DELTA_V_TOT);
fprintf(['Note : this local minimum is the best one yet. But it is important to understand\n' ...
    'the arrival date related to this departure.\nFor mission constraints, it cannot be later' ...
    ' than 01/01/2060 00:00:00.'])
fprintf('\nTotal TOF:\nTOF1: %d days\n',MISSION3_REF3.MostEfficient.TOF1);
fprintf('TOF2: %d days\n',...
    MISSION3_REF3.MostEfficient.TOF2);

% Compute the arrival date for this last region (MJD2000)
MISSION3_REF3.MostEfficient.t_arr = MISSION3_REF3.MostEfficient.t_dep + ...
                                    MISSION3_REF3.MostEfficient.TOF1 + ...
                                    MISSION3_REF3.MostEfficient.TOF2;
% Convert to datetime
t_arr_dt = datetime(mjd20002date(MISSION3_REF3.MostEfficient.t_arr), ...
                    'Format','yyyy-MM-dd HH:mm:ss');
fprintf('The arrival date for the last region is: %s\n', string(t_arr_dt));
fprintf('THIS REGION CANNOT BE CONSIDERED DUE TO SPECIFICATIONS.')
fprintf(['\nUnluckily, we have to discard this region for the optimization despite\n' ...
    'it was the most promising one of the whole analysis.\n'])

% ---------- MISSION OPTIMIZATION - 2050-2060 ------------

fprintf('\n----INTERPLANETARY TRANSFER MISSION OPTIMIZATION 2050-2060----\n')

% The optimization algorithm with fmincon is considered taking
% always into account the constraints on v_launcher and on the minimum
% altitude for the powered flyby.

% Defining the struct MISSION3_OPT  (only one case for these last 10
% years, the case of the refiniment 2)
MISSION3_OPT.LEG1.Data = MISSION3_REF2.LEG1.Data;
MISSION3_OPT.LEG2.Data = MISSION3_REF2.LEG2.Data;
MISSION3_OPT.FLYBY.Data = MISSION3_REF2.FLYBY.Data;
MISSION3_OPT.LEG1.TimeWindow = MISSION3_REF2.LEG1.TimeWindow;
MISSION3_OPT.LEG2.TimeWindow = MISSION3_REF2.LEG2.TimeWindow;

% Setting up the initial guess for the optimization tool to start from,
% chosen as x0 = [t_dep, TOF1, TOF2] in MJD2000 format (3 DOFs).
% The initial guess is found in MISSION3_REF2.MostEfficient field
x0 = [MISSION3_REF2.MostEfficient.t_dep, ...
      MISSION3_REF2.MostEfficient.TOF1,...
      MISSION3_REF2.MostEfficient.TOF2];

% Call optimization function
[MISSION3_OPT, ~, ~, ~, ~] = deltaV_mission_fmincon_3B(MISSION3_OPT, x0);

% Display Optimal solution found
disp(MISSION3_OPT.MostEfficient)

% The analysis for the last 10 years is concluded as the entire
% interplanetary transfer (with powered GA) mission design 
% Storing the Optimal solution found
DV_STORE.OPT3 = MISSION3_OPT.MostEfficient;

%% -------- FINAL COMPARISON ----------

fprintf('\n----INTERPLANETARY TRANSFER MISSION FINAL RESULTS----\n')

% At this point, we can finally compare all the results stored in DV_STORE
% that explore all the possibilities (according to our analysis) for the
% interplanetary mission with departure from Mars, powered flyby at the
% Earth and arrival to the asteroid N. 470864, with early departure
% constraints (1/1/2030 , 00:00:00) and latest arrival constraints
% (1/1/2060 00:00:00).
% In our analysis, we also considered a constraint on the launcher velocity
% at the manouvring point to the Lamber's transfer arc to Earth (without
% considering the injection orbit from Mars), and lastly a constraint to
% the minimum altitude of the perigee of the hyperbola at the close
% encounter with Earth.

% The best possible solution is found taking the minimum, IN TERMS OF
% DELTA_V_TOT, of all the three candidates up to now.
% That is, the 'best' has to be intended only from the total cost of the
% mission point of view.
% It is also important to underline, as a final remark, that the 'global'
% minimum found is an optimized solution according to our analysis, our
% algorithm and discretization. Indeed, it could be found another minimum
% (maybe slightly different) using other type of algorithms and
% discretizations.

% Find the 'global' minima
DELTA_V_TOT_MIN = min([MISSION3_OPT.MostEfficient.DELTA_V_TOT, ...
                       MISSION2_OPT.MostEfficient.DELTA_V_TOT,...
                       MISSION1_OPT.MostEfficient.DELTA_V_TOT]);
fprintf('\nThe global minimum found corresponds to a total cost of the mission of: %.4f',...
     DELTA_V_TOT_MIN);
fprintf('\nThis minimum corresponds to the first good region refined in the 2040-2050 years of MISSION2_OPT\n')

% Define the final struct RESULTS
RESULTS.Solution = MISSION2_OPT.MostEfficient;
disp(RESULTS.Solution)


%% ---------- INTERPLANETARY TRANSFER MISSION : RESULTS -------------- 

% Completing the struct RESULTS containing all the data 
RESULTS.LEG1.Data.mu_sun = astroConstants(4); %[km^3/s^2]
RESULTS.LEG1.Data.dep_planet = 4; % uplanet ID Mars
RESULTS.LEG1.Data.target_planet = 3; % uplanet ID Earth
RESULTS.LEG2.Data.mu_sun = astroConstants(4); %[km^3/s^2]
RESULTS.LEG2.Data.dep_planet = 3; % uplanet ID Earth
RESULTS.LEG2.Data.target_asteroid = 470864; % ephAsteroid ID
RESULTS.FLYBY.Data.flyby_planet =  3; % uplanet ID Earth
RESULTS.FLYBY.Data.R_planet = 6378; % Equatorial radius of Earth, [km]
RESULTS.FLYBY.Data.mu_planet = astroConstants(13); % [km^3/s^2]
RESULTS.FLYBY.Data.h_min =  400; % Minium altitude of flyby [km]

% Creating time vectors
% LEG1 
RESULTS.LEG1.TimeWindow.Departure_v = mjd20002date(RESULTS.Solution.t_dep); 
RESULTS.LEG1.TimeWindow.Departure = datetime(RESULTS.LEG1.TimeWindow.Departure_v);
RESULTS.LEG1.TimeWindow.t_FB_v = mjd20002date( (RESULTS.Solution.t_dep +...
                                                    RESULTS.Solution.TOF1));
RESULTS.LEG1.TimeWindow.t_FB = datetime(RESULTS.LEG1.TimeWindow.t_FB_v);
% Matching flyby
RESULTS.FLYBY.TimeWindow.t_FB_v = RESULTS.LEG1.TimeWindow.t_FB_v;
RESULTS.FLYBY.TimeWindow.t_FB = RESULTS.LEG1.TimeWindow.t_FB;
% LEG2
RESULTS.LEG2.TimeWindow.t_FB_v = RESULTS.FLYBY.TimeWindow.t_FB_v;
RESULTS.LEG2.TimeWindow.t_FB = RESULTS.FLYBY.TimeWindow.t_FB;
RESULTS.LEG2.TimeWindow.Arrival_v = mjd20002date( (RESULTS.Solution.t_dep +...
                                                    RESULTS.Solution.TOF1 +...
                                                    RESULTS.Solution.TOF2));
RESULTS.LEG2.TimeWindow.Arrival = datetime(RESULTS.LEG2.TimeWindow.Arrival_v);


% Reconstructing the entire mission
% --------- LEG1 ------------
[kep_M,~] = uplanet(RESULTS.Solution.t_dep, RESULTS.LEG1.Data.dep_planet);
[kep_E_FB,~] = uplanet(RESULTS.Solution.t_dep + RESULTS.Solution.TOF1, ...
                       RESULTS.FLYBY.Data.flyby_planet);

kep_M(3:6)    = rad2deg(kep_M(3:6));
kep_E_FB(3:6) = rad2deg(kep_E_FB(3:6));

[sM,~] = kep2car(kep_M(1),kep_M(2),kep_M(3), ...
                 kep_M(4),kep_M(5),kep_M(6), RESULTS.LEG1.Data.mu_sun);

[sE_FB,~] = kep2car(kep_E_FB(1),kep_E_FB(2),kep_E_FB(3), ...
                    kep_E_FB(4),kep_E_FB(5),kep_E_FB(6), RESULTS.LEG1.Data.mu_sun);

RESULTS.LEG1.State.dep_planet.kep_M = kep_M;
RESULTS.LEG1.State.dep_planet.car_M = sM;

RESULTS.LEG1.State.flyby_planet.kep_E_FB = kep_E_FB;
RESULTS.LEG1.State.flyby_planet.car_E_FB = sE_FB;

% Lambert arc
[~,~,~,~,vM_t,vE_t1,~,~] = lambertMR( ...
    sM(1:3), sE_FB(1:3), seconds(days(RESULTS.Solution.TOF1)), RESULTS.LEG1.Data.mu_sun, 0);

vM_t  = vM_t';
vE_t1 = vE_t1';

RESULTS.LEG1.Transfer.v_Mt     = vM_t;
RESULTS.LEG1.Transfer.vE_t1     = vE_t1;
RESULTS.LEG1.Transfer.DELTA_V1   = RESULTS.Solution.DELTA_V1;
RESULTS.LEG1.Transfer.TOF1       = RESULTS.Solution.TOF1;

% --------- LEG2 ------------
[kep_E,~] = uplanet(RESULTS.Solution.t_dep + RESULTS.Solution.TOF1, RESULTS.LEG2.Data.dep_planet);
[kep_A,~] = ephAsteroids(RESULTS.Solution.t_dep + RESULTS.Solution.TOF1 + RESULTS.Solution.TOF2, ...
                       RESULTS.LEG2.Data.target_asteroid);

kep_E(3:6)    = rad2deg(kep_E(3:6));
kep_A(3:6) = rad2deg(kep_A(3:6));

[sE,~] = kep2car(kep_E(1),kep_E(2),kep_E(3), ...
                 kep_E(4),kep_E(5),kep_E(6), RESULTS.LEG2.Data.mu_sun);

[sA,~] = kep2car(kep_A(1),kep_A(2),kep_A(3), ...
                    kep_A(4),kep_A(5),kep_A(6), RESULTS.LEG2.Data.mu_sun);

RESULTS.LEG2.State.dep_planet.kep_E = kep_E;
RESULTS.LEG2.State.dep_planet.car_E = sE;

RESULTS.LEG2.State.target_asteroid.kep_A = kep_A;
RESULTS.LEG2.State.target_asteroid.car_A = sA;

% Lambert arc
[~,~,~,~,vE_t2,vA_t,~,~] = lambertMR( ...
    sE(1:3), sA(1:3), seconds(days(RESULTS.Solution.TOF2)), RESULTS.LEG2.Data.mu_sun, 0);

vE_t2  = vE_t2';
vA_t = vA_t';

RESULTS.LEG2.Transfer.vE_t2     = vE_t2;
RESULTS.LEG2.Transfer.vA_t     = vA_t;
RESULTS.LEG2.Transfer.DELTA_V2   = RESULTS.Solution.DELTA_V2;
RESULTS.LEG2.Transfer.TOF2       = RESULTS.Solution.TOF2;

% --------- FLYBY -------------
vE = sE_FB(4:6);
% Incoming hyperbolic excess speed
vinf_min  = vE_t1 - vE;
% Outgoing hyperbolic excess speed
vinf_plus = vE_t2 - vE;   

RESULTS.FLYBY.Solution = flyby_powered(vinf_min, vinf_plus, ...
                              RESULTS.FLYBY.Data.mu_planet, ...
                              RESULTS.FLYBY.Data.R_planet, ...
                              RESULTS.FLYBY.Data.h_min );
% Checking correspondance 
fprintf('\nMission Design DELTAV_FB: %.8f\n', RESULTS.Solution.DELTA_VFB);
fprintf('\nReconstructed DELTAV_FB : %.8f\n',RESULTS.FLYBY.Solution.Deltavp);


% TIME OF PERMANENCE IN EARTH'S SOI 
% (for a powered trailing side fly-by)
RESULTS.FLYBY.Data.r_E = astroConstants(2); % Earth's mean orbital radius around Sun = 1 astronomical unit [Km]
RESULTS.FLYBY.Data.r_SOI = RESULTS.FLYBY.Data.r_E * (RESULTS.FLYBY.Data.mu_planet / RESULTS.LEG2.Data.mu_sun)^(2/5); % radius sphere of influence of Earth wrt Sun 


% 1)
% Conic equation in order to get the true anomaly at the exit of the SOI 
RESULTS.FLYBY.Solution.HYP2.p_hyp = norm(RESULTS.FLYBY.Solution.HYP2.a_hyp)*...
                                          (RESULTS.FLYBY.Solution.HYP2.e_hyp^2-1); 
RESULTS.FLYBY.Solution.HYP2.theta_exitSOI = acos((1/RESULTS.FLYBY.Solution.HYP2.e_hyp)...
    *(RESULTS.FLYBY.Solution.HYP2.p_hyp/RESULTS.FLYBY.Data.r_SOI - 1)); %select the positive result at exit 

% HYP2.e_hyp = eccentricity of the exit hyperbola of the fly-by 
% HYP2.a_hyp = semi-major axis of the exit hyperbola of the fly-by

% Hyperbolic anomaly at the exit of the SOI
RESULTS.FLYBY.Solution.HYP2.F_exitSOI = 2*atanh(sqrt((RESULTS.FLYBY.Solution.HYP2.e_hyp-1)...
    /(RESULTS.FLYBY.Solution.HYP2.e_hyp+1)) * tan(RESULTS.FLYBY.Solution.HYP2.theta_exitSOI/2)); 

% Hyperbolic time law at the exit of the SOI 
RESULTS.FLYBY.Solution.HYP2.n = sqrt(RESULTS.FLYBY.Data.mu_planet/...
                (-RESULTS.FLYBY.Solution.HYP2.a_hyp^3)); % Angular velocity [Hz]
RESULTS.FLYBY.Solution.HYP2.deltat_p_to_exit = 1/RESULTS.FLYBY.Solution.HYP2.n... 
    * (RESULTS.FLYBY.Solution.HYP2.e_hyp * sinh(RESULTS.FLYBY.Solution.HYP2.F_exitSOI) ...
    - RESULTS.FLYBY.Solution.HYP2.F_exitSOI); % Delta t between perigee and exit of SOI

% 2)
% Conic equation in order to get the true anomaly at the enter of the SOI 
RESULTS.FLYBY.Solution.HYP1.p_hyp = norm(RESULTS.FLYBY.Solution.HYP1.a_hyp)...
                                    *(RESULTS.FLYBY.Solution.HYP1.e_hyp^2-1); 
RESULTS.FLYBY.Solution.HYP1.theta_entrySOI = -acos((1/RESULTS.FLYBY.Solution.HYP1.e_hyp)...
                                *(RESULTS.FLYBY.Solution.HYP1.p_hyp/...
                            RESULTS.FLYBY.Data.r_SOI - 1)); %select the negative result at entry 
% HYP1.e_hyp = eccentricity of the entry hyperbola of the fly-by 
% HYP1.a_hyp = semi-major axis of the entry hyperbola of the fly-by

% Hyperbolic anomaly at the entry of the SOI
RESULTS.FLYBY.Solution.HYP1.F_entrySOI = 2*atanh(sqrt((RESULTS.FLYBY.Solution.HYP1.e_hyp-1)...
                                        /(RESULTS.FLYBY.Solution.HYP1.e_hyp+1)) * ...
                                    tan(RESULTS.FLYBY.Solution.HYP1.theta_entrySOI/2)); 

% Hyperbolic time law at the entry of the SOI 
RESULTS.FLYBY.Solution.HYP1.n = sqrt(RESULTS.FLYBY.Data.mu_planet/...
        (-RESULTS.FLYBY.Solution.HYP1.a_hyp^3)); % Angular velocity [Hz]
RESULTS.FLYBY.Solution.HYP1.deltat_entry_to_p = 1/RESULTS.FLYBY.Solution.HYP1.n...
    * (-RESULTS.FLYBY.Solution.HYP1.e_hyp *...
    sinh(RESULTS.FLYBY.Solution.HYP1.F_entrySOI)...
    + RESULTS.FLYBY.Solution.HYP1.F_entrySOI); % Delta t between entry of SOI and perigee

% 3)
% Total delta_t (since the flyby is powered, it's no more symmetric and the two contributions need to be summed up)
RESULTS.FLYBY.Solution.delta_t_tot = RESULTS.FLYBY.Solution.HYP2.deltat_p_to_exit ...
                        + RESULTS.FLYBY.Solution.HYP1.deltat_entry_to_p;


%% ---------- INTERPLANETARY TRANSFER MISSION : PLOTS -------------- 

AU = astroConstants(2); % [km]

% Loading Textures
% Earth Texture 
Re = astroConstants(23)/ AU;
earth_img = imread('EarthTexture.jpg');
earth_img = flipud(earth_img);
[Nx_e, Ny_e] = deal(100, 200); 
[xe, ye, ze] = sphere(Nx_e);
scale_E = 1000;
Xe = xe * Re * scale_E;
Ye = ye * Re * scale_E;
Ze = ze * Re * scale_E;

% Mars Texture
Rm = astroConstants(24)/AU;
mars_img = imread('MarsTexture.jpg');
mars_img = flipud(mars_img);
[Nx_m, Ny_m] = deal(100, 200); 
[xm, ym, zm] = sphere(Nx_m);
scale_M = 2000;
Xm = xm * Rm * scale_M;
Ym = ym * Rm * scale_M;
Zm = zm * Rm * scale_M;

% Sun Texture
Rs = astroConstants(3)/AU;
sun_img = imread('SunTexture.jpg');
sun_img = flipud(sun_img);
[Nx_s, Ny_s] = deal(100, 200); 
[xs, ys, zs] = sphere(Nx_s);
scale_S = 20;
Xs = xs * Rs * scale_S;
Ys = ys * Rs * scale_S;
Zs = zs * Rs * scale_S;


% Key epochs (MJD2000, days)
t_dep = RESULTS.Solution.t_dep;
t_fb  = RESULTS.Solution.t_dep + RESULTS.Solution.TOF1;
t_arr = RESULTS.Solution.t_dep + RESULTS.Solution.TOF1 + RESULTS.Solution.TOF2;

% --- Spacecraft heliocentric Lambert arcs (propagated) ---
% LEG1 (Mars -> Earth flyby)
t_lamb1 = linspace(0, seconds(days(RESULTS.Solution.TOF1)), 1500); % [s]
s0_leg1 = [sM(1:3); vM_t]; % [km; km/s]
st1 = twobody(s0_leg1, RESULTS.LEG1.Data.mu_sun, t_lamb1);

% LEG2 (Earth flyby -> Asteroid)
t_lamb2 = linspace(0, seconds(days(RESULTS.Solution.TOF2)), 1500); % [s]
s0_leg2 = [sE_FB(1:3); vE_t2]; % [km; km/s]
st2 = twobody(s0_leg2, RESULTS.LEG2.Data.mu_sun, t_lamb2);

% --- Heliocentric orbits of the three bodies (sampled over the mission time span) ---
Norb = 2500;
t_orb = linspace(t_dep, t_arr, Norb); % [MJD2000 days]

rM_orb = zeros(Norb,3);
rE_orb = zeros(Norb,3);
rA_orb = zeros(Norb,3);

for k = 1:Norb
    % Mars
    kepM = uplanet(t_orb(k), RESULTS.LEG1.Data.dep_planet);
    kepM(3:6) = rad2deg(kepM(3:6));
    sTmp = kep2car(kepM(1),kepM(2),kepM(3),kepM(4),kepM(5),kepM(6),RESULTS.LEG1.Data.mu_sun);
    rM_orb(k,:) = sTmp(1:3);

    % Earth
    kepE = uplanet(t_orb(k), RESULTS.FLYBY.Data.flyby_planet);
    kepE(3:6) = rad2deg(kepE(3:6));
    sTmp = kep2car(kepE(1),kepE(2),kepE(3),kepE(4),kepE(5),kepE(6),RESULTS.LEG1.Data.mu_sun);
    rE_orb(k,:) = sTmp(1:3);

    % Asteroid
    kepA = ephAsteroids(t_orb(k), RESULTS.LEG2.Data.target_asteroid);
    kepA(3:6) = rad2deg(kepA(3:6));
    sTmp = kep2car(kepA(1),kepA(2),kepA(3),kepA(4),kepA(5),kepA(6),RESULTS.LEG1.Data.mu_sun);
    rA_orb(k,:) = sTmp(1:3);
end

%%--- Positions of ALL the three bodies at departure, flyby and arrival epochs ---
% Departure epoch
kepM_dep = uplanet(t_dep, RESULTS.LEG1.Data.dep_planet);
kepE_dep = uplanet(t_dep, RESULTS.FLYBY.Data.flyby_planet);
kepA_dep = ephAsteroids(t_dep, RESULTS.LEG2.Data.target_asteroid);

kepM_dep(3:6) = rad2deg(kepM_dep(3:6));
kepE_dep(3:6) = rad2deg(kepE_dep(3:6));
kepA_dep(3:6) = rad2deg(kepA_dep(3:6));

sM_dep = kep2car(kepM_dep(1),kepM_dep(2),kepM_dep(3),kepM_dep(4),kepM_dep(5),kepM_dep(6),RESULTS.LEG1.Data.mu_sun);
sE_dep = kep2car(kepE_dep(1),kepE_dep(2),kepE_dep(3),kepE_dep(4),kepE_dep(5),kepE_dep(6),RESULTS.LEG1.Data.mu_sun);
sA_dep = kep2car(kepA_dep(1),kepA_dep(2),kepA_dep(3),kepA_dep(4),kepA_dep(5),kepA_dep(6),RESULTS.LEG1.Data.mu_sun);

% Flyby epoch
kepM_fb = uplanet(t_fb, RESULTS.LEG1.Data.dep_planet);
kepE_fb = uplanet(t_fb, RESULTS.FLYBY.Data.flyby_planet);
kepA_fb = ephAsteroids(t_fb, RESULTS.LEG2.Data.target_asteroid);

kepM_fb(3:6) = rad2deg(kepM_fb(3:6));
kepE_fb(3:6) = rad2deg(kepE_fb(3:6));
kepA_fb(3:6) = rad2deg(kepA_fb(3:6));

sM_fb = kep2car(kepM_fb(1),kepM_fb(2),kepM_fb(3),kepM_fb(4),kepM_fb(5),kepM_fb(6),RESULTS.LEG1.Data.mu_sun);
sE_fb = kep2car(kepE_fb(1),kepE_fb(2),kepE_fb(3),kepE_fb(4),kepE_fb(5),kepE_fb(6),RESULTS.LEG1.Data.mu_sun);
sA_fb = kep2car(kepA_fb(1),kepA_fb(2),kepA_fb(3),kepA_fb(4),kepA_fb(5),kepA_fb(6),RESULTS.LEG1.Data.mu_sun);

% Arrival epoch
kepM_arr = uplanet(t_arr, RESULTS.LEG1.Data.dep_planet);
kepE_arr = uplanet(t_arr, RESULTS.FLYBY.Data.flyby_planet);
kepA_arr = ephAsteroids(t_arr, RESULTS.LEG2.Data.target_asteroid);

kepM_arr(3:6) = rad2deg(kepM_arr(3:6));
kepE_arr(3:6) = rad2deg(kepE_arr(3:6));
kepA_arr(3:6) = rad2deg(kepA_arr(3:6));

sM_arr = kep2car(kepM_arr(1),kepM_arr(2),kepM_arr(3),kepM_arr(4),kepM_arr(5),kepM_arr(6),RESULTS.LEG1.Data.mu_sun);
sE_arr = kep2car(kepE_arr(1),kepE_arr(2),kepE_arr(3),kepE_arr(4),kepE_arr(5),kepE_arr(6),RESULTS.LEG1.Data.mu_sun);
sA_arr = kep2car(kepA_arr(1),kepA_arr(2),kepA_arr(3),kepA_arr(4),kepA_arr(5),kepA_arr(6),RESULTS.LEG1.Data.mu_sun);

% --- Plot: heliocentric trajectory + orbits + body positions at key epochs ---
figure('Name','HELIOCENTRIC - INTERPLANETARY MISSION WITH POWERED GA','NumberTitle','off')
hold on

% Orbits 
hM = plot3(rM_orb(:,1)/AU, rM_orb(:,2)/AU, rM_orb(:,3)/AU, '--', ...
    'LineWidth', 2, 'DisplayName', 'Mars orbit');
hE = plot3(rE_orb(:,1)/AU, rE_orb(:,2)/AU, rE_orb(:,3)/AU, '--', ...
    'LineWidth', 2, 'DisplayName', 'Earth orbit');
hA = plot3(rA_orb(:,1)/AU, rA_orb(:,2)/AU, rA_orb(:,3)/AU, '--', ...
    'LineWidth', 2, 'DisplayName', 'Asteroid orbit');

% Lambert transfer arcs 
plot3(st1(:,1)/AU, st1(:,2)/AU, st1(:,3)/AU, '-', ...
    'LineWidth', 3, 'DisplayName', 'Lambert arc LEG1');
plot3(st2(:,1)/AU, st2(:,2)/AU, st2(:,3)/AU, '-', ...
    'LineWidth', 3, 'DisplayName', 'Lambert arc LEG2');


% --- Body positions at DEPARTURE, FLYBY, ARRIVAL ---

% Sun 
surf(Xs, Ys, Zs, ...
    'CData', sun_img, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'HandleVisibility','off');

% Mars
surf(Xm + sM_dep(1)/AU,Ym + sM_dep(2)/AU,Zm + sM_dep(3)/AU, ...
        'CData', mars_img, 'FaceColor', 'texturemap', 'EdgeColor', 'none', ...
        'HandleVisibility','off');
surf(Xm + sM_fb(1)/AU,Ym + sM_fb(2)/AU,Zm + sM_fb(3)/AU, ...
        'CData', mars_img, 'FaceColor', 'texturemap', 'EdgeColor', 'none', ...
        'HandleVisibility','off');
surf(Xm + sM_arr(1)/AU,Ym + sM_arr(2)/AU,Zm + sM_arr(3)/AU, ...
        'CData', mars_img, 'FaceColor', 'texturemap', 'EdgeColor', 'none', ...
        'HandleVisibility','off');

% Earth
surf(Xe + sE_dep(1)/AU,Ye + sE_dep(2)/AU,Ze + sE_dep(3)/AU, ...
        'CData', earth_img, 'FaceColor', 'texturemap', 'EdgeColor', 'none', ...
        'HandleVisibility','off');
surf(Xe + sE_fb(1)/AU,Ye + sE_fb(2)/AU,Ze + sE_fb(3)/AU, ...
        'CData', earth_img, 'FaceColor', 'texturemap', 'EdgeColor', 'none', ...
        'HandleVisibility','off');
surf(Xe + sE_arr(1)/AU,Ye + sE_arr(2)/AU,Ze + sE_arr(3)/AU, ...
        'CData', earth_img, 'FaceColor', 'texturemap', 'EdgeColor', 'none', ...
        'HandleVisibility','off');

% Asteroid
plot3(sA_dep(1)/AU, sA_dep(2)/AU, sA_dep(3)/AU, 'o', ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.6 0.6 0.6], 'Color', [0.6 0.6 0.6], 'HandleVisibility','off');
plot3(sA_fb(1)/AU,  sA_fb(2)/AU,  sA_fb(3)/AU,  'o', ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.6 0.6 0.6], 'Color', [0.6 0.6 0.6], 'HandleVisibility','off');
plot3(sA_arr(1)/AU, sA_arr(2)/AU, sA_arr(3)/AU, 'o', ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.6 0.6 0.6], 'Color', [0.6 0.6 0.6], 'HandleVisibility','off');

% Optional text labels
text(sM_dep(1)/AU, sM_dep(2)/AU, sM_dep(3)/AU, '  M_{dep}', 'Color', hM.Color);
text(sM_fb(1)/AU,  sM_fb(2)/AU,  sM_fb(3)/AU,  '  M_{fb}',  'Color', hM.Color);
text(sM_arr(1)/AU, sM_arr(2)/AU, sM_arr(3)/AU, '  M_{arr}', 'Color', hM.Color);

text(sE_dep(1)/AU, sE_dep(2)/AU, sE_dep(3)/AU, '  E_{dep}', 'Color', hE.Color);
text(sE_fb(1)/AU,  sE_fb(2)/AU,  sE_fb(3)/AU,  '  E_{fb}',  'Color', hE.Color);
text(sE_arr(1)/AU, sE_arr(2)/AU, sE_arr(3)/AU, '  E_{arr}', 'Color', hE.Color);

text(sA_dep(1)/AU, sA_dep(2)/AU, sA_dep(3)/AU, '  A_{dep}', 'Color', hA.Color);
text(sA_fb(1)/AU,  sA_fb(2)/AU,  sA_fb(3)/AU,  '  A_{fb}',  'Color', hA.Color);
text(sA_arr(1)/AU, sA_arr(2)/AU, sA_arr(3)/AU, '  A_{arr}', 'Color', hA.Color);

grid on
axis equal
xlabel('x [AU]');
ylabel('y [AU]'); 
zlabel('z [AU]');
title('INTERPLANETARY MISSION (heliocentric): Mars – Earth Flyby – Asteroid');
legend('Location','best');

% --- Plot: Earth flyby (planetocentric) - incoming/outgoing hyperbola arcs ---
mu_E = RESULTS.FLYBY.Data.mu_planet;  % [km^3/s^2]
R_E  = RESULTS.FLYBY.Data.R_planet;   % [km]

t_fb_plot = linspace(0, seconds(hours(3)), 1000); % [s] (for visualization)

s0_flyby_in  = [RESULTS.FLYBY.Solution.rp_v; RESULTS.FLYBY.Solution.HYP1.vp_hyp_v];
s0_flyby_out = [RESULTS.FLYBY.Solution.rp_v; RESULTS.FLYBY.Solution.HYP2.vp_hyp_v];

s_flyby_in  = twobody(s0_flyby_in,  mu_E, -t_fb_plot);
s_flyby_out = twobody(s0_flyby_out, mu_E,  t_fb_plot);

figure('Name','PLANETOCENTRIC - EARTH POWERED GA','NumberTitle','off')
hold on
plot3(s_flyby_in(:,1)/R_E,  s_flyby_in(:,2)/R_E,  s_flyby_in(:,3)/R_E, ...
    'LineWidth', 3, 'DisplayName', 'Incoming hyperbola');
plot3(s_flyby_out(:,1)/R_E, s_flyby_out(:,2)/R_E, s_flyby_out(:,3)/R_E, ...
    'LineWidth', 3, 'DisplayName', 'Outgoing hyperbola');

% Pericentre point
plot3(RESULTS.FLYBY.Solution.rp_v(1)/R_E, RESULTS.FLYBY.Solution.rp_v(2)/R_E, RESULTS.FLYBY.Solution.rp_v(3)/R_E, ...
    'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'DisplayName','Pericentre');

% Plot of the asympthotes and apse line
% Calculate the asymptotes and apse line for the hyperbolic trajectory
e_hat = RESULTS.FLYBY.Solution.rp_v / norm(RESULTS.FLYBY.Solution.rp_v);
apse_line = e_hat * linspace(-10*R_E,10*R_E,1000);
hold on
plot3(apse_line(1,:)/R_E, apse_line(2,:)/R_E, apse_line(3,:)/R_E,'k--',...
    'HandleVisibility', 'off','LineWidth',2);
% Center of asymptotes of hyperbolas
C1 = (RESULTS.FLYBY.Solution.rp - RESULTS.FLYBY.Solution.HYP1.a_hyp)*e_hat;
C2 = (RESULTS.FLYBY.Solution.rp - RESULTS.FLYBY.Solution.HYP2.a_hyp)*e_hat;
% Computing Delta Vector 
t_hat_in = vinf_min /  norm(vinf_min);
t_hat_out = vinf_plus / norm(vinf_plus);
% - because u_hat point in the plane
n_hat_in = cross(-RESULTS.FLYBY.Solution.u_hat,t_hat_in);
n_hat_out = cross(-RESULTS.FLYBY.Solution.u_hat,t_hat_out);
Deltav_in = RESULTS.FLYBY.Solution.HYP1.Delta1 * n_hat_in;
Deltav_out = RESULTS.FLYBY.Solution.HYP2.Delta2 * n_hat_out;
asymptote_in = (Deltav_in - C1) * linspace(0,4,1000);
asymptote_out = (Deltav_out - C2) * linspace(0,8,1000);
hold on 
plot3( (C1(1)+asymptote_in(1,:))/R_E,(C1(2)+asymptote_in(2,:))/R_E,(C1(3)+asymptote_in(3,:))/R_E,...
        'LineWidth',2,'DisplayName','Incoming asymptote')
hold on
plot3( (C2(1)+asymptote_out(1,:))/R_E,(C2(2)+asymptote_out(2,:))/R_E,(C2(3)+asymptote_out(3,:))/R_E,...
        'LineWidth',2,'DisplayName','Outcoming asymptote')

% Earth's velocity at flyby
hold on
quiver3(0,0,0,sE_fb(4)/10,sE_fb(5)/10,sE_fb(6)/10,'Color','k','LineWidth',3,'DisplayName','Earth velocity vector');

% Earth (unit radius, plotted in R_E)
[xe, ye, ze] = sphere(100);
earth_img = imread('EarthTexture.jpg');
earth_img = flipud(earth_img);
surf(xe, ye, ze, ...
'CData', earth_img, ...
'FaceColor', 'texturemap', ...
'EdgeColor', 'none', ...
'HandleVisibility', 'off');

grid on
axis equal
xlabel('x [R_E]');
ylabel('y [R_E]');
zlabel('z [R_E]');
title('EARTH FLYBY (planetocentric): incoming/outgoing hyperbola arcs');
legend('Location','best');
