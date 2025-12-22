%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%                                                         %%%%%%%%%
%%%%%%%%%                   ORBITAL MECHANICS PROJECT             %%%%%%%%%
%%%%%%%%%                       A.A. 2025-2026                    %%%%%%%%%
%%%%%%%%%                                                         %%%%%%%%%
%%%%%%%%%                          GROUP 2519                     %%%%%%%%%
%%%%%%%%%                       ASSIGNMENT 1-2                    %%%%%%%%%
%%%%%%%%%                                                         %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CONTRIBUTORS
% Amura Fabio, 10830618, fabio.amura@mail.polimi.it
% Crisanti Costanza, 10911209, costanza.crisanti@mail.polimi.it
% Deli Luca 
% Tomas (non so i cognomi)

% REQUIREMENTS : Assignment 1
% Asteroid N.470864 
% Departure date : 00:00:00 01/01/2030 
% Arrival date: 00:00:00 01/01/2060 
% Departure planet : Mars
% Flyby planet : Earth
% 
% 
% REQUIREMENTS : Assignment 2
% Central planet Earth 
%  a = 2.4979e4 [km] 
%  e = 0.5162 [-]
%  i = 59 [deg]
%  w = 43 [deg]
%  Raan = 53 [deg]
%  k = 1
%  m = 5
%  Perturbations: J2, Moon 
%  Initial Date: 2041-6-28-22-40-59

clear
close all
clc

%% --------------- ASSIGNMENT 1 ---------------

% --------- LEG 1 --------
% Defining a struct containing all the data of the LEG1
LEG1.Data.mu_sun = astroConstants(4); %[km^3/s^2]
LEG1.Data.dep_planet = 4; % uplanet ID Mars
LEG1.Data.target_planet = 3; % uplanet ID Earth
% ARBITRARLY CHOOSE THE DATES FOR NOW
LEG1.TimeWindow.Departure = [datetime(2030,1,1,0,0,0); datetime(2060,1,1,0,0,0)];
LEG1.TimeWindow.Departure_v = [2030,1,1,0,0,0; 2060,1,1,0,0,0];
LEG1.TimeWindow.Arrival = [datetime(2030,1,1,0,0,0); datetime(2060,1,1,0,0,0)];
LEG1.TimeWindow.Arrival_v = [2030,1,1,0,0,0; 2060,1,1,0,0,0];
% Eventual constraint on v_launch : UNCONSTRAINED
LEG1.Data.v_launch = inf; %[km/s]

% Computing the DeltaV of the first leg 
time_step_case = 7; % 1 time step each 14 days
[LEG1] = deltaV_mission(LEG1,time_step_case);
% Refine solution using fmincon : initial guess is the one contained in the
% MostEfficient field of LEG1
x0 = [LEG1.MostEfficient.t1, LEG1.MostEfficient.t2];
[LEG1, ~, ~, exitflag, ~] = deltaV_mission_fmincon(LEG1, x0);

% Minimum deltaV transfer for LEG1
disp(LEG1.MostEfficient);

% ------- LEG 2 ----------
% Defining a struct containing all the data of the LEG1
LEG2.Data.mu_sun = astroConstants(4); %[km^3/s^2]
LEG2.Data.dep_planet = 3; % uplanet ID Earth
LEG2.Data.target_asteroid = 470864; % uplanet ID Earth
% ARBITRARLY CHOOSE THE DATES FOR NOW
LEG2.TimeWindow.Departure = [datetime(2030,1,1,0,0,0); datetime(2060,1,1,0,0,0)];
LEG2.TimeWindow.Departure_v = [2030,1,1,0,0,0; 2060,1,1,0,0,0];
LEG2.TimeWindow.Arrival = [datetime(2030,1,1,0,0,0); datetime(2060,1,1,0,0,0)];
LEG2.TimeWindow.Arrival_v = [2030,1,1,0,0,0; 2060,1,1,0,0,0];
% Eventual constraint on v_launch : UNCONSTRAINED
LEG2.Data.v_launch = inf; %[km/s]

% Computing the DeltaV of the first leg 
time_step_case = 7; % 1 time step each 14 days
[LEG2] = deltaV_mission_asteroid(LEG2,time_step_case);
% Refine solution using fmincon : initial guess is the one contained in the
% MostEfficient field of LEG1
x0 = [LEG2.MostEfficient.t1, LEG2.MostEfficient.t2];
[LEG2, ~, ~, exitflag, ~] = deltaV_mission_fmincon_asteroid(LEG2, x0);

% Minimum deltaV transfer for LEG
disp(LEG2.MostEfficient);

% -------- POWERED GRAVITY ASSIST DESIGN -------
% Defining the already known data
FLYBY.Data.V_minus = LEG1.MostEfficient.v_t2_tr; % S/c absolute veolcity wrt Sun at entry of flyby
FLYBY.Data.V_planet = LEG1.MostEfficient.v_t2; % Planet absolute veolcity wrt Sun
FLYBY.Data.V_plus = LEG2.MostEfficient.v_t1_tr; % S/c absolute velocity wrt Sun at exit of flyby
if norm(FLYBY.Data.V_minus)>norm(FLYBY.Data.V_plus)
    FLYBY.Data.type = "leading"; 
else
    FLYBY.Data.type = "trailing"; 
end 
FLYBY.Data.mu_E = astroConstants(13);

% Finding the relative velocities wrt the planet
FLYBY.Velocities.v_inf_minus = FLYBY.Data.V_minus - FLYBY.Data.V_planet; 
FLYBY.Velocities.v_inf_plus = FLYBY.Data.V_plus - FLYBY.Data.V_planet;
FLYBY.Velocities.Delta_v_flyby = FLYBY.Velocities.v_inf_plus - FLYBY.Velocities.v_inf_minus; % Total change of velocity
FLYBY.Velocities.Delta_v_flyby_norm = norm(FLYBY.Velocities.Delta_v_flyby); 

% Total turning angle 
FLYBY.Arcs.delta = acos(dot(FLYBY.Velocities.v_inf_minus, FLYBY.Velocities.v_inf_plus) / (norm(FLYBY.Velocities.v_inf_minus) * norm(FLYBY.Velocities.v_inf_plus))); % Total turinig angle [rad] 

% Solving for radius of pericentre, using equation solver
FLYBY.Arcs.e_0 = 1/sin(FLYBY.Arcs.delta/2);
FLYBY.Arcs.a_0 = -FLYBY.Data.mu_E /(norm(FLYBY.Velocities.v_inf_plus))^2; 
FLYBY.Arcs.rp_0 = FLYBY.Arcs.a_0*(1-FLYBY.Arcs.e_0);% First guess  
fun = @(rp) asin(1/(1+rp*(norm(FLYBY.Velocities.v_inf_minus))^2/FLYBY.Data.mu_E)) +asin(1/(1+rp*(norm(FLYBY.Velocities.v_inf_plus))^2/FLYBY.Data.mu_E)) - FLYBY.Arcs.delta; 
FLYBY.Arcs.rp = fzero(fun, FLYBY.Arcs.rp_0); 
FLYBY.Data.R_E = astroConstants(23); 
FLYBY.Data.h_atm = 100; %[km] Karman Line
%if FLYBY.Arcs.rp > (FLYBY.Data.R_E + FLYBY.Data.h_atm)
%   fprintf("valid rp = %f\n", FLYBY.Arcs.rp); 
%else
%    error("Impact"); 
%end 
% ATTENTION: NOW WE HAVE IMPACT -> ADJUST DELTA_V (THIS IS A CONSTRAINT)

% Arcs characteristics
FLYBY.Arcs.e_plus = 1+FLYBY.Arcs.rp*(norm(FLYBY.Velocities.v_inf_plus))^2/FLYBY.Data.mu_E; 
FLYBY.Arcs.a_plus = -FLYBY.Data.mu_E /(norm(FLYBY.Velocities.v_inf_plus))^2; 
FLYBY.Arcs.delta_plus =2*asin(1/FLYBY.Arcs.e_plus);
FLYBY.Arcs.e_minus = 1+FLYBY.Arcs.rp*(norm(FLYBY.Velocities.v_inf_minus))^2/FLYBY.Data.mu_E;
FLYBY.Arcs.a_minus = -FLYBY.Data.mu_E /(norm(FLYBY.Velocities.v_inf_minus))^2; 
FLYBY.Arcs.delta_minus =2*asin(1/FLYBY.Arcs.e_minus);

% Velocities at pericentre and cost of the gravity assist
FLYBY.Arcs.u = cross(FLYBY.Velocities.v_inf_minus,FLYBY.Velocities.v_inf_plus);   
FLYBY.Arcs.u = FLYBY.Arcs.u / norm(FLYBY.Arcs.u); % normal to plane of arcs (APPROXIMATION - in realty, v_inf_plus and v_inf_minus are NOT on the same plane)
FLYBY.Arcs.direction = Rodrigues(FLYBY.Velocities.v_inf_minus,FLYBY.Arcs.u, FLYBY.Arcs.delta_minus/2); 
FLYBY.Arcs.direction = FLYBY.Arcs.direction/norm(FLYBY.Arcs.direction);
FLYBY.Velocities.vp_minus = sqrt(norm(FLYBY.Velocities.v_inf_minus)^2+2*FLYBY.Data.mu_E/FLYBY.Arcs.rp)*FLYBY.Arcs.direction; 
FLYBY.Velocities.vp_plus = sqrt(norm(FLYBY.Velocities.v_inf_plus)^2+2*FLYBY.Data.mu_E/FLYBY.Arcs.rp)*FLYBY.Arcs.direction; 
FLYBY.Velocities.Delta_v_ga = FLYBY.Velocities.vp_plus - FLYBY.Velocities.vp_minus; % Delta V = manoeuvre to be given at pericentre
FLYBY.Velocities.Delta_v_ga_norm = norm (FLYBY.Velocities.Delta_v_ga); 
% CHECK THE APPROXIMATION USED FOR VECTORS DIRECTIONS


%% -------------- PLOTS : ASSIGNMENT 1 -----------------

% Contour plot (porkchop plot) of LEG1 
porkchop_plot(LEG1);

% Contour plot (porkchop plot) of LEG2
porkchop_plot(LEG2)

% Propagate the orbits, initial conditions chosen as t1 and t2 of the most
% efficient transfer
% Converting into cartesian state for initial condition state
[s0_1,~] = kep2car( LEG1.MostEfficient.Eph_dep_planet(1),...
                    LEG1.MostEfficient.Eph_dep_planet(2),...
                    rad2deg(LEG1.MostEfficient.Eph_dep_planet(3)),...
                    rad2deg(LEG1.MostEfficient.Eph_dep_planet(4)),...
                    rad2deg(LEG1.MostEfficient.Eph_dep_planet(5)),...
                    rad2deg(LEG1.MostEfficient.Eph_dep_planet(6)),...
                    LEG1.Data.mu_sun);
[s0_2,~] = kep2car( LEG1.MostEfficient.Eph_target_planet(1),...
                    LEG1.MostEfficient.Eph_target_planet(2),...
                    rad2deg(LEG1.MostEfficient.Eph_target_planet(3)),...
                    rad2deg(LEG1.MostEfficient.Eph_target_planet(4)),...
                    rad2deg(LEG1.MostEfficient.Eph_target_planet(5)),...
                    rad2deg(LEG1.MostEfficient.Eph_target_planet(6)),...
                    LEG1.Data.mu_sun);
% Define the period of the orbit 1 (Mars)
T_1 = 2*pi*sqrt(LEG1.MostEfficient.Eph_dep_planet(1)^3/LEG1.Data.mu_sun);      % [s]
% Time vector 1
t_v_1 = linspace(0,T_1,1000);
% Define the period of the orbit 2 (Earth)
T_2 = 2*pi*sqrt(LEG1.MostEfficient.Eph_target_planet(1)^3/LEG1.Data.mu_sun);      % [s]
% Time vector 2
t_v_2 = linspace(0,T_2,1000);

% Solving Lambert's Problem for the transfer arc
[a_t,~,~,ERROR,v_t1,v_t2,~,~] = lambertMR(s0_1(1:3),s0_2(1:3),...
                                    LEG1.MostEfficient.TOF,LEG1.Data.mu_sun,0);
% Define the period of the transfer orbit
T_t = 2*pi*sqrt(a_t^3/LEG1.Data.mu_sun);     % [s]
% Time vector of the transfer orbit
t_v_t = linspace(0,LEG1.MostEfficient.TOF,1000);
% Initial state of the transfer orbit
s0_t = [s0_1(1:3); v_t1'];

% Propagation of the orbits 
[s1] = twobody(s0_1,LEG1.Data.mu_sun,t_v_1);     % Full Orbit 1 
[s1_t] = twobody(s0_1,LEG1.Data.mu_sun,t_v_t);   % Planet 1 during transfer 
% PROPAGATION BACKWARD IN TIME FOR ORBIT 2
[s2] = twobody(s0_2, LEG1.Data.mu_sun, -t_v_2);   % Full Orbit 2
[s2_t] = twobody(s0_2, LEG1.Data.mu_sun, -t_v_t); % Planet 2 during transfer
[st] = twobody(s0_t,LEG1.Data.mu_sun,t_v_t);  

% Orbit plots 
AU = astroConstants(2);

% Earth Texture 
Re = astroConstants(23)/ AU;
earth_img = imread('EarthTexture.jpg');
earth_img = flipud(earth_img);
[Nx_e, Ny_e] = deal(100, 200); 
[xe, ye, ze] = sphere(Nx_e);
scale_E = 1000;
% Initial position of the Earth
Xe_0 = xe * Re * scale_E + s2_t(end,1)/AU; 
Ye_0 = ye * Re * scale_E + s2_t(end,2)/AU;
Ze_0 = ze * Re * scale_E + s2_t(end,3)/AU;
% Final position of the Earth
Xe_f = xe * Re * scale_E + s0_2(1)/AU; 
Ye_f = ye * Re * scale_E + s0_2(2)/AU;
Ze_f = ze * Re * scale_E + s0_2(3)/AU;

% Mars Texture
Rm = astroConstants(24)/AU;
mars_img = imread('MarsTexture.jpg');
mars_img = flipud(mars_img);
[Nx_m, Ny_m] = deal(100, 200); 
[xm, ym, zm] = sphere(Nx_m);
scale_M = 2000;
% Initial position of Mars
Xm_0 = xm * Rm * scale_M + s0_1(1)/AU; 
Ym_0 = ym * Rm * scale_M + s0_1(2)/AU;
Zm_0 = zm * Rm * scale_M + s0_1(3)/AU;
% Final position of  Mars 
Xm_f = xm * Rm * scale_M + s1_t(end,1)/AU; 
Ym_f = ym * Rm * scale_M + s1_t(end,2)/AU;
Zm_f = zm * Rm * scale_M + s1_t(end,3)/AU;

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

figure
% Departure planet full orbit
E1 = plot3(s1(:,1)/AU, s1(:,2)/AU, s1(:,3)/AU,...
        '--', 'LineWidth',3, 'DisplayName','Mars Orbit');
% Departure planet orbit during transfer
hold on
plot3(s1_t(:,1)/AU, s1_t(:,2)/AU, s1_t(:,3)/AU,...
        'Color',E1.Color, 'LineWidth',3,'DisplayName','Mars Motion during Transfer')
clear E1
% Target planet full orbit
hold on
E2 = plot3(s2(:,1)/AU, s2(:,2)/AU, s2(:,3)/AU,...
        '--', 'LineWidth',3,'DisplayName','Earth Orbit');
% Target planet orbit during transfer
hold on
plot3(s2_t(:,1)/AU, s2_t(:,2)/AU, s2_t(:,3)/AU,...
         'Color',E2.Color,'LineWidth',3,'DisplayName','Earth Motion during Transfer')
clear E2

% Transfer Arc
hold on
plot3(st(:,1)/AU, st(:,2)/AU, st(:,3)/AU,...
        'LineWidth',3,'DisplayName','Transfer Arc');
% Earth Initial Position
hold on 
surf(Xe_0, Ye_0, Ze_0, ...
    'CData', earth_img, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'HandleVisibility','off');
% Earth Final Poition
surf(Xe_f, Ye_f, Ze_f, ...
    'CData', earth_img, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'HandleVisibility','off');
% Mars Initial Position
hold on 
surf(Xm_0, Ym_0, Zm_0, ...
    'CData', mars_img, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'HandleVisibility','off');
% Mars Final Poition
surf(Xm_f, Ym_f, Zm_f, ...
    'CData', mars_img, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'HandleVisibility','off');
% Sun 
surf(Xs, Ys, Zs, ...
    'CData', sun_img, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'HandleVisibility','off');
title('MARS-EARTH TRANSFER LEG1')
grid on
xlabel('x [AU]')
ylabel('y [AU]')
zlabel('z [AU]')
legend('Location','northeast')
axis equal

% Projection onto x-y plane
figure
% Departure planet full orbit
E1 = plot3(s1(:,1)/AU, s1(:,2)/AU, s1(:,3)/AU,...
        '--', 'LineWidth',3, 'DisplayName','Mars Orbit');
% Departure planet orbit during transfer
hold on
plot3(s1_t(:,1)/AU, s1_t(:,2)/AU, s1_t(:,3)/AU,...
        'Color',E1.Color, 'LineWidth',3,'DisplayName','Mars Motion during Transfer')
clear E1
% Target planet full orbit
hold on
E2 = plot3(s2(:,1)/AU, s2(:,2)/AU, s2(:,3)/AU,...
        '--', 'LineWidth',3,'DisplayName','Earth Orbit');
% Target planet orbit during transfer
hold on
plot3(s2_t(:,1)/AU, s2_t(:,2)/AU, s2_t(:,3)/AU,...
         'Color',E2.Color,'LineWidth',3,'DisplayName','Earth Motion during Transfer')
clear E2

% Transfer Arc
hold on
plot3(st(:,1)/AU, st(:,2)/AU, st(:,3)/AU,...
        'LineWidth',3,'DisplayName','Transfer Arc');
% Earth Initial Position
hold on 
surf(Xe_0, Ye_0, Ze_0, ...
    'CData', earth_img, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'HandleVisibility','off');
% Earth Final Poition
surf(Xe_f, Ye_f, Ze_f, ...
    'CData', earth_img, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'HandleVisibility','off');
% Mars Initial Position
hold on 
surf(Xm_0, Ym_0, Zm_0, ...
    'CData', mars_img, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'HandleVisibility','off');
% Mars Final Poition
surf(Xm_f, Ym_f, Zm_f, ...
    'CData', mars_img, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'HandleVisibility','off');
% Sun 
surf(Xs, Ys, Zs, ...
    'CData', sun_img, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'HandleVisibility','off');
title('MARS-EARTH TRANSFER LEG1')
grid on
xlabel('x [AU]')
ylabel('y [AU]')
zlabel('z [AU]')
legend('Location','northeast')
axis equal
view(2)