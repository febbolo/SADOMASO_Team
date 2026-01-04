clear
clc
close all

% _________________ Data __________________

% REQUIREMENTS : Assignment 2
% Central planet Earth 
%  a = 2.4979e4 [km] 
%  e = 0.5162 [-]
%  i = 59 [deg]
%  w = 43 [deg]
%  Raan = 53 [deg]
%  f0 = 0 [deg] always, as written in the ass.2
%  k = 1
%  m = 5
%  Perturbations: J2, Moon 
%  Initial Date: 2041-6-28-22-40-59

% declaring costants of the primary using astroconstants, OBLATENESS for
% the report

Primary.mu = astroConstants(13);
Primary.Radius = astroConstants(23);
Primary.mass = Primary.mu / astroConstants(1);
Primary.J2 = astroConstants(9);
Primary.Initial_date_Mjd = date2mjd2000([2041, 6, 28, 22, 40, 59]); % converting in Julian the initial date to exploit it for ephmoon
Primary.w = deg2rad(15.04)/3600;

% _________________ For earth plot and planisphere __________________

earth_img = imread('EarthTexture.jpg');
earth_img = flipud(earth_img);
[Nx,Ny] = deal(100,200);

[xe, ye, ze] = sphere(Nx);
xe = xe * Primary.Radius;
ye = ye * Primary.Radius;
ze = ze * Primary.Radius;

% Loading the planisphere image and defining its longitude/latitude axes
earth = imread('EarthTexture.jpg');

% Angle values contained between -180 and 180 degrees for longitude
lon_img = linspace(-180, 180, size(earth,2));

% Angle values contained between -90 and 90 degrees for latitude
lat_img = linspace(90, -90, size(earth,1));

% _________________ Data from sheet __________________

a     = 2.4979e4; % km
e     = 0.5162;
i     = deg2rad(59);
Raan  = deg2rad(53);
omega = deg2rad(43);
f0    = 0;
kep = [a, e, i, Raan, omega, f0];

% For repeating Ground Track
k = 1; 
m = 5;

% _________________ Ground Track __________________

GST_deg = 280.46061837 + 360.98564736629 * Primary.Initial_date_Mjd;
ThetaG0 = mod(GST_deg * pi/180, 2*pi); 

% __________________ Nominal (non-repeating) orbit __________________

a_nom = a;    % Nominal semi-major axis [km]
n_nom = sqrt(Primary.mu/(a_nom^3));
T_nom = 2*pi/n_nom;

% Nominal Keplerian vector (same order as Kep2Car_vec: [a e i Om om theta])
kep_nom = [a_nom, e, i, Raan, omega, f0];

% Converting Keplerian to Cartesian for propagation
[r0_nom, v0_nom] = Kep2Car_vec(kep_nom, Primary.mu);
y0_nom = [r0_nom; v0_nom];

% __________________ Repeating Ground Track orbit __________________

% Computing the orbital period required for a repeating ground track
T_rgt = T_RGT(k, m);
n_rgt = (2*pi)/T_rgt;
a_rgt = (Primary.mu/(n_rgt^2))^(1/3);    % Semi-major axis for RGT orbit

% Updating Keplerian elements for RGT orbit
kep_rgt    = kep;
kep_rgt(1) = a_rgt;

% Converting
[r0_rgt, v0_rgt] = Kep2Car_vec(kep_rgt, Primary.mu);

% Building initial conditions vector for the propagation
y0_rgt = [r0_rgt; v0_rgt];


% % __________________ SINGLE 2BP propagation for both __________________
% 
% options = odeset('RelTol',1e-9,'AbsTol',1e-10);
% 
% % 2BP propagation for nominal (non-repeating) orbit
% tspan_nom = linspace(0, 3*T_rgt, 50000);
% [t_nom, y_nom] = ode113(@(t,y) fun_2bp(t,y,Primary.mu), tspan_nom, y0_nom, options);
% r_nom = y_nom(:,1:3);
% 
% % 2BP propagation for repeating ground track orbit
% tspan_rgt = linspace(0, 3* T_rgt, 50000);   % <-- must be T_rgt, not T_nom
% [t_rgt, y_rgt] = ode113(@(t,y) fun_2bp(t,y,Primary.mu), tspan_rgt, y0_rgt, options);
% r_rgt = y_rgt(:,1:3);
% 
% %% __________________ 3D orbit plot __________________
% 
% figure
% hold on
% 
% surf(xe, ye, ze,'CData', earth_img,'FaceColor','texturemap', 'EdgeColor','none', 'HandleVisibility','off');
% 
% plot3(r_nom(:,1), r_nom(:,2), r_nom(:,3), 'LineWidth', 1.2, 'Color', 'g', 'DisplayName', 'Nominal orbit');
% plot3(r_rgt(:,1), r_rgt(:,2), r_rgt(:,3),'LineWidth', 1.2, 'Color', 'r', 'DisplayName', 'RGT orbit');
% 
% xlabel('x [km]')
% ylabel('y [km]')
% zlabel('z [km]')
% legend show
% axis equal
% grid on
% view(3);
% 
% % __________________ Ground track computation __________________
% 
% [alpha_nom, delta_nom, lon_nom, lat_nom] = Ground_track(t_nom, r_nom, ThetaG0, t_nom(1), Primary.w);
% [alpha_rgt, delta_rgt, lon_rgt, lat_rgt] = Ground_track(t_rgt, r_rgt, ThetaG0, t_rgt(1), Primary.w);
% 
% % Unwrap longitude to keep continuity, then convert to degrees
% 
% % __________________ Ground track plot __________________
% 
% figure
% hold on
% image(lon_img, lat_img, earth);
% set(gca,'YDir','normal');
% 
% plot(lon_nom, lat_nom, 'g', 'LineStyle','none', 'Marker','.', 'Color','g', 'MarkerSize',2,'DisplayName','Nominal ground track');
% plot(lon_rgt, lat_rgt, 'r', 'LineStyle','none', 'Marker','.', 'Color','r', 'MarkerSize',2, 'DisplayName', 'RGT ground track');
% plot(lon_nom(1),  lat_nom(1),  'go', 'MarkerSize',6, 'DisplayName','Nominal start');
% plot(lon_nom(end),lat_nom(end),'gx', 'MarkerSize',6, 'DisplayName','Nominal end');
% plot(lon_rgt(1),  lat_rgt(1),  'ro', 'MarkerSize',6, 'DisplayName','RGT start');
% plot(lon_rgt(end),lat_rgt(end),'rx', 'MarkerSize',6, 'DisplayName','RGT end');
% 
% xlabel('Longitude [deg]')
% ylabel('Latitude [deg]')
% title('Ground Track')
% legend ("show", 'Location', 'northeastoutside')
% xlim([-180 180]); ylim([-90 90]);
% grid on
% axis equal tight
 
%% ________________________________________________________________________
% This selection is if I want to run the 4 following case or If I want to
% run one of the four directly to see the results, I only need to change
% the index in case I want to run only one, or change the string to all If
% I want to run them all, for tha animation instead, I chose to switch on
% or off the animation during the plot

% _________________ Animation settings __________________
% ANIMATE_3D = true;
% ANIMATE_GT = true;
ans_3d = input('Enable ANIMATE_3D? (true/false) [default: true]: ', 's');
if isempty(ans_3d), ANIMATE_3D = true; else, ANIMATE_3D = strcmpi(ans_3d, 'true'); end

ans_gt = input('Enable ANIMATE_GT? (true/false) [default: true]: ', 's');
if isempty(ans_gt), ANIMATE_GT = true; else, ANIMATE_GT = strcmpi(ans_gt, 'true'); end

STEP_ANIM = 2000;      % bigger = faster animation

% _________________ Case selection __________________
% RUN_MODE = 'all'    -> run all 4 cases
% RUN_MODE = 'single' -> run one case (CASE_ID = 1..4)
% RUN_MODE = 'all';
% CASE_ID  = 1;

RUN_MODE = input('RUN_MODE (single/all) [default: all]: ', 's');
if isempty(RUN_MODE)
    RUN_MODE = 'all'; 
end

CASE_ID = 1; % Default Case ID
if strcmpi(RUN_MODE, 'single')
    sel_id = input('CASE_ID (\n1: T_nom, \n2: 10*T_Earth, \n3: 2T_rgt, \n4: T_Earth, \n5: 5T_Earth) \n[default: 1]: ', 's');
    if ~isempty(sel_id)
        CASE_ID = str2double(sel_id);
    end
end

T_Earth = 2*pi / Primary.w;     % Earth rotation period (from wE in Radians)
case_names = { ...
    '1 nominal orbital period (T_{nom})', ...
    '10 Earth revolutions (T_{rgt})', ...
    '2 RGT orbital periods (2T_{rgt})', ...
    '1 Earth rotation (TEarth)', ...
    '5 Eart revolutions (5TEarth)'};

% Same time span for BOTH orbits, depending on the selected case
case_Tend = [ T_nom, 10*T_Earth, 5*T_rgt, T_Earth, 5*T_Earth ];

if strcmpi(RUN_MODE,'all')
    case_list = 1:5;
else
    case_list = CASE_ID;
end

% Number of points for propagation
Npts = 200000;

% _________________ Loop over selected cases __________________

for icase = case_list

    Tend = case_Tend(icase);   % <-- same Tend (Plotting_period) for nominal and RGT

    % _________________ 2BP propagation for both __________________

    options = odeset('RelTol',1e-9,'AbsTol',1e-10);

    tspan_nom = linspace(0, Tend, Npts);
    [t_nom, y_nom] = ode113(@(t,y) fun_2bp(t,y,Primary.mu), tspan_nom, y0_nom, options);
    r_nom = y_nom(:,1:3);

    tspan_rgt = linspace(0, Tend, Npts);
    [t_rgt, y_rgt] = ode113(@(t,y) fun_2bp(t,y,Primary.mu), tspan_rgt, y0_rgt, options);
    r_rgt = y_rgt(:,1:3);

    % _________________ Ground track computation (INSIDE CASE LOOP) __________________
    % Ground_track already returns lon/lat in degrees and wrapped, so no extra unwrap/mod here.

    [alpha_nom, delta_nom, lon_nom, lat_nom] = Ground_track(t_nom, r_nom, ThetaG0, t_nom(1), Primary.w);
    [alpha_rgt, delta_rgt, lon_rgt, lat_rgt] = Ground_track(t_rgt, r_rgt, ThetaG0, t_rgt(1), Primary.w);

    % _________________ 3D orbit plot __________________

    close all

    fig3D = figure;
    hold on

    surf(xe, ye, ze, 'CData', earth_img, 'FaceColor','texturemap', 'EdgeColor','none', 'HandleVisibility','off');

    if ANIMATE_3D  % I decide each time if I want to plot or not the animation

        % h_nom/h_rgt are the animated line objects (orbit history)
        h_nom = plot3(NaN,NaN,NaN, 'g', 'LineWidth', 1.2, 'DisplayName', 'Nominal orbit');
        h_rgt = plot3(NaN,NaN,NaN, 'r', 'LineWidth', 1.2, 'DisplayName', 'RGT orbit');

        % hDot_nom/hDot_rgt are the moving markers (current spacecraft position)
        hDot_nom = plot3(NaN,NaN,NaN, 'go', 'MarkerSize', 6, 'LineWidth', 1.5, 'HandleVisibility','off');
        hDot_rgt = plot3(NaN,NaN,NaN, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5, 'HandleVisibility','off');

    else
        plot3(r_nom(:,1), r_nom(:,2), r_nom(:,3),'LineWidth', 5, 'Color', 'g', 'DisplayName', 'Nominal orbit')
        % plot3(r_rgt(:,1), r_rgt(:,2), r_rgt(:,3),'LineWidth', 5, 'Color', 'r', 'DisplayName', 'RGT orbit')
    end

    xlabel('x [km]')
    ylabel('y [km]')
    zlabel('z [km]')
    legend show
    axis equal
    grid on
    view(3);
    title(['3D Orbits - ' case_names{icase}])

    if ANIMATE_3D

        % Nanim is the number of animation steps (limited by the shortest history)
        Nanim = min(size(r_nom,1), size(r_rgt,1));

        for kk = 1:STEP_ANIM:Nanim

            set(h_nom, 'XData', r_nom(1:kk,1), 'YData', r_nom(1:kk,2), 'ZData', r_nom(1:kk,3));
            set(h_rgt, 'XData', r_rgt(1:kk,1), 'YData', r_rgt(1:kk,2), 'ZData', r_rgt(1:kk,3));

            % Update moving dots at current index
            set(hDot_nom, 'XData', r_nom(kk,1), 'YData', r_nom(kk,2), 'ZData', r_nom(kk,3));
            set(hDot_rgt, 'XData', r_rgt(kk,1), 'YData', r_rgt(kk,2), 'ZData', r_rgt(kk,3));

            drawnow
        end
    end

 % _________________ Ground track plot __________________

    figGT = figure;
    hold on
    image(lon_img, lat_img, earth);
    set(gca,'YDir','normal');
    
    if ANIMATE_GT
    
        % hGT_nom / hGT_rgt are the animated POINT CLOUD objects (ground track history)
        hGT_nom = plot(NaN,NaN, 'LineStyle','none', 'Marker','.', 'Color','g', 'MarkerSize',2,'DisplayName','Nominal ground track');
        hGT_rgt = plot(NaN,NaN, 'LineStyle','none', 'Marker','.', 'Color','r', 'MarkerSize',2,'DisplayName','RGT ground track');
    
        % hDotGT_nom / hDotGT_rgt are the moving markers (current sub-satellite point)
        hDotGT_nom = plot(NaN,NaN, 'go', 'MarkerSize',10, 'LineWidth',3, 'HandleVisibility','off');
        hDotGT_rgt = plot(NaN,NaN, 'ro', 'MarkerSize',12, 'LineWidth',4, 'HandleVisibility','off');
    
    else
        plot(lon_nom, lat_nom, 'LineStyle','none', 'Marker','.', 'Color','g', 'MarkerSize',2,'DisplayName','Nominal ground track');
        plot(lon_rgt, lat_rgt, 'LineStyle','none', 'Marker','.', 'Color','r', 'MarkerSize',2, 'DisplayName','RGT ground track');
    end
    
    % start/end markers
    plot(lon_nom(1),   lat_nom(1),   'go', 'MarkerSize',10, 'LineWidth', 3 ,'DisplayName','Nominal start');
    plot(lon_nom(end), lat_nom(end), 'gx', 'MarkerSize',12, 'LineWidth',4 , 'DisplayName','Nominal end');
    
    plot(lon_rgt(1),   lat_rgt(1),   'ro', 'MarkerSize',8,'LineWidth',3,'DisplayName','RGT start');
    plot(lon_rgt(end), lat_rgt(end), 'rx', 'MarkerSize',12,'LineWidth',4,'DisplayName','RGT end');
    
    xlabel('Longitude [deg]')
    ylabel('Latitude [deg]')
    title(['Ground Track - ' case_names{icase}])
    legend("show", 'Location','northeast')
    xlim([-180 180]); ylim([-90 90]);
    grid on
    axis equal tight

    if ANIMATE_GT
    
        NanimGT = min(length(lon_nom), length(lon_rgt));
    
        % IMPORTANT: forcing the animation to include the LAST point (NanimGT)
        frames = unique([1:STEP_ANIM:NanimGT, NanimGT]);
    
        for kk = frames
    
            set(hGT_nom, 'XData', lon_nom(1:kk), 'YData', lat_nom(1:kk));
            set(hGT_rgt, 'XData', lon_rgt(1:kk), 'YData', lat_rgt(1:kk));
    
            set(hDotGT_nom, 'XData', lon_nom(kk), 'YData', lat_nom(kk));
            set(hDotGT_rgt, 'XData', lon_rgt(kk), 'YData', lat_rgt(kk));
    
            drawnow
        end
    end
end

%% ___________________ PROPAGATION FOR THE PERTURBED CASE BY MOON AND J2____________________________

%_____________________ SETTINGS ___________________________

% 'cart'  -> propagate [r;v] with 2BP + acc_pert_fun (J2+Moon)
% 'gauss' -> propagate kep with Gauss eq (uses acc_pert_fun internally)
comando_utente = input('How you want to propagate the orbit? (default Cartesian)\nWrite -Gauss- for a Gaussian propagation\nWrite -Cart- for a cartesian propagation:\n', 's');
if isempty(comando_utente)
    comando_utente = 'Cart';
end

PROP_MODE = comando_utente;
options = odeset('RelTol',1e-13,'AbsTol',1e-13);


% _________________ GROUND TRACK __________________

% Parameters vector for perturbations
% [J2, mu, mass, Radius, t0_mjd2000]

parameters = [Primary.J2, Primary.mu, Primary.mass, Primary.Radius, Primary.Initial_date_Mjd];

% _________________ Animation settings __________________
Animation_3D = input('Animation 3D plot? (default off)\nWrite -true- for an animated time propagation\nWrite -false- for a normal plot:\n', 's');
if isempty(Animation_3D)
    Animation_3D = 'false';
end
Animation_GT = input('Animation Ground Track plot? (default off)\nWrite -true- for an animated time propagation\nWrite -false- for a normal plot:\n', 's');
if isempty(Animation_GT)
    Animation_GT = 'false';
end
ANIMATE_3D = strcmp(Animation_3D, 'true');
ANIMATE_GT = strcmp(Animation_GT, 'true');
STEP_ANIM  = 20000;

% _________________ Case selection __________________
Running = input('do you want to plot all the case or only just one?\nWrite -single- or -all-: ', 's');
if isempty(Running)
    Running = 'single';
end
RUN_MODE = Running;


T_Earth = 2*pi / Primary.w;
Pert.case_names = { ...
    '1 nominal orbital period (T_{nom})', ...
    '2 RGT orbital period (T_{rgt})', ...
    '50 RGT orbital periods (2T_{rgt})', ...
    '1 Earth rotation (T_Earth)'};
Pert.case_Tend = [ T_nom, 2*T_rgt, 50*T_rgt, T_Earth ]; 

if strcmp(RUN_MODE, 'single')

    choice = input(['Which case you want to plot?\n' ...
        '1. 1 x nominal orbital period (T_{nom})\n' ...
        '2. 2 x RGT orbital period (T_{rgt})\n' ...
        '3. 50 x RGT orbital periods (2T_{rgt})\n' ...
        '4. 1 x Earth rotation (T_Earth)\n' ...
        'Selection (default 1): '], 's');
    
    if isempty(choice)
        CASE_ID = 1; % default
    else
        CASE_ID = str2double(choice); %converting string into double
    end
    
    Pert.case_list = CASE_ID;
else
    % If RUN_MODE i 'all'
    Pert.case_list = 1:4;
end

% Simulation parameters
Npts = 200000;
options = odeset('RelTol',1e-9,'AbsTol',1e-10);

for icase = Pert.case_list

Tend = Pert.case_Tend(icase);
tspan = linspace(0, Tend, Npts);

% _________________ PROPAGATION (BOTH ORBITS PERTURBED) __________________
    
    if strcmpi(PROP_MODE,'cart')
    
        % __________________ NOMINAL ORBIT (CART + J2/MOON) __________________
        [Pert.t_nom, Pert.y_nom] = ode113(@(t,s) fun_2bp_pert(t, s, @acc_pert_fun, parameters), tspan, y0_nom, options);
        Pert.r_nom = Pert.y_nom(:,1:3);
    
        % __________________ RGT ORBIT (CART + J2/MOON) __________________
        [Pert.t_rgt, Pert.y_rgt] = ode113(@(t,s) fun_2bp_pert(t, s, @acc_pert_fun, parameters), tspan, y0_rgt, options);
        Pert.r_rgt = Pert.y_rgt(:,1:3);
    
    elseif strcmpi(PROP_MODE,'gauss')
    
        % __________________ NOMINAL ORBIT (GAUSS + J2/MOON) __________________
        kep0_nom = kep_nom.'; 
        [Pert.t_nom, Pert.kep_nom_hist] = ode113(@(t,kep) eq_motion_Gauss_J2_Moon(t, kep, @acc_pert_fun, parameters), tspan, kep0_nom, options);
    
        Pert.r_nom = zeros(length(Pert.t_nom),3);
        for iStep = 1:length(Pert.t_nom)
            [r_tmp, ~] = Kep2Car_vec(Pert.kep_nom_hist(iStep,:), Primary.mu);
            Pert.r_nom(iStep,:) = r_tmp.';
        end
    
        % __________________ RGT ORBIT (GAUSS + J2/MOON) __________________
        kep0_rgt = kep_rgt.'; 
        [Pert.t_rgt, Pert.kep_rgt_hist] = ode113(@(t,kep) eq_motion_Gauss_J2_Moon(t, kep, @acc_pert_fun, parameters), tspan, kep0_rgt, options);
    
        Pert.r_rgt = zeros(length(Pert.t_rgt),3);
        for iStep = 1:length(Pert.t_rgt)
            [r_tmp, ~] = Kep2Car_vec(Pert.kep_rgt_hist(iStep,:), Primary.mu);
            Pert.r_rgt(iStep,:) = r_tmp.';
        end
    
    else
        error('\nPROP_MODE must be ''gauss'' or ''cart'' ')
    end  

    % _________________ Ground track computation __________________
    % Ground_track returns lon/lat already in degrees and wrapped, so DO NOT unwrap/mod/jump/NaN here.
    [~, ~, Pert.lon_nom, Pert.lat_nom] = Ground_track(Pert.t_nom, Pert.r_nom, ThetaG0, Pert.t_nom(1), Primary.w);
    [~, ~, Pert.lon_rgt, Pert.lat_rgt] = Ground_track(Pert.t_rgt, Pert.r_rgt, ThetaG0, Pert.t_rgt(1), Primary.w);



   % ________________ PLOTS IN SYNC ___________________
    % Cleaning and imposing the plots I want in base on the CASE value

    % close all
    figSYNC = figure;
    set(figSYNC, 'Name', ['Case - ' Pert.case_names{icase}]);
    
    % __________________ FIGURE LAYOUT (TWO AXES MANUAL) ____________________
    
    ax3D = axes('Parent', figSYNC);
    axGT = axes('Parent', figSYNC);
    
    set(ax3D, 'Units','normalized', 'Position',[0.05 0.20 0.45 0.75]);   % [x y w h]
    set(axGT, 'Units','normalized', 'Position',[0.55 0.08 0.42 0.87]);   % [x y w h]
    
    % ________________________ 3D ORBITS (LEFT PLOT) ___________________________
    
    axes(ax3D);
    hold on
    surf(xe, ye, ze, 'CData', earth_img, 'FaceColor','texturemap', 'EdgeColor','none', 'HandleVisibility','off');    
    xlabel('x [km]')
    ylabel('y [km]')
    zlabel('z [km]')
    axis equal
    grid on
    view(3);
    title(['3D Orbits - ' Pert.case_names{icase}])
    
    % ___________________ NOMINAL ORBIT (TIME COLORED) ___________________
    if ANIMATE_3D
        h_nom = surface(NaN,NaN,NaN,NaN, 'FaceColor','none', 'EdgeColor','interp', 'LineWidth', 1.2, 'DisplayName', 'Nominal orbit (small perturbed)');
        hDot_nom = plot3(NaN,NaN,NaN, 'go', 'MarkerSize', 6, 'LineWidth', 1.5, 'HandleVisibility','off');    
    else   
        Xn = [Pert.r_nom(:,1), Pert.r_nom(:,1)];
        Yn = [Pert.r_nom(:,2), Pert.r_nom(:,2)];
        Zn = [Pert.r_nom(:,3), Pert.r_nom(:,3)];
        Cn = [Pert.t_nom, Pert.t_nom];
        surface(Xn, Yn, Zn, Cn, 'FaceColor','none', 'EdgeColor','interp', 'LineWidth', 1.2, 'DisplayName','Nominal orbit (small perturbed)');
    end
       
    % ___________________ RGT ORBIT (PERTURBED) in parula form with surface ___________________
 
    if ANIMATE_3D    
        hSurfRGT = surface(NaN,NaN,NaN,NaN, 'FaceColor','none', 'EdgeColor','interp', 'LineWidth', 2, 'DisplayName','RGT orbit (perturbed) - time colored');
        hDot_rgt = plot3(NaN,NaN,NaN, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5, 'HandleVisibility','off');  
    else
    
        Xr = [Pert.r_rgt(:,1), Pert.r_rgt(:,1)];
        Yr = [Pert.r_rgt(:,2), Pert.r_rgt(:,2)];
        Zr = [Pert.r_rgt(:,3), Pert.r_rgt(:,3)];
        Cr = [Pert.t_rgt, Pert.t_rgt];
        surface(Xr, Yr, Zr, Cr, 'FaceColor','none', 'EdgeColor','interp', 'LineWidth',2, 'DisplayName','RGT orbit (perturbed) - time colored');
    
    end    
    colormap(ax3D, "parula");    
    cb = colorbar(ax3D, 'southoutside');
    cb.Label.String = 'time [s]';
    cb.Position = [0.08 0.10 0.39 0.03]; % forcing the position    
    legend(ax3D, 'show', 'Location','northoutside')
        
    % ______________________ GROUND TRACK (RIGHT PLOT) ________________________
    
    axes(axGT);
    hold on
    
    image(lon_img, lat_img, earth);
    set(gca,'YDir','normal');
    
    xlabel('Longitude [deg]')
    ylabel('Latitude [deg]')
    title(['Ground Track - ' Pert.case_names{icase}])
    xlim([-180 180]); ylim([-90 90]);
    grid on
    axis equal tight
    
    if ANIMATE_GT
    
        hGT_nom = plot(NaN,NaN,'LineStyle','none', 'Marker','.', 'Color','g', 'MarkerSize',2,'DisplayName','Nominal ground track (perturbed)');
        hGT_rgt = plot(NaN,NaN,'LineStyle','none', 'Marker','.', 'Color','r', 'MarkerSize',2,'DisplayName','RGT ground track (perturbed)');
        hDotGT_nom = plot(NaN,NaN, 'go', 'MarkerSize',6, 'LineWidth',1.5, 'HandleVisibility','off');
        hDotGT_rgt = plot(NaN,NaN, 'ro', 'MarkerSize',6, 'LineWidth',1.5, 'HandleVisibility','off');
    
    else
    
        plot(Pert.lon_nom, Pert.lat_nom,'LineStyle','none', 'Marker','.', 'Color','g', 'MarkerSize',2,'DisplayName','Nominal ground track (perturbed)');
        plot(Pert.lon_rgt, Pert.lat_rgt,'LineStyle','none', 'Marker','.', 'Color','r', 'MarkerSize',2,'DisplayName','RGT ground track (perturbed)');
    
    end
    
    plot(Pert.lon_nom(1), Pert.lat_nom(1), 'go', 'MarkerSize',10, 'LineWidth', 3,'DisplayName','Nominal start');
    plot(Pert.lon_nom(end), Pert.lat_nom(end), 'gx', 'MarkerSize',12, 'LineWidth', 4,'DisplayName','Nominal end');
    
    plot(Pert.lon_rgt(1), Pert.lat_rgt(1), 'ro', 'MarkerSize',10,'LineWidth', 3, 'DisplayName','RGT start');
    plot(Pert.lon_rgt(end), Pert.lat_rgt(end), 'rx', 'MarkerSize',12, 'LineWidth', 4, 'DisplayName','RGT end');
    
    legend(axGT, "show", 'Location','northoutside')
        
    % =========================== SYNC ANIMATION ============================  I want them going together
    
    if ANIMATE_3D || ANIMATE_GT
    
        Nanim3D = min(size(Pert.r_nom,1), size(Pert.r_rgt,1));
        NanimGT = min(length(Pert.lon_nom), length(Pert.lon_rgt));
        NanimALL = min(Nanim3D, NanimGT);
    
        frames = unique([1:STEP_ANIM:NanimALL, NanimALL]);
    
        for kk = frames
    
            if ANIMATE_3D
                % --- Aggiornamento ORBITA NOMINALE ---
                Xn = [Pert.r_nom(1:kk,1), Pert.r_nom(1:kk,1)];
                Yn = [Pert.r_nom(1:kk,2), Pert.r_nom(1:kk,2)];
                Zn = [Pert.r_nom(1:kk,3), Pert.r_nom(1:kk,3)];
                Cn = [Pert.t_nom(1:kk), Pert.t_nom(1:kk)];
                set(h_nom, 'XData', Xn, 'YData', Yn, 'ZData', Zn, 'CData', Cn); % Aggiorno la superficie nominale
                
                set(hDot_nom, 'XData', Pert.r_nom(kk,1), 'YData', Pert.r_nom(kk,2), 'ZData', Pert.r_nom(kk,3));
            
                % --- Aggiornamento ORBITA RGT ---
                Xr = [Pert.r_rgt(1:kk,1), Pert.r_rgt(1:kk,1)];
                Yr = [Pert.r_rgt(1:kk,2), Pert.r_rgt(1:kk,2)];
                Zr = [Pert.r_rgt(1:kk,3), Pert.r_rgt(1:kk,3)];
                Cr = [Pert.t_rgt(1:kk), Pert.t_rgt(1:kk)];
                set(hSurfRGT, 'XData', Xr, 'YData', Yr, 'ZData', Zr, 'CData', Cr);
            
                set(hDot_rgt, 'XData', Pert.r_rgt(kk,1), 'YData', Pert.r_rgt(kk,2), 'ZData', Pert.r_rgt(kk,3));
            end

            if ANIMATE_GT
    
                set(hGT_nom, 'XData', Pert.lon_nom(1:kk), 'YData', Pert.lat_nom(1:kk));
                set(hGT_rgt, 'XData', Pert.lon_rgt(1:kk), 'YData', Pert.lat_rgt(1:kk));
                set(hDotGT_nom, 'XData', Pert.lon_nom(kk), 'YData', Pert.lat_nom(kk));
                set(hDotGT_rgt, 'XData', Pert.lon_rgt(kk), 'YData', Pert.lat_rgt(kk));
            end 
            drawnow % imposing the drawing while calculating instead then at teh end
        end
    end
end



%% _______________ ANALISING KEPLERIAN ELEMENT PLOTTED WITH PERTURBATIONS___________________

% _________________ defining the parameters vector to pass to the functions

parameters = [Primary.J2, Primary.mu, Primary.mass, Primary.Radius, Primary.Initial_date_Mjd];
options = odeset('RelTol',1e-13,'AbsTol',1e-13); % ODE solver options

% I find values for the state from keplerian elements
[r0, v0] = Kep2Car_vec(kep, Primary.mu);

% Initial State Cartesian vector
s0 = [r0; v0];

n_points = 200000;
t0 = 0;
n_periods = 200;

T = 2*pi*sqrt(kep(1)^3/Primary.mu);
Tspan = linspace(t0, n_periods*T, n_points);
kep0 = kep; 

% I want to calculate both form, in cartesian and keplerian form, so that i
% can find the differences in keplerian parameters perturbations

[t, cart.y] = ode89(@(t,s) fun_2bp_pert(t, s, @acc_pert_fun, parameters), Tspan, s0, options);
[t, gauss.kep] = ode89( @(t,kep) eq_motion_Gauss_J2_Moon(t, kep, @acc_pert_fun, parameters), Tspan, kep0, options);


tau = t / T;

% Post processing for the verifications

cart.r= cart.y(:,1:3);
cart.v = cart.y(:,4:6);

% Normal plot in cartesian form

X = [cart.r(:,1), cart.r(:,1)];
Y = [cart.r(:,2), cart.r(:,2)];
Z = [cart.r(:,3), cart.r(:,3)];
C = [tau, tau];

figure
hold on
surf(xe, ye, ze,'CData', earth_img,'FaceColor','texturemap',EdgeColor='none',HandleVisibility='off');
surface(X, Y, Z, C, 'FaceColor','none','EdgeColor','interp','LineWidth',2)
colormap("parula");
colorbar('southoutside');
% Forcing the position for better Latex Plots
ax  = gca;
clb = get(ax,'Colorbar');
ax.Position = [0.08 0.15 0.86 0.78];
clb.Position(4) = 0.03;                   
clb.Position(3) = 0.50;                   
clb.Position(1) = 0.25;                   
clb.Position(2) = 0.06;
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
title('Propagation in cartesian coordinates')
axis equal 
grid on
view(3)

% Plot with cartesian form propagation of J2

r_plot_try = zeros (3,length(t));
v_plot_kep_try = zeros (3,length(t));
for h = 1 : length(t)
    [r_plot_try(:,h), v_plot_kep_try(:,h)] = Kep2Car_vec(gauss.kep(h,:), Primary.mu);
end

% Building Axis to plot

X = [r_plot_try(1,:).', r_plot_try(1,:).'];  % N×2
Y = [r_plot_try(2,:).', r_plot_try(2,:).'];  % N×2
Z = [r_plot_try(3,:).', r_plot_try(3,:).'];  % N×2
C = [tau, tau];                                            % N×2

figure
surface(X, Y, Z, C, 'FaceColor','none','EdgeColor','interp','LineWidth',2);
hold on
surf(xe, ye, ze,'CData', earth_img,'FaceColor','texturemap',EdgeColor='none',HandleVisibility='off');
colormap("parula");
colorbar('southoutside');
% Forcing the position for better Latex Plots
ax  = gca;
clb = get(ax,'Colorbar');
ax.Position = [0.08 0.15 0.86 0.78];
clb.Position(4) = 0.03;                   
clb.Position(3) = 0.50;                   
clb.Position(1) = 0.25;                   
clb.Position(2) = 0.06;
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
title('Propagation with Keplerian element')
axis equal 
grid on
view(3)

% -------------------- Post processing : plots of what We found --------------------

% Converting cartesian into keplerian
kep_from_cart = zeros(length(t),6);

for k = 1:length(t)
    kep_from_cart(k,:) = car2kep(cart.r(k,:)', cart.v(k,:)', Primary.mu).';
end

% For filtering

N = length(tau);
ptsPerOrbit = round(N / n_periods);
w = ptsPerOrbit;

% wrapping function for angular measures
wrapDiff = @(A,B) atan2(sin(A-B), cos(A-B));

% measures respect which I normalize in plots
a0 = kep0(1);
normVec = [a0, 2*pi, 2*pi, 2*pi, 2*pi, 2*pi];

% Differences in abs for relative error
d = zeros(N,6);
d(:,1) = abs(kep_from_cart(:,1) - gauss.kep(:,1));
d(:,2) = abs(kep_from_cart(:,2) - gauss.kep(:,2));
d(:,3) = abs(kep_from_cart(:,3) - gauss.kep(:,3));
d(:,4) = abs(wrapDiff(kep_from_cart(:,4), gauss.kep(:,4)));
d(:,5) = abs(wrapDiff(kep_from_cart(:,5), gauss.kep(:,5)));
d(:,6) = abs(wrapDiff(kep_from_cart(:,6), gauss.kep(:,6)));

errRel = d ./ normVec; % Calculate relative errors for the differences

%build plot arrays (coerenti in unità)
kep_from_cart_plot = zeros(N,6);
gauss_kep_plot = zeros(N,6);

kep_from_cart_plot(:,1:2) = kep_from_cart(:,1:2);
gauss_kep_plot(:,1:2) = gauss.kep(:,1:2);

for j=3:6
    kep_from_cart_plot(:,j) = rad2deg(unwrap(kep_from_cart(:,j)));
    gauss_kep_plot(:,j) = rad2deg(unwrap(gauss.kep(:,j)));
end

% filter in the SAME units used for plotting
kep_from_cart_filt_plot = zeros(N,6);
for j=1:6
    kep_from_cart_filt_plot(:,j) = movmean(kep_from_cart_plot(:,j), w, 'Endpoints','shrink');
end

names  = {'a','e','i','Raan','omega','f'};
yfull  = {'a [km]','e [-]','i [deg]','Raan [deg]','omega [deg]','f [deg]'};
yerr   = {'|a_{cart}-a_{Gauss}|/a_0 [-]', ...
          '|e_{cart}-e_{Gauss}|/(2\pi) [-]', ...
          '|i_{cart}-i_{Gauss}|/(2\pi) [-]', ...
          '|\Delta\Omega|/(2\pi) [-]', ...
          '|\Delta\omega|/(2\pi) [-]', ...
          '|\Delta f|/(2\pi) [-]'};

for j = 1:6

    % --- error (semilog)
    figure('Name', ['Sample results - error - ' names{j}]);
    semilogy(tau, errRel(:,j), 'LineWidth', 1.2);
    grid on;
    axis tight;
    xlabel('time [T]');
    ylabel(yerr{j});
    title(['Error on ' names{j}]);

    % --- full evolution (0–100T)
    figure('Name', ['Sample results - ' names{j}]);
    hold on;
    plot(tau, gauss_kep_plot(:,j), 'LineWidth', 1.0);
    plot(tau, kep_from_cart_plot(:,j), 'LineWidth', 1.0);
    plot(tau, kep_from_cart_filt_plot(:,j), 'LineWidth', 1.2);
    % --- plotting the trend exploiting the interpolated polynomial, to
    % better show
    % ---- secular trend line (linear fit) ----
    p = polyfit(tau(:), kep_from_cart_filt_plot(:,j), 1);
    trend = polyval(p, tau);
    plot(tau, trend, 'LineWidth', 1);
    grid on;
    xlabel('time [T]');
    ylabel(yfull{j});
    title(['Evolution of ' names{j} ' (0–200T)']);
    legend('Gauss equations','Cartesian','Secular (filtered)','Trend','Location','south');

    % --- zoom 0–10T
    figure('Name', ['Sample results - zoom - ' names{j}]);
    hold on;
    idx10 = tau <= 10;
    plot(tau(idx10), gauss_kep_plot(idx10,j), 'LineWidth', 1.0);
    plot(tau(idx10), kep_from_cart_plot(idx10,j), 'LineWidth', 1.0);
    plot(tau(idx10), kep_from_cart_filt_plot(idx10,j), 'LineWidth', 1.2);
    grid on;
    xlabel('time [T]');
    ylabel(yfull{j});
    title(['Evolution of ' names{j} ' (0–10T)']);
    legend('Gauss equations','Cartesian','Secular (filtered)','Location','south');

end

%% SATELLITE EPJHEMERIDES ANALYSIS

% Under here the two debris/satellite used for the analisis of the
% perturbation.
% Star-Track was used to find the TLES's of the objects
% Nasa JPL Horizon System was used to track the ephemerides between the
% departures data of our project and 

% _____________________ OBJECT 1 LEO DESCRIPTION and EXPLANATION OF THE CHOICE ____________________

% Object Name: COSMOS 2251 DEB
% Parent Satellite: COSMOS 2251
% Orbit Class: Low Earth Orbit (LEO)
% Operational Status: Uncontrolled debris
% Perturbation Sensitivity: High

% COSMOS 2251_DEB refers to debris generated from the catastrophic collision of the Russian military satellite COSMOS 2251 
% with the Iridium-33 satellite in 2009. The debris cloud spans a wide range of LEO orbital parameters, 
% with objects typically characterized by:

% ---> Moderate eccentricity,
% ---> Inclinations near polar or high-inclination regimes,
% ---> Rapid orbital evolution due to environmental effects.

% COSMOS 2251 debris fragments are completely 
% uncontrolled and subject to strong perturbations, including:

% Atmospheric drag, which dominates long-term decay and semi-major axis reduction,
% J2 perturbation, producing rapid nodal regression,
% Third-body perturbations, which are less dominant than drag but still present,
% Solar radiation pressure (for high area-to-mass ratio fragments).

% These objects are particularly useful as a contrast case with respect to MEO satellites:
% They highlight the importance of non-conservative forces,
% They exhibit fast divergence of Keplerian elements,
% They are ideal for demonstrating the limits of simplified perturbation models.


% _____________________ OBJECT 2 MEO DESCRIPTION and EXPLANATION OF THE CHOICE _____________________

% NORAD ID: 25933
% Common Name: GPS BIIR-10 (NAVSTAR)
% Orbit Class: Medium Earth Orbit (MEO)
% Operational Status: Nominally controlled during mission lifetime, currently suitable as a reference MEO object
% Inclination: ~55°
% Semi-major Axis: ~26,560 km
% Eccentricity: Near-circular
% 
% Description and relevance
% 
% NORAD 25933 corresponds to the GPS Block IIR-10 satellite, part of the NAVSTAR Global Positioning System constellation. 
% It operates in a Medium Earth Orbit (MEO) with a semi-major axis of approximately 26,560 km and an inclination close to 55°, 
% which is characteristic of the GPS constellation geometry.
% 
% From the point of view of orbital perturbation analysis, this satellite represents an excellent reference case:
% Atmospheric drag is negligible, making it ideal for isolating gravitational perturbations.
% 
% The dominant long-term perturbations are due to:
% Earth's oblateness (J2), responsible for secular precession of RAAN and argument of perigee.
% Third-body effects, primarily from the Moon (and secondarily the Sun), 
% which introduce long-period and secular variations in inclination and eccentricity.
% The relatively high altitude allows perturbation effects to accumulate over long timescales, 
% making Keplerian element variations clearly observable.
% 
% For these reasons, NORAD 25933 is well suited for:
%  --> Validation of Gauss planetary equations,
%  --> Comparison between Cartesian propagation and Keplerian propagation with perturbations,
%  --> Analysis of repeated ground track degradation under non-Keplerian dynamics.


% _____________________ OBJECT 1 LEO DESCRIPTION and EXPLANATION OF THE CHOICE ____________________

A = importdata('COSMOS2251_2010.txt', ',');

COSMOS.evec=A.data(:,1);
COSMOS.ivec=deg2rad(A.data(:,3));
COSMOS.Raanvec=deg2rad(A.data(:,4));
COSMOS.omegavec=deg2rad(A.data(:,5));
COSMOS.theta0vec=deg2rad(A.data(:,9));
COSMOS.avec=A.data(:,10);
COSMOS.e0=A.data(1,1);
COSMOS.i0=deg2rad(A.data(1,3));
COSMOS.Raan0=deg2rad(A.data(1,4));
COSMOS.omega0=deg2rad(A.data(1,5));
COSMOS.theta00=deg2rad(A.data(1,9));
COSMOS.a0=A.data(1,10);
COSMOS.kep0 = [COSMOS.a0, COSMOS.e0, COSMOS.i0, COSMOS.Raan0, COSMOS.omega0, COSMOS.theta00];

[COSMOS.r0,COSMOS.v0] = Kep2Car_vec(COSMOS.kep0, Primary.mu);

r=[];
v=[];
% plotting on one hour

% converting also in cartesian, to do for both the comparation of orbital
% elements and orbit in Gaussian form

N = length(COSMOS.evec);

dt = 3600;
COSMOS.tvec = (0:N-1) * dt;

COSMOS.r = zeros(N,3);
COSMOS.v = zeros(N,3);

for j=1:N
    [COSMOS.r(j,:), COSMOS.v(j,:)] = Kep2Car(COSMOS.avec(j), COSMOS.evec(j), COSMOS.ivec(j),COSMOS.Raanvec(j), COSMOS.omegavec(j), COSMOS.theta0vec(j), Primary.mu);
end

X = [COSMOS.r(:,1), COSMOS.r(:,1)];
Y = [COSMOS.r(:,2), COSMOS.r(:,2)];
Z = [COSMOS.r(:,3), COSMOS.r(:,3)];

T = 2*pi*sqrt(COSMOS.avec(1)^3/Primary.mu);

C = repmat(COSMOS.tvec(:)/T, 1, 2);

figure
hold on
surf(xe, ye, ze,'CData', earth_img,'FaceColor','texturemap',EdgeColor='none',HandleVisibility='off');
surface(X, Y, Z, C, 'FaceColor','none','EdgeColor','interp','LineWidth',2)

colormap("parula");
cb = colorbar;
cb.Label.String = 'Numero di periodi (t/T)';

xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
title('Propagation in cartesian coordinates')
axis equal
grid on
view(3)

kepCOSMOS = [COSMOS.avec, COSMOS.evec, COSMOS.ivec, COSMOS.Raanvec, COSMOS.omegavec, COSMOS.theta0vec];
plotKepFiltered(COSMOS.tvec, kepCOSMOS, Primary.mu, "secular");




%% _____________________ OBJECT 2 MEO DESCRIPTION and EXPLANATION OF THE CHOICE _____________________


A = importdata('NORAD25933.txt', ',');

GPS_BIIR.evec      = A.data(:,1);
GPS_BIIR.ivec      = deg2rad(A.data(:,3));
GPS_BIIR.Raanvec   = deg2rad(A.data(:,4));
GPS_BIIR.omegavec  = deg2rad(A.data(:,5));
GPS_BIIR.theta0vec = deg2rad(A.data(:,9));
GPS_BIIR.avec      = A.data(:,10);

GPS_BIIR.e0      = A.data(1,1);
GPS_BIIR.i0      = deg2rad(A.data(1,3));
GPS_BIIR.Raan0   = deg2rad(A.data(1,4));
GPS_BIIR.omega0  = deg2rad(A.data(1,5));
GPS_BIIR.theta00 = deg2rad(A.data(1,9));
GPS_BIIR.a0      = A.data(1,10);

GPS_BIIR.kep0 = [GPS_BIIR.a0, GPS_BIIR.e0, GPS_BIIR.i0, GPS_BIIR.Raan0, GPS_BIIR.omega0, GPS_BIIR.theta00];

[GPS_BIIR.r0,GPS_BIIR.v0] = Kep2Car_vec(GPS_BIIR.kep0, Primary.mu);

r=[];
v=[];

dt = 3600;
GPS_BIIR.tvec = (0:length(GPS_BIIR.evec)-1) * dt;

for j=1:length(GPS_BIIR.evec)
     [GPS_BIIR.r(j,:), GPS_BIIR.v(j,:)] = Kep2Car(GPS_BIIR.avec(j,1), GPS_BIIR.evec(j,1), GPS_BIIR.ivec(j,1), GPS_BIIR.Raanvec(j,1), GPS_BIIR.omegavec(j,1), GPS_BIIR.theta0vec(j,1), Primary.mu);
end

X = [GPS_BIIR.r(:,1), GPS_BIIR.r(:,1)];
Y = [GPS_BIIR.r(:,2), GPS_BIIR.r(:,2)];
Z = [GPS_BIIR.r(:,3), GPS_BIIR.r(:,3)];

T = 2*pi *sqrt(GPS_BIIR.avec(1)^3/Primary.mu);
C = [GPS_BIIR.tvec'/T GPS_BIIR.tvec'/T];

figure
hold on
surf(xe, ye, ze,'CData', earth_img,'FaceColor','texturemap',EdgeColor='none',HandleVisibility='off');
surface(X, Y, Z, C, 'FaceColor','none','EdgeColor','interp','LineWidth',2)
colormap("parula");
cb = colorbar;
cb.Label.String = 'Numero di periodi (t/T)';
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
title('Propagation in cartesian coordinates')
axis equal 
grid on
view(3)

kepGPS_BIIR = [GPS_BIIR.avec, GPS_BIIR.evec, GPS_BIIR.ivec, GPS_BIIR.Raanvec, GPS_BIIR.omegavec, GPS_BIIR.theta0vec];
plotKepFiltered(GPS_BIIR.tvec, kepGPS_BIIR, Primary.mu, "secular");


% Sec
% NORAD 25933
% (GPS BIIR-2, ormai spesso usato come caso studio)




% ______________ SATELLITE SCENARIO FOR PLOT TO IMPLEMENT ______________


% Implementing the satellite scenario

% %% ----------- SATELLITE SCENARIO ----------------
% 
% %Conversion a km -> m
% a = a*1000;         %[m]
% startTime = 
% 
% %Satellite Scenario 
% stopTime = startTime + seconds(T);
% sampleTime = 60;
% sc = satelliteScenario(startTime,stopTime,sampleTime);
% viewer = satelliteScenarioViewer(sc,"CameraReferenceFrame","Inertial","Dimension","3D");
% sat = satellite(sc,a,e,incl,raan,w,theta);
% show(sat)
% sat.Visual3DModel = "SmallSat.glb";
% coordinateAxes(sat, Scale=2); % red = x_B; green = y_B; blue = z_B
% camtarget(viewer, sat);
% groundTrack(sat,"LeadTime",3600,"LeadLineColor",[0 1 0],"TrailLineColor",[0 1 0]);
% play(sc,PlaybackSpeedmultiplier=500)
% 
% %Conversion a m -> km
% a = a/1000;         %[km] 


