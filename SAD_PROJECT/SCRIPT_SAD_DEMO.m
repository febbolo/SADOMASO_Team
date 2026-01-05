%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%                                                         %%%%%%%%%
%%%%%%%%%           SPACECRAFT ATTITUDE DYNAMICS PROJECT          %%%%%%%%%
%%%%%%%%%                       A.A. 2025-2026                    %%%%%%%%%
%%%%%%%%%                                                         %%%%%%%%%
%%%%%%%%%                  GROUP 44, PROJECT N.151                %%%%%%%%%
%%%%%%%%%                                                         %%%%%%%%%
%%%%%%%%%                                                         %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CONTRIBUTORS
% Amura Fabio, 10830618, fabio.amura@mail.polimi.it
% Crisanti Costanza, 10911209, costanza.crisanti@mail.polimi.it
% Deli Luca 
% Tomás Fadista, 11027292, tomas.nascimentodasilva@mail.polimi.it

%% ------------ DATA : Earth-S/C - Sun synchronous Orbit ------------

clear
clc
close all

% S/C Inertia 
% 6U CubeSat inertia
J_depl = [32.08, 0, 0;...
          0, 14.21, 0;...
          0, 0, 44.67]*1e-2; %kg*m^2

% Constants
G = 6.67e-20; % [km^3/(kg*s^2)]
Mt = 5.97e24;   % [kg]
Re = astroConstants(23); %[km]
J2 = 0.00108263; % Second Zonal harmonic
mu = astroConstants(13); %[km^3/s^2]

% Date Time 
% Choose the spring equinox so that the Earth sees the Sun along the x
% inertial direction, almost same of x_t axis of the target reference frame
startTime = datetime(2025,3,21,0,0,0);   % Initial time epoch of project

% Orbital Parameters
e = 0.005; % Low eccentricity orbit
draan = 1.99096871e-7; % Sun-synchronous rate of precession of raan [rad/s]
zp = 600; %[km]
rp = zp + Re; %[km]
a = rp * (1 - e);  % Semi-major axis [km]
raan = 90;   %[deg]
w = 0;     %[deg]
theta = 0;  %[deg]
incl = rad2deg( acos( -2/3 * (1-e^2)^2*a^(7/2) * draan / ( sqrt(mu)*J2*Re^2 ) )); %[deg]

% Mean motion
n = sqrt( (mu/ (a)^3));   %[rad/s]

%Orbital period
T = 2*pi/n;     %[s]

% Initial Conditions
w0 = [0; 0; 0];  %[rad/s]

% Creating initial condition Keplerian elements vector
kep = [a,e,incl,raan,w,theta];

%% ---------------- Sun-Earth -----------------

% Revolution Solar Year Period
T_sun = 365*24*60*60;        % [s]

% Computing mean motion of Earth around the Sun
n_sun = 2*pi/T_sun;         %[rad/s]

% Obliquity of Equatorial Plane
eps_E = deg2rad(astroConstants(63)); %[rad]

% Astronomic unit
AU = astroConstants(2);   %[km]

% Distance of the Sun
r_sun = 1*AU;       %[km]



%% ------------- Solar Radiation Pressure (SRP) Torque ----------

F_e = 1358; %[W/m^2]
c = astroConstants(5)*1e3;  %[m/s]

SRP.body1.N_hat = [1,0,0];
SRP.body1.rho_s = 0.5;
SRP.body1.rho_d = 0.1;
SRP.body1.A = 6e-2;            %[m^2]
SRP.body1.r_F = [10,0,0]*1e-2; %[m]

SRP.body2.N_hat = [0,1,0];
SRP.body2.rho_s = 0.5;
SRP.body2.rho_d = 0.1;
SRP.body2.A = 6e-2;            %[m^2]
SRP.body2.r_F = [0,10,0]*1e-2; %[m]

SRP.body3.N_hat = [-1,0,0];
SRP.body3.rho_s = 0.5;
SRP.body3.rho_d = 0.1;
SRP.body3.A = 6e-2;            %[m^2]
SRP.body3.r_F = [-10,0,0]*1e-2; %[m]

SRP.body4.N_hat = [0,-1,0];
SRP.body4.rho_s = 0.5;
SRP.body4.rho_d = 0.1;
SRP.body4.A = 6e-2;            %[m^2]
SRP.body4.r_F = [0,-10,0]*1e-2;       %[m]

SRP.body5.N_hat = [0,0,1];
SRP.body5.rho_s = 0.5;
SRP.body5.rho_d = 0.1;
SRP.body5.A = 4e-2;            %[m^2]
SRP.body5.r_F = [0,0,+15]*1e-2;       %[m]

SRP.body6.N_hat = [0,0,-1];
SRP.body6.rho_s = 0.5;
SRP.body6.rho_d = 0.1;
SRP.body6.A = 4e-2;            %[m^2]
SRP.body6.r_F = [0,0,-15]*1e-2; %[m]

SRP.panel1.N_hat = [+1,0,0];
SRP.panel1.rho_s = 0.1;
SRP.panel1.rho_d = 0.1;
SRP.panel1.A = 12e-2;            %[m^2]
SRP.panel1.r_F = [0,45,0]*1e-2; %[m]

SRP.panel2.N_hat = [-1,0,0];
SRP.panel2.rho_s = 0.1;
SRP.panel2.rho_d = 0.1;
SRP.panel2.A = 12e-2;            %[m^2]
SRP.panel2.r_F = [0,45,0]*1e-2; %[m]

SRP.panel3.N_hat = [+1,0,0];
SRP.panel3.rho_s = 0.1;
SRP.panel3.rho_d = 0.1;
SRP.panel3.A = 12e-2;            %[m^2]
SRP.panel3.r_F = [0,-45,0]*1e-2; %[m]

SRP.panel4.N_hat = [-1,0,0];
SRP.panel4.rho_s = 0.1;
SRP.panel4.rho_d = 0.1;
SRP.panel4.A = 12e-2;            %[m^2]
SRP.panel4.r_F = [0,-45,0]*1e-2; %[m]

% Note that SRP torque DOES NOT DEPEND on initial conditions on omega, or
% satellite properties in terms of inertia. It only depends on the distance
% from the barycenter and on the extension of the cross section.
% If we consider a symmetric S/C in terms of reflective surfaces, panels ecc.. we obtain a zero torque
% The only way to change the torque is to consider a non-symmetric S/C.

%% ---------  ALbedo / Earth Radiation Models - DEMO -----------

P=F_e/c;

Nmat=[
    1 0 0;
    0 1 0;
    -1 0 0;
    0 -1 0;
    0 0 1;
    0 0 -1;
    1 0 0;
    -1 0 0;
    1 0 0;
    -1 0 0;
];


rhos=zeros(10,1);
rhos(1:6,1)=0.5;
rhos(7:10)=0.1;

rhod=0.1*ones(10,1);

Aree=zeros(10,1);
Aree(1:4)=6e-2;
Aree(5:6)=4e-2;
Aree(7:10)=12e-2;


rf=[
    10e-2 0 0;
    0 10e-2 0;
    -10e-2 0 0;
    0 -10e-2 0;
    0 0 15e-2;
    0 0 -15e-2;
    0 45e-2 0;
    0 45e-2 0;
    0 -45e-2 0;
    0 -45e-2 0;
];

F0 = 237; %[W/m^2]

%% ---------- GG torque (for inclined and elliptical orbits) ------------

% Adimensional coefficients for stability
I_moments = diag(J_depl);
K_yaw = (I_moments(3)-I_moments(2))/I_moments(1)
K_roll = (I_moments(3)-I_moments(1))/I_moments(2)
K_pitch = (I_moments(2)-I_moments(1))/I_moments(3)


% STABILITY CONDITIONS HERE ARE FOR NADIR POINTING CASE MAYBE : VERIFY!!!!
% Verify Stability conditions
if((K_roll*K_yaw>0) && ( (1+3*K_roll+K_roll*K_yaw)^2 > 16*K_yaw*K_roll) )
    disp('Stability conditions are satisfied.');
else
    disp('Stability conditions are not satisfied.');
end


%% -------- Magnetic Torque (IGRF 2025 Model) - Optimization Cycle for B ------

g_tab = readmatrix("IGRF14_g_coeffs_2025.csv");
h_tab = readmatrix("IGRF14_h_coeffs_2025.csv");

% Cycle from N = 2 
for N = 1:13
    
end 

%% ------------- Magnetic Torque (IGRF 2025 Model) ----------------

JD0 = juliandate(startTime);            % Conversion in Julian Date

j_B = [0.01; 0.05; 0.01]; % [A*m^2], dipole vector of S/C
w_E = deg2rad(15.04)/3600; %[rad/s], average angular velocity of the Earth 

g_tab = readmatrix("IGRF14_g_coeffs_2025.csv");
h_tab = readmatrix("IGRF14_h_coeffs_2025.csv");
N = 13;  % Order of IGRF Model

g = zeros(N, N+1);   % g coefficients
for ii = 1:size(g_tab,1)
    nn = g_tab(ii,1);
    mm = g_tab(ii,2);
    g(nn, mm+1) = g_tab(ii,3);
end

h = zeros(N, N+1);   % h coefficients
for ii = 1:size(h_tab,1)
    nn = h_tab(ii,1);
    mm = h_tab(ii,2);
    h(nn, mm+1) = h_tab(ii,3);
end

clear nn
clear mm
clear g_tab
clear h_tab
t0 = 0;      %[s]

%% -------------- SENSORS ------------------

% Seed 
base_seed = 12345;

% Sun Sensor : Solar MEMS nanoSSOC-D60 
% Sampling Time from data sheet 
Sun_sensor.f = 5;     %Sampling frequency [Hz]
Sun_sensor.Ts = 1/Sun_sensor.f; %Sampling time [s]
Sun_sensor.R = [1, deg2rad(0.1), -deg2rad(0.1);
                -deg2rad(0.1), 1, deg2rad(0.1);
                deg2rad(0.1), -deg2rad(0.1), 1]; % Misalignment Error Matrix [rad]
Sun_sensor.D = 0.1; % Noise density [u/sqrt(Hz)]
Sun_sensor.b = 0.05; %[deg]
Sun_sensor.AD_nbit = 16; % Number of bits of A/D 
Sun_sensor.FS = 120;     % [deg]
Sun_sensor.LSB = Sun_sensor.FS/ (2^(Sun_sensor.AD_nbit));  %[levels]
Sun_sensor.FOV = [-60 60]; %[deg]
Sun_sensor.lat = 1;      %[sample]

% Magnetometer : BOSCH BMM350DS-001
Magmeter.f = 5;         %Sampling frequency [Hz]
Magmeter.Ts = 1/Magmeter.f; %Sampling time [s]
Magmeter.O = [1, deg2rad(0.1), -deg2rad(0.1);
                -deg2rad(0.1), 1, deg2rad(0.1);
                deg2rad(0.1), -deg2rad(0.1), 1]; % Non-orthogonality Error Matrix [rad]
Magmeter.R = [1, deg2rad(0.1), -deg2rad(0.1);
                -deg2rad(0.1), 1, deg2rad(0.1);
                deg2rad(0.1), -deg2rad(0.1), 1]; % Misalignment Error Matrix [rad]
Magmeter.Dx = 27e-9; % Noise density on x [nT/sqrt(Hz)] (data sheet)
Magmeter.Dy = 27e-9; % Noise density on y [nT/sqrt(Hz)] (data sheet)
Magmeter.Dz = 64e-9; % Noise density on x [nT/sqrt(Hz)] (data sheet)
Magmeter.b = 2e-6; % Bias Error [ from muT -> to T] (data sheet)
Magmeter.AD_nbit = 16; % Number of bits of A/D (data sheet)
Magmeter.FS = [-2000e-6 2000e-6]; % Maximum range [from muT -> to T] (data sheet)
Magmeter.LSB = ( Magmeter.FS(2) - Magmeter.FS(1)) / (2^(Magmeter.AD_nbit)); % Least Significant Bit [levels]
Magmeter.lat = 5; % Latency [sample]
Magmeter.SFNx = 2.5; % Scale factor Nonlinearity on x [1/T]
Magmeter.SFNy = 2.5; % Scale factor Nonlinearity on y [1/T]
Magmeter.SFNz = 5; % Scale factor Nonlinearity on z [1/T]

% Gyroscope 
Gyro.f = 10; %[Hz]
% Full scale range +- 250 deg/s
Gyro.Ts = 1 / Gyro.f;  % Calculate the sampling time based on frequency
% Full scale range
Gyro.half_scale_range = 250;  % Full scale range in deg/s
Gyro.FS = Gyro.half_scale_range*2;
Gyro.AD_nbit = 16 ;
% Nonlinearity, in percentage (0.001 * value)
Gyro.NL = 0.1;
% Noise Density Gyro
Gyro.D = 0.0028;  % Noise density [1deg/s/sqrt(Hz)]
% Bias error
Gyro.b = 0.05;  %[deg/s] Bias gyro 
% Correlation Time for Allan variance (bias instability, Brown Noise)
Gyro.ctime = 125; % [s], from literature for cubesat MEMS
% Misalignement angles
Gyro.alpha = deg2rad(0.8); %[rad]
Gyro.beta = deg2rad(0.5);  %[rad]
Gyro.gamma = deg2rad(0.9); %[rad]
% Misalignement rotation matrix around x
Gyro.Rx= [1, 0, 0;...
          0, cos(Gyro.alpha), sin(Gyro.alpha);...
          0, -sin(Gyro.alpha), cos(Gyro.alpha)];
% Misalignement rotation matrix around y
Gyro.Ry = [cos(Gyro.beta), 0, sin(Gyro.beta);...
            0, 1, 0;... 
           -sin(Gyro.beta), 0, cos(Gyro.beta)];
% Misalignement rotation matrix around z
Gyro.Rz = [cos(Gyro.gamma), sin(Gyro.gamma), 0;...
           -sin(Gyro.gamma), cos(Gyro.gamma) 0;...
           0, 0, 1];
% Misalignement DCM Error Matrix [rad]
Gyro.R = Gyro.Rx * Gyro.Ry* Gyro.Rz;
% Non Orthogonality Error Matrix [rad]
Gyro.O = [1, deg2rad(0.002), deg2rad(0.0025);...
          deg2rad(0.001), 1, deg2rad(0.0022);...
          deg2rad(0.0018), deg2rad(0.002), 1];
% Gyro Delay (neglected since very very small)
Gyro.delay = 0;
% Gyro LSB 
Gyro.LSB = Gyro.FS / (2^Gyro.AD_nbit);


%% ------------- ATTITUDE DETERMINATION : q-method ----------

Magmeter.accuracy = 0.00475;        % [deg]
Magmeter.weight = 1/Magmeter.accuracy;

Sun_sensor.accuracy = 0.01*pi/180;  % [deg]
Sun_sensor.weight = 1/Sun_sensor.accuracy;

Magmeter.alpha = Magmeter.weight / (Magmeter.weight + Sun_sensor.weight); %[-]
Sun_sensor.alpha = Sun_sensor.weight / (Magmeter.weight + Sun_sensor.weight); %[-]


q.maxsensor = max(Sun_sensor.Ts, Magmeter.Ts);
q.Ts = max(q.maxsensor);

% --------------- FILTER ON GYROSCOPE --------------------------

% Filtering in the AD in order to not accumulate errors on control

% GYROSCOPE FILTER
% Sampling time : Ts = 0.1 s
% Sampling frequency : fc = 1/Ts = 10 Hz
% In rad/s : omega_c = 10 Hz * 2pi = 62.8 rad/s
% Nyquist frequency in rad/s : omega_N <= omega_c/2 : omega_N = 31.4 rad/s
% Frequency of cut-off (taglio) : omega_t = omega_N / k, k simulates the
% real filter and security factor . k = 2 
% --> omega_t = 15.7 rad/s

% Define a varying cutoff frequency for different phasis of mission
% For detumbling, high angular velocities --> conservative cutoff freq
Filter.DETUMBLING.omega_t = 10; %[rad/s]
% For pointing, small or almost null angular velocities --> more aggressive
% cutoff frequency to try to avoid all noises and disturbances
Filter.POINTING.omega_t = 1; %[rad/s]

% The OBC tries to identify the phases of the mission setting a treshold on
% the angular velocities : as a limit 2-3 deg/s can be considered tumbling
% In rad/s : 3*pi/180 = 0.05
% If the maximum angular rate component detected is above the treshold : satellite is
% tumbling, otherwise can be considered as slew-pointing phase

%% ----------- CONTROL : LQR ----------

% From state space model
A = [0, 0, 0, 0, 0, 0;...
     0, 0, 0, 0, 0, 0;...
     0, 0, 0, 0, 0, 0;...
     1, 0, 0, 0, 0, 0;...
     0, 1, 0, 0, 0, 0;...
     0, 0, 1, 0, 0, 0];
 % From state space model
B = [1/J_depl(1), 0, 0;...
     0, 1/J_depl(5), 0;...
     0, 0, 1/J_depl(9);...
     0, 0, 0;...
     0, 0, 0;...
     0, 0, 0];

% Checking controllability 
C = ctrb(A,B);
% Display the rank of matrix C
rank_C = rank(C);
disp(['The rank of matrix C is: ', num2str(rank_C),', equal to the number of states.' ...
    ' So the system is controllable.']);

% Checking observability
C_obs = eye(6); %Output matrix C_obs
Ob = obsv(A, C_obs);
% Display the rank of matrix Ob
rank_Ob = rank(Ob);
n_states = size(A, 1);
fprintf('\n--- Observability ---\n');
fprintf('Numero of states: %d\n', n_states);
fprintf('Rank of observability matrix: %d\n', rank_Ob);
if rank_Ob == n_states
    disp('Completely observable system');
else
    disp('Non observable system');
    disp('Verify availabilty of sensors or C_obs matrix definition');
end

% DEFINING THRESHOLDS FOR THE DIFFERENT CASES 
Control.w_tumbling = deg2rad(5); % Treshold for tumbling, [rad]
Control.w_pointing = deg2rad(0.5);     % Treshold for pointing, [rad]
% In the middle we have slew manouvre 

% ----- Angular velocities thresholds for State Flow ------
w_enter = Control.w_tumbling; % To exit DETUMBLING and enter in SLEW
w_exit = deg2rad(4);  % If exceeded, from SLEW come back to DETUMBLING
w_enter_2 = Control.w_pointing; % To exit SLEW and enter POINTING
w_exit_2 = deg2rad(1.5); % If exceeded, from POINTING come back to SLEW

% ------- Thresholds for State Flow ------
err_max    = deg2rad(15);   % angles error for POINTING
eps_w      = deg2rad(0.1);           % tolerance on |ω|
eps_att    = deg2rad(0.2);   % tolerance on attitude error
eps_h_r = 0.001;              % tolerance on h_r [Nms]

% ------- Dwell time ---------
dwell_time = 30; % [s] time before switching to another MODE if conditions are satisfied


% ------------ DE-TUMBLING ---------------
% B_dot method
% B_dot_m = -w_est x B_m
% D_required = -k_b*B_dot_m
% M_c = D_required x B = (-k_b*B_dot_m) x B_m

% Gain of the B_dot : TUNING  
k_b = 200000;

% -------- STABILIZATION / TRACKING -------------

% Bryson Method to define Q,R
% Arbitrary choice of maximum acceptable error : ~15 deg
alpha_max = deg2rad(15); %[rad]
% Since the slew manouvre starts right after the slew, we consider as
% omega_max the treshold on the slew case
% Defining desired torques:
M_RW_max = 8e-3; % [Nm]
M_M_max = 1e-4; % [Nm]

x_max = [Control.w_pointing,Control.w_pointing, Control.w_pointing, alpha_max, alpha_max, alpha_max]; 
u_max = [M_M_max, M_M_max, M_RW_max];

Q = diag(1./(x_max).^2); 
R = diag(1./(u_max).^2);

% Considering a long time of transient of P(t) matrix, P(t)≃ cost ->
% algebraic Riccati equation 
% The control is working at discrete time -> dlqr command is needed
% To avoid K matrix too small (1e-20 order), the information about the
% sampling time should given through the system state sys = ss(A, B, C, D, ts)
Ts = 0.2; % Sampling time from sensors information
sysc = ss(A, B, [], []); 
sysd = c2d(sysc, Ts, 'zoh'); 
[K,S,P] = dlqr(sysd.A,sysd.B,Q,R);    
% [K,S,P] = lqr(A,B,Q,R); 
 

% -------SLEW MANOEUVRE / RE-POINTING------------
% At the beginning, there are big angles variations -> an extension of
% linear control is needed, using THE SAME GAIN MATRIX found before

Kp = K(:, 4:6); % Proportional gain matrix dim.(3x3)
Kd = K(:, 1:3); % Derivative gain matrix dim.(3x3)

% The only effectively different part is in the simulink model, where the state vector components 
% are substituted with the correct terms due to non-small angles. 


% -------- RW DESATURATION ----------
h_r_max = 0.05;             % RW max angular momentum [Nms]
h_r_safe = 0.001;              % RW safe angular momentum [Nms]

%% ------------ ACTUATORS -----------------

% Magnetorquers - CR0020 (X-Y axis), CR0010 (Z axis) (CubeSpace)
MTQ.D_max_XY = 2; % Am^2
MTQ.D_max_Z = 1; % Am^2

% Reaction Wheel - RW400 (AAC Clyde Space)
RW.T_max = 0.008; % N*m
RW.h_max = 0.05; % N*m*s
RW.omega_max = 5000 * (2*pi/60); % rad/s (5000rpm)
RW.Iw_est = RW.h_max / RW.omega_max; % kg*m^2

% A Matrix for RW (vector, 1 RW on z axis)
RW.A = [0;0;1];


%% -------------- SIMULATION  ------------

simout = sim('SIM_SAD_DEMO');
time = simout.tout;
discrete_time_row = time(1):q.Ts:time(end);
discrete_time = discrete_time_row';
A_B_N = simout.A_B_N;
A_B_N_ortho = simout.A_B_N_ortho;
A_B_LVLH = simout.A_B_LVLH;
w_B_LVLH = simout.w_B_LVLH;

% Real Angular velocity vector (continuous) 
w_B = simout.w_B;
% Estimated angular velcity vector (discrete)
w_est = simout.w_est; 
% Estimated Attitude Error Matrix (discrete)  
A_e = simout.A_e;
% State Space vector (discrete)
x = simout.x;
% Ideal control torque (discrete) 
u = simout.u;
% Actuators real torque (continuous)
M_act = simout.M_act;
% Quaternion's attitude error (discrete)
q_err = simout.q_err;

%% ---------- PLOTS -------------

% Pre-allocate Q, error and norm_error
Q = zeros(3, 3, length(time));
error = zeros(3, 3, length(time));
norm_error = zeros(1, length(time));

% Plot the attitude error matrix components over time
figure('Name','Attitude Error Matrix Components');
for i = 1:3
    for j = 1:3
        subplot(3, 3, (i-1)*3 + j);
        plot(discrete_time, squeeze(A_e(i, j, :)), 'LineWidth', 1);
        title(['A_{e,' num2str(i) num2str(j) '} over Time']);
        xlabel('Time (s)');
        ylabel(['A_{e,' num2str(i) num2str(j) '}']);
        grid on;
    end
end

% Testing Orthonormality of A_e
for i = 1:length(discrete_time)
    Q(:,:,i) = A_e(:,:,i)'*A_e(:,:,i);
    error(:,:,i) = abs(eye(3)-Q(:,:,i));
    norm_error(i) = norm( Q(:,:,i) - eye(3), 'fro' );
end
% Plot orthonormality error of A_e
figure('Name','Orthonormality error of A_e')
plot(time, norm_error)


% Plot the real and estimated angular velocity components over time
figure('Name','Real and Estimated Angular Velocity Components');
subplot(2, 1, 1);
hold on;
plot(time, w_B(1, :), 'LineWidth', 1, 'DisplayName', 'Real w_x');
plot(time, w_B(2, :), 'LineWidth', 1, 'DisplayName', 'Real w_y');
plot(time, w_B(3, :), 'LineWidth', 1, 'DisplayName', 'Real w_z');
title('Real Angular Velocity Components');
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend('show');
grid on;
hold off;
subplot(2, 1, 2);
hold on;
plot(time, w_est(1, :), 'LineWidth', 1, 'DisplayName', 'Estimated w_x');
plot(time, w_est(2, :), 'LineWidth', 1, 'DisplayName', 'Estimated w_y');
plot(time, w_est(3, :), 'LineWidth', 1, 'DisplayName', 'Estimated w_z');
title('Estimated Angular Velocity Components');
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend('show');
grid on;
hold off;


% Plot the state vector components over time
figure('Name','State Vector Components');
for i = 1:3
    subplot(2, 3, i);
    plot(discrete_time, x(i, :), 'LineWidth', 1);
    title(['Angular Velocity Component \omega_{' num2str(i) '} over Time']);
    xlabel('Time (s)');
    ylabel(['\omega_{' num2str(i) '} (rad/s)']);
    grid on;
end
for i = 4:6
    subplot(2, 3, i);
    plot(discrete_time, x(i, :), 'LineWidth', 1);
    title(['Angle Component \theta_{' num2str(i-3) '} over Time']);
    xlabel('Time (s)');
    ylabel(['\theta_{' num2str(i-3) '} (rad)']);
    grid on;
end

% Plot the ideal control torque components over time
figure('Name','Control Ideal Torque Components');
% u_x
subplot(3,1,1);
plot(discrete_time, u(:,1), 'LineWidth', 1);
title('Ideal Torque u_x');
xlabel('Time (s)');
ylabel('Torque (N·m)');
grid on;
% u_y
subplot(3,1,2);
plot(discrete_time, u(:,2), 'LineWidth', 1);
title('Ideal Torque u_y');
xlabel('Time (s)');
ylabel('Torque (N·m)');
grid on;
% u_z
subplot(3,1,3);
plot(discrete_time, u(:,3), 'LineWidth', 1);
title('Ideal Torque u_z');
xlabel('Time (s)');
ylabel('Torque (N·m)');
grid on;

% Plot actuator torque components over time
figure('Name','Actuator Torque Components');
% M_x
subplot(3,1,1);
plot(time, squeeze(M_act(1,:)), 'LineWidth', 1);
title('Actuator Torque M_x');
xlabel('Time (s)');
ylabel('Torque (N·m)');
grid on;
% M_y
subplot(3,1,2);
plot(time, squeeze(M_act(2,:)), 'LineWidth', 1);
title('Actuator Torque M_y');
xlabel('Time (s)');
ylabel('Torque (N·m)');
grid on;
% M_z
subplot(3,1,3);
plot(time, squeeze(M_act(3,:)), 'LineWidth', 1);
title('Actuator Torque M_z');
xlabel('Time (s)');
ylabel('Torque (N·m)');
grid on;

% Plot poiting error over time
q4_err = q_err(:,4);
theta_deg = rad2deg(2*acos(q4_err)) - 180;

figure('Name','Pointing error over time');
plot(discrete_time, theta_deg, 'LineWidth', 1);
title('Pointing Error vs Time');
xlabel('Time (s)');
ylabel('\Delta\theta_{point} [deg]');
grid on;


%% ---------- MONTE CARLO ANALYSIS ----------

N_MC = 50;

% Preallocate results
% w_MC(sim, axis, time)
w_MC = cell(N_MC,1);
time_MC = cell(N_MC,1);

% Save nominal values
J_nom   = J_depl;
w0_nom  = w0;
a_nom   = a;
e_nom   = e;
incl_nom = incl;
raan_nom = raan;
SRP_nom = SRP;

for k = 1:N_MC

    fprintf('Monte Carlo run %d / %d \n', k, N_MC);

    % Random seed
    rng(k);
    base_seed_MC = base_seed;
    base_seed = base_seed_MC + k; 

    % INERTIA (±10% sui principali, ±5% cross)
    J_depl = J_nom .* diag(1 + 0.10*randn(3,1));

    % ANGULAR VELOCITIES ±30%
    w0 = w0_nom .* (1 + 0.30*randn(3,1));

    % ORBITS UNCERTAINTIES
    a     = a_nom    * (1 + 0.005*randn);     % ±0.5%
    e     = max(0, e_nom + 0.001*randn);      % small variation
    incl  = incl_nom + 0.1*randn;             % ±0.1 deg
    raan  = raan_nom + 0.2*randn;             % ±0.2 deg

    kep = [a,e,incl,raan,w,theta];

    % SRP DISTURBANCES UNCERTAINTY
    fields = fieldnames(SRP_nom);
    for ii = 1:length(fields)
        SRP.(fields{ii}).A = SRP_nom.(fields{ii}).A * (1 + 0.05*randn);
        SRP.(fields{ii}).r_F = SRP_nom.(fields{ii}).r_F .* (1 + 0.10*randn(1,3));
    end

    % Simulation
    simout = sim('SIM_SAD_DEMO','FastRestart','off');

    % Saving results
    w_MC{k} = simout.w_B;
    time_MC{k} = simout.tout;    
end

% Creating matrix for angular velocities
Nt = length(time_MC{1});
w_MC_mat = zeros(N_MC,3,Nt);
for k = 1:N_MC
    w_MC_mat(k,:,:) = rad2deg(w_MC{k});
end

%% --------- PLOT : MONTE CARLO ANALYSIS ---------- 

figure('Name','Monte Carlo Angular Velocities');
hold on;
N_MC = size(w_MC_mat,1);
Nt   = size(w_MC_mat,3);
t = time_MC{1};
for k = 1:N_MC
    plot(t, squeeze(w_MC_mat(k,1,:)),'LineWidth',1); % w_x
    plot(t, squeeze(w_MC_mat(k,2,:)),'LineWidth',1); % w_y
    plot(t, squeeze(w_MC_mat(k,3,:)),'LineWidth',1); % w_z
end
xlabel('Time [s]');
ylabel('Angular velocity [deg/s]');
title('Monte Carlo Angular Velocities (w_x, w_y, w_z)');
grid on;
legend({'\omega_x','\omega_y','\omega_z'});
hold off;


% Restoring nominal values
J_depl = J_nom;
w0     = w0_nom;
a      = a_nom;
e      = e_nom;
incl   = incl_nom;
raan   = raan_nom;
SRP    = SRP_nom;



%% ----------- SATELLITE SCENARIO ----------------

%Conversion a km -> m
a = a*1000;         %[m]

%Satellite Scenario 
stopTime = startTime + seconds(2*T);
sampleTime = 60;
sc = satelliteScenario(startTime,stopTime,sampleTime);
viewer = satelliteScenarioViewer(sc,"CameraReferenceFrame","Inertial","Dimension","3D");
sat = satellite(sc,a,e,incl,raan,w,theta);
show(sat)
sat.Visual3DModel = "SmallSat.glb";
coordinateAxes(sat, Scale=2); % red = x_B; green = y_B; blue = z_B
camtarget(viewer, sat);
groundTrack(sat,"LeadTime",3600,"LeadLineColor",[0 1 0],"TrailLineColor",[0 1 0]);
play(sc,PlaybackSpeedMultiplier=500)

%Conversion a m -> km
a = a/1000;  %[km] 

