%% -------------- TO DO -----------------
% ONLY SCRIPT / SIMULINK STUFF
% 1) Cycle for B_IGRF and optimization
% 2) (Ordering Albedo and Earth Radiation Pressure)
% 3) Pink noise (Lowpass Filter) NOT DOING MAYBE
% 5) Control
% 6) Actuators

%% ------------ DATA : Earth-S/C - Sun synchronous Orbit ------------

clear
clc
close all

% S/C Inertia 
% J_depl = [100.9, 0, 0;...
%           0, 25.1, 0;...
%           0, 0, 91.6]*1e-2; %kg*m^2
% 6U CubeSat inertia
J_depl = [18, 0, 0;...
          0, 12, 0;...
          0, 0, 22]*1e-2; %kg*m^2

% Constants
G = 6.67e-20; % [km^3/(kg*s^2)]
Mt = 5.97e24;   % [kg]
Re = astroConstants(23); %[km]
J2 = 0.00108263; % Second Zonal harmonic
mu = astroConstants(13); %[km^3/s^2]

% Date Time 
startTime = datetime(2025,1,1,0,0,0);   % Initial time epoch of project (TO BE DEFINED)

% Orbital Parameters
e = 0.005; % Low eccentricity orbit
draan = 1.99096871e-7; % Sun-synchronous rate of precession of raan [rad/s]
zp = 600; %[km]
rp = zp + Re; %[km]
a = rp * (1 - e);  % Semi-major axis [km]
raan = 0;   %[deg]
w = 0;     %[deg]
theta = 0;  %[deg]
incl = rad2deg( acos( -2/3 * (1-e^2)^2*a^(7/2) * draan / ( sqrt(mu)*J2*Re^2 ) )); %[deg]

% Mean motion
n = sqrt( (mu/ (a)^3));   %[rad/s]

%Orbital period
T = 2*pi/n;     %[s]

% Initial Conditions
w0 = [0.5; 0.7; 0.6];  %[rad/s]
% w0 = [1e-6; 1e-6; n];  %[rad/s]

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
Dt = str2double(get_param('SIM_SAD_DEMO','FixedStep'));

q.maxsensor = min(Sun_sensor.Ts, Magmeter.Ts);
q.Ts = min(q.maxsensor, Dt);

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

% DEFINING THRESHOLDS FOR THE DIFFERENT CASES 
Control.w_tumbling = deg2rad(2); % Treshold for tumbling, [rad]
Control.w_pointing = deg2rad(0.5);     % Treshold for pointing, [rad]
% In the middle we have slew manouvre 

w_enter = deg2rad(2);
w_exit = deg2rad(5);
w_enter_2 = deg2rad(0.5);
w_exit_2 = deg2rad(1);

% ------------ DE-TUMBLING ---------------
% B_dot method
% B_dot_m = -w_est x B_m
% D_required = -k_b*B_dot_m
% M_c = D_required x B = (-k_b*B_dot_m) x B_m

% Gain of the B_dot : TUNING  
k_b = 200000;

% -------- SLEW / RE-POINTING -------------

% Bryson Method to define Q,R
% Arbitrary choice of maximum acceptable error : ~0.5 deg
alpha_max = deg2rad(0.5); %[rad]
% Since the slew manouvre starts right after the detumbling we consider as
% omega_max the treshold on the detumbling case
% Defining maximum torque for RW:  TO BE DEFINED!!!
M_RW_max = 3e-3; %[Nm]

x_max = [Control.w_pointing,Control.w_pointing, Control.w_pointing, alpha_max, alpha_max, alpha_max]; 
u_max = [M_RW_max, M_RW_max, M_RW_max];

Q = diag(1./(x_max).^2); 
R = diag(1./(u_max).^2);
   
% Considering a long time of transient of P(t) matrix, P(t)≃ cost ->
% algebraic Riccati equation 
% The control is working at discrete time -> dlqr command is needed
% To avoid K matrix too small (1e-20 order), the information about the
% sampling time should given through the system state sys = ss(A, B, C, D, ts)
Ts = 0.2; % Sampling time fron sensors information
sysc = ss(A, B, [], []); 
sysd = c2d(sysc, Ts, 'zoh'); 
[K,S,P] = dlqr(sysd.A,sysd.B,Q,R);    
% [K,S,P] = lqr(A,B,Q,R); 
 
% A the beginning, there are big angles variations -> an extension of
% linear control is needed, using THE SAME GAIN MATRIX found before

Kp = K(:, 4:6); % Proportional gain matrix dim.(3x3)
Kd = K(:, 1:3); % Derivative gain matrix dim.(3x3)

% The only effectively different part is in the simulink model, where the state vector components 
% are substituted with the correct terms due to non-small angles. 


% %% PID pole placement for single-axis rotational dynamics: I*theta_ddot = T
% % Plant: G(s) = 1/(I*s^2)
% % Controller (PID): C(s) = Kp + Ki/s + Kd*s
% % Closed-loop characteristic (unity feedback):
% %   1 + C(s)G(s) = 0  ->  I s^3 + Kd s^2 + Kp s + Ki = 0
% % Divide by I:
% %   s^3 + (Kd/I)s^2 + (Kp/I)s + (Ki/I) = 0
% %
% % Choose desired poles p1,p2,p3, build desired polynomial:
% %   (s-p1)(s-p2)(s-p3) = s^3 + a2 s^2 + a1 s + a0
% % Match coefficients:
% %   Kd = I*a2,  Kp = I*a1,  Ki = I*a0
% 
% clear; clc; close all;
% 
% Ixx = 0.85163;	Ixy = 0.00000;	Ixz = 0.00000;
% Iyx = 0.00000;	Iyy = 6.82684;	Iyz = 0.00000;
% Izx = 0.00000;	Izy = 0.00000;	Izz = 7.18354;
% 
% % Parameters
% I = Izz;   % [kg*m^2] inertia (edit)
% 
% % Desired poles (edit these)
% % Option A: directly set 3 poles (must be stable: negative real parts)
% p1 = -0.005;
% p2 = -0.0050;
% p3 = -0.0012;
% 
% % Option B (comment Option A and use this): 2nd-order shape + extra pole
% % wn   = 10;        % natural frequency [rad/s]
% % zeta = 0.7;       % damping ratio
% % p3   = -5*wn;     % extra real pole, typically 3-10x faster than dominant
% % p1   = -zeta*wn + 1j*wn*sqrt(1-zeta^2);
% % p2   = -zeta*wn - 1j*wn*sqrt(1-zeta^2);
% 
% % Compute PID gains by coefficient matching
% desired_poly = poly([p1 p2 p3]);  % returns [1 a2 a1 a0]
% % Reminder:
% %   POLY(V), when V is a vector, is a vector whose elements are
% %   the coefficients of the polynomial whose roots are the
% %   elements of V. For vectors, ROOTS and POLY are inverse
% %   functions of each other, up to ordering, scaling, and
% %   roundoff error.
% 
% a2 = desired_poly(2);
% a1 = desired_poly(3);
% a0 = desired_poly(4);
% 
% Kd = I*a2;
% Kp = I*a1;
% Ki = I*a0;
% 
% fprintf('Desired poles: [%s]\n', num2str([p1 p2 p3]));
% fprintf('PID gains:\n  Kp = %.6g\n  Ki = %.6g\n  Kd = %.6g\n', Kp, Ki, Kd);
% 
% % Build transfer functions and check
% s = tf('s');
% G = 1/(I*s^2);
% C = Kp + Ki/s + Kd*s;
% 
% L = C*G;
% Tcl = feedback(L, 1);   % closed-loop theta/ref
% 
% % Verify characteristic polynomial numerically:
% % Denominator of Tcl should be I*s^3 + Kd*s^2 + Kp*s + Ki (up to scaling)
% disp('Closed-loop denominator (Tcl):');
% disp('Denominator of Tcl should be I*s^3 + Kd*s^2 + Kp*s + Ki (up to scaling)')
% Tcl_den = Tcl.Denominator{1};
% disp(Tcl_den);
% 
% % Step response
% figure; step(Tcl);
% grid on;
% title('Closed-loop \theta response to step reference');
% 
% % Control effort transfer function (torque / ref)
% % Torque T = C(s)*(ref - theta)
% % With unity feedback: T/ref = C(s)/(1 + C(s)G(s))
% Tu = minreal( C / (1 + C*G) );



%% ------------ ACTUATORS -----------------

B_z_min = 1e-6;

% Magnetorquers - RWp100 (Blue Canyon Technologies)
mag_perm = 0;
n_wind = 0;
coil_area = 0; 

D_max_XY = 0.4; % Am^2
D_max_Z = 0.5; % Am^2

I_max_XY = 0.12; % A
I_max_Z = 0.18; % A

% Reaction Wheel - RW400 (AAC Clyde Space)
T_max = 0.008; % N*m
h_max = 0.05; % N*m*s
omega_max = 5000 * (2*pi/60); % rad/s (5000rpm)
Iw_est = h_max / omega_max; % kg*m^2











% %% -------------- CHECKING ORTHONORMALITY && ATTITUDE --------
% 
% 
% 
% simout = sim('SIM_SAD_DEMO');
% time = simout.tout; 
% A_B_N = simout.A_B_N;
% A_B_N_ortho = simout.A_B_N_ortho;
% A_B_LVLH = simout.A_B_LVLH;
% w_B_LVLH = simout.w_B_LVLH;
% 
% 
% 
% % Pre-allocate Q, error and norm_error
% Q = zeros(3, 3, length(time));
% error = zeros(3, 3, length(time));
% norm_error = zeros(1, length(time));
% 
% % Testing Orthonormality of A
% for i = 1:length(time)
%     Q(:,:,i) =A_B_N(:,:,i)'*A_B_N(:,:,i);
%     error(:,:,i) = abs(eye(3)-Q(:,:,i));
%     norm_error(i) = norm( Q(:,:,i) - eye(3), 'fro' );
% end
% 
% % Plot orthonormality error of A
% figure('Name','Orthonormality error of A')
% plot(time, norm_error)
% 
% % Testing Orthonormality of A_ortho
% for i = 1:length(time)
%     Q(:,:,i) = A_B_N_ortho(:,:,i)'*A_B_N_ortho(:,:,i);
%     error(:,:,i) = abs(eye(3)-Q(:,:,i));
%     norm_error(i) = norm( Q(:,:,i) - eye(3), 'fro' );
% end
% 
% % Plot orthonormality error of A_ortho
% figure('Name','Orthonormality error of A after orthonormalisation')
% plot(time, norm_error)
% 
% % Plot the attitude error between Body Frame and Uniformly rotating LVLH
% % Frame
% figure('Name','Attitude Matrix A_B_LVLH Components over Time');
% for i = 1:3
%     for j = 1:3
%         subplot(3, 3, (i-1)*3 + j);
%         plot(time, squeeze(A_B_LVLH(i, j, :)), 'LineWidth', 1.5);
%         title(['A_{' num2str(i) num2str(j) '} over Time']);
%         xlabel('Time (s)');
%         ylabel(['A_{' num2str(i) num2str(j) '}']);
%         grid on;
%     end
% end
% 
% % Plot the error of the angular velocities (in body frame) of the absolute
% % angular velocity wrt the LVLH angular velocity
% figure('Name','Error on w (B wrt LVLH) over time')
% subplot(3, 1, 1);
% plot(time, w_B_LVLH(1,:), 'b', 'LineWidth', 1.5);
% title('Angular Velocity in X Direction');
% xlabel('Time (s)');
% ylabel('error on w_x (rad/s)');
% grid on;
% subplot(3, 1, 2);
% plot(time, w_B_LVLH(2,:), 'r', 'LineWidth', 1.5);
% title('Angular Velocity in Y Direction');
% xlabel('Time (s)');
% ylabel('error on w_y (rad/s)');
% grid on;
% subplot(3, 1, 3);
% plot(time, w_B_LVLH(3,:), 'g', 'LineWidth', 1.5);
% title('Angular Velocity in Z Direction');
% xlabel('Time (s)');
% ylabel('error on w_z (rad/s)');
% grid on;

% %% ----------- SATELLITE SCENARIO ----------------
% 
% %Conversion a km -> m
% a = a*1000;         %[m]
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
% play(sc,PlaybackSpeedMultiplier=500)
% 
% %Conversion a m -> km
% a = a/1000;         %[km] 


