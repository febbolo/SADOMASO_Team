%% Earth-S/C

clear
clc;
close all

% S/C Inertia 
J_depl = [100.9, 0, 0;...
          0, 25.1, 0;...
          0, 0, 91.6]*1e-2; %kg*m^2
% J_depl = [8, 0, 0;...
%           0, 4, 0;...
%           0, 0, 6]*1e-2; %kg*m^2

% Orbital Parameters for LEO orbit (Ex. 3 Lab 2 Orbital Mechanics)
e = 0.1976; % Low eccentricity orbit : for plot/debug/simulation
% e = 0;
% incl = astroConstants(63);   %[deg]
incl = 0;   %[deg]
raan = 0;   %[deg]
w = 0;     %[deg]
theta = 0;  %[deg]
mu = astroConstants(13); %[km^3/s^2]
a = 8350 * 1000; %[m]

% Orbital Data
G = 6.67e-20; %km^3/(kg*s^2)
Mt = 5.97e24;   %kg
Re = astroConstants(23); %km

% Mean motion
n = sqrt((mu/ ((a/1000)^3)));   %[rad/s]

%Orbital period
T = 2*pi/n;     %[s]

% Initial Conditions
w0 = [0; 0; n];  %[rad/s]
% w0 = [1e-6; 1e-6; n];  %[rad/s]

% Satellite Scenario 
% startTime = datetime(2025,3,21,8,0,0);
% stopTime = startTime + seconds(T);
% sampleTime = 60;
% sc = satelliteScenario(startTime,stopTime,sampleTime);
% viewer = satelliteScenarioViewer(sc,"Camera–ReferenceFrame","Inertial","Dimension","3D");
% sat = satellite(sc,a,e,incl,raan,w,theta);
% show(sat)
% groundTrack(sat,"LeadTime",3600,"LeadLineColor",[0 1 0],"TrailLineColor",[0 1 0]);
% play(sc,PlaybackSpeedMultiplier=500)

% Conversion a m -> km
a = a/1000;         %[km]

% Creating initial condition Keplerian elements vector
kep = [a,e,incl,raan,w,theta];

% Plot for debugging (r)
% [s0,Q] = kep2car(a,e,incl,raan,w,theta,mu);
% N = 5000;
% [s] = twobody(s0,mu,linspace(0,5*T,N));
% r = vecnorm(s(:,1:3),2,2);
% figure
% plot(linspace(0,5*T,N),r)
% [E_v,t_v,T] = kepler_timelaw_vec(e,a,mu,0,0,5,N);
% theta_v = 2*atan2d( (sqrt(1+e) * tand(E_v/2) ), sqrt(1-e) );
% r = r(:)';
% r_N = [r .* (cosd(raan).*cosd(w+theta_v) - sind(raan).*sind(w+theta_v).*cosd(incl)); ...
%        r .* (sind(raan).*cosd(w+theta_v) + cosd(raan).*sind(w+theta_v).*cosd(incl)); ...
%        r .* (sind(w+theta_v).*sind(incl))];
% figure
% plot(linspace(0,5*T,N),r_N)
% legend


%% Sun-Earth

% Revolution Solar Year Period
T_sun = 365*24*60*60;        % [s]
% Computing mean motion of Earth around the Sun
n_sun = 2*pi/T_sun;         %[rad/s]
% Obliquity of Equatorial Plane
eps_E = deg2rad(astroConstants(63)); %[rad]

% Astronomic unit
ua = astroConstants(2);   %[km]

% Distance of the Sun
r_sun = 1*ua;       %[km]


%% Solar Radiation Pressure (SRP) Torque

F_e = 1358; %[W/m^2]
F0=237;
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

%% GG torque (for inclined and elliptical orbits)

% Adimensional coefficients for stability
I_moments = diag(J_depl);
K_yaw = (I_moments(3)-I_moments(2))/I_moments(1)
K_roll = (I_moments(3)-I_moments(1))/I_moments(2)
K_pitch = (I_moments(2)-I_moments(1))/I_moments(3)

% Verify Stability conditions
if((K_roll*K_yaw>0) && ( (1+3*K_roll+K_roll*K_yaw)^2 > 16*K_yaw*K_roll) )
    disp('Stability conditions are satisfied.');
else
    disp('Stability conditions are not satisfied.');
end


%% Magnetic Torque

j_B = [0.01; 0.05; 0.01]; % [A*m^2], dipole vector of S/C
w_E = deg2rad(15.04)/3600; %[rad/s], average angular velocity of the Earth 

% Coefficients for H_0, Simple Dipole Model (DGRF 2000)
g0_1 = -29404.8e-9; %[T = N/(A*m)]
g1_1 = -1450.9e-9;  %T
h1_1 = 4652.5e-9;   %T

% Initial conditions for ra_GST
theta_G_0 = 0;  %[deg]
t0 = 0;         %[s]



% %% CHECKING ORTHONORMALITY && ATTITUDE
% 
% simout = sim('sim_lab7_disturbances');
% time = simout.tout; 
% A_B_N = simout.A_B_N;
% A_B_N_ortho = simout.A_B_N_ortho;
% A_B_LVLH = simout.A_B_LVLH;
% w_B_LVLH = simout.w_B_LVLH;
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
% plot(time, w_B_LVLH(:,1), 'b', 'LineWidth', 1.5);
% title('Angular Velocity in X Direction');
% xlabel('Time (s)');
% ylabel('error on w_x (rad/s)');
% grid on;
% subplot(3, 1, 2);
% plot(time, w_B_LVLH(:,2), 'r', 'LineWidth', 1.5);
% title('Angular Velocity in Y Direction');
% xlabel('Time (s)');
% ylabel('error on w_y (rad/s)');
% grid on;
% subplot(3, 1, 3);
% plot(time, w_B_LVLH(:,3), 'g', 'LineWidth', 1.5);
% title('Angular Velocity in Z Direction');
% xlabel('Time (s)');
% ylabel('error on w_z (rad/s)');
% grid on;

%% CHECKING FFT
% Bx_fft = abs(fft(Bx));
% By_fft = abs(fft(By));
% f = (0:length(Bx_fft)-1)/length(Bx_fft);
% 
% figure
% plot(f,Bx_fft,'Marker','o','LineWidth',2)
% hold on
% plot(f, n/(2*pi)*ones(size(f)), 'Marker','x', 'LineWidth', 2, 'DisplayName', 'n/(2*pi)');
% legend('show');
% title('FFT of Bx')
% xlim([0,4]);
% figure
% plot(f,By_fft,'Marker','o','LineWidth',2)
% hold on
% plot(f, n/(2*pi)*ones(size(f)), 'Marker','x', 'LineWidth', 2, 'DisplayName', 'n/(2*pi)');
% legend('show');
% 
% title('FFT of By')
% xlim([0,4]);


%% Per albedo try vectors

I_moments = diag(J_depl);
K_yaw = (I_moments(3)-I_moments(2))/I_moments(1);
K_roll = (I_moments(3)-I_moments(1))/I_moments(2);
K_pitch = (I_moments(2)-I_moments(1))/I_moments(3);

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



