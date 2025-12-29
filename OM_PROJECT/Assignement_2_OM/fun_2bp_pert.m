% Define the function for the two-body problem
function dydt = fun_2bp_pert(t, s, acc_pert_fun, parameters)
  
% fun_2bp_pert is a function used to describe a perturbed keplerian motion
% for two body with acceleration as input as handle function
% parameters = [Primary.J2, Primary.mu, Primary.mass, Primary.Radius, t0_mjd2000]


% PROTOTYPE
% dy= fun_2bp(t, y, mu);
%
% INPUT
%       t[1]    Time( here can be omitted)
%       y[6x1]  State of the body (rx, ry, rz, vx, vy, vz)
%       acc_pert_fun[3x1] column vector with ECI perturbation
%       parameters   constants
% 
% OUTPUT:
%       dy[6x1] Derivative of the state of the body
% 
% CONTRIBUTORS
%       Luca Deli
%
% VERSIONS
%       3/10/2025

    % Fattore comune
    
r = s(1:3);
rnorm = norm(r);

% Define the function for the two-body problem

v = s(4:6); % velocity vector

% Accelerazione centrale
a_c = -parameters(2) * r / rnorm^3;
a_tot = a_c + acc_pert_fun(t, s, parameters);

% Derivate del sistema
dydt = [v; a_tot];
    
    % Combine the derivatives into a single vector
end