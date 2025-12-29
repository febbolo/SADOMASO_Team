% Define the function for the two-body problem
function dydt = fun_2bp(t, y, Mu_E)
    
% fun_2bp is a function used to describe a keplerian motion for two body

% PROTOTYPE
% dy= fun_2bp(t, y, mu);
%
% INPUT
%       t[1]    Time( here can be omitted)
%       y[6x1]  State of the body (rx, ry, rz, vx, vy, vz)
%       mu[1]   Specific energy of the body
% 
% OUTPUT:
%       dy[6x1] Derivative of the state of the body
% 
% CONTRIBUTORS
%       Luca Deli
%
% VERSIONS
%       3/10/2025


% Define the function for the two-body problem
    r = y(1:3); % position vector
    v = y(4:6); % velocity vector
    
    % Calculate the acceleration due to gravity
    a = -Mu_E * r / norm(r)^3;
    
    % Combine the derivatives into a single vector
    dydt = [v; a];
end

