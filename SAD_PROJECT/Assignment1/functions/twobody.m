function [s] = twobody(s0,mu,t_v)
%TWOBODY Solve the two-body problem given the initial state vector
%
% INPUTS: s0 : initial state vector, written in the geocentric equatorial frame 
%               columns: [rx ry rz rdot_x rdot_y rdot_z], rows: one for each time
%               step
%         mu: planetary constant [km^3/s^2]
%         t_v : time vector for the integration
%
% OUTPUTS :  s : solution of the two body problem, written in the
%                geocentric equatorial frame 
%               columns: [rx ry rz rdot_x rdot_y rdot_z], rows: one for each time
%               step
% 
% 
% AUTHORS :     Amura Fabio
%

% Choosing tolerances 
options = odeset( 'RelTol', 1e-13,...
                  'AbsTol',1e-14);

% Solving the two-body problem
f = @(t,s) [s(4:6); -mu*s(1:3)/(norm(s(1:3)))^3]; %s(4:6) velocities rdot_x,rdot_y,rdot_z
[~,s] = ode113(f,t_v,s0,options);

end