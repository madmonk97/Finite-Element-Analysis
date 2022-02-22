% Given Data

dia     = 0.05;         % Diameter of aluminium rod in meters
E       = 70*10^9;      % Youngs modulus in N/m2
rho     = 2700;         % Density in kg/m^3
m       = 1;            % mass in kg
L       = 1;            % Length of the rod in m
A       = pi*(dia^2)/4; % Cross sectional area of rod in m^2
N       = 2000;         % RPM
omega   = 2*pi*N/60;    % angular velocity in rad/s

% Displacement 'u'


r = 1
u = (rho*omega^2/(2*E))*(2*m*L*r/(rho*A) + L^2*r - r^3/3)

