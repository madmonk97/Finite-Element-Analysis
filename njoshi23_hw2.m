function [conn, K] = njoshi23_hw2(element_lengths, element_areas) 

% Specify default values for testing
% rho = 8050 in kg/m^3
E = 200000000000; % in Pa

if nargin == 0 
    element_lengths = [10;20;30;40];   % (in mm) 
    element_areas   = pi*[2;4;2;4].^2; % (in mm^2) 
end 

% Conversion of length to m and area to m^2 respectively 
element_lengths = element_lengths/(10^3);
element_areas   = element_areas/(10^6);

% Code goes here to construct connectivity and stiffness matrices.  
ne   = length(element_lengths);

conn = [1:1:ne;2:1:ne+1];
K    = zeros(ne+1);


for i = 1:ne
    c = conn(:,i);
    Ke = E*element_areas(i)/element_lengths(i)*[1 -1; -1 1];
    K(c,c) = K(c,c) + Ke
end
