% Henkle.m
%
% This script solves for the amplitide of 
% sinusoidal pressure change ?p using the relation described Henkle, 1970: 
% The Role of Waves in Causing Submarine Landslides. 
%
% Input Parameters
% 	
% L	Wavelength (m)
% d	Waveheight (m)
% k 	Cu/Gamma'z 	
% Gamma Density (g/ml)
% z	Water Depth (m)
% y	Longitude (m)
% HL	Half Length - Henkle's variable x (m)
%
% The inputs Y,Z are taken as three columns in the array "data()"
%
% Calculated Parameters
%
% Beta	Slope Angle (�) 
% A  	sin(alpha)-(alpha)cos(alpha)
% B  	sin(theta)-(theta)cos(theta)/sin^3(theta)
%
% Output Parameters
%
% P	Amplitude of Sinusoidal Pressure Change


%% Define Constant Values
k=0.04; 
gamma=0.7;
wave_period = 5; %Seconds - Typical of shelf...
deepwater_wave_height = 1; % Meters... Just guessing!


%% Import YZ data from .csv file 

% Calculate slope angle 'Beta' (rate of change of Y relative to Z)
% 1st derivative delta_y / delta_z
%data = importdata('ew_n_2015_new.txt');
data = importdata('2015_central_ns_transect.txt');
%data = data.data;
y = data(:,2);
z = data(:,3);

%% Calculate Wavelength

% Note:  doc file mentions just assuming a 20m wavelength. 
% To do that, skip this and just make L = 20;

% Run wavelength function to calculate L for a given depth 
L = zeros(size(z));
for i = 1:size(z)
    [l,e,c] = wavelength(wave_period, z(i),30*pi/180,0.25,150*pi/180);
    L(i) = l;
end

%% Calculate Waveheight

% Run waveheight function to calculate d for a given depth z and wavelength
d = zeros(size(z));
for i = 1:size(z)
    d(i) = w_height(L(i), z(i), deepwater_wave_height);
end

%% Define x

% Half-length for calculations
% This is the half length of the area that's failing ("x" in Henkel's
% paper)
% The maximum delta-p will occur when this is half of a wavelength. We need
% to calculate the maximum delta-p, so we'll set it to half the surface
% wavelength.
x = L / 2;


%% Calculate Theta

% Now we need to know the depth to which the sliding surface extends. This
% is "d" in Henkel's paper, and is required to calculate "theta".  We'll
% assume that the sliding surface extends at the same dip as the seafloor
% over the entire length we've defined for HL above (to maximize things,
% it's one-half the surface wavelength).
% Unfortunately, Henkel defines things in such a way that it's very
% difficult to calculate "theta" from the equations given.  It's much
% easier to use an alternate formulation calculated from the fact that
% "theta", "x", and "d" describe a circular segment. 
% E.g. see http://en.wikipedia.org/wiki/Circular_segment (Note that "d"
% here is Henkel's d. It's equivalent to "h" in the wikipedia article.)

% Using gradient instead of diff so that the output is the same size as the
% input.
beta = abs(gradient(z) ./ gradient(y));
d = beta .* x; 
d(isnan(d)) = 0; % If y-values are identical, set gradient to 0

% Finishing theta calculation based on circular segment
r = d / 2 + L.^2 ./ (8 * d);
theta = 2 .* asin(L ./ (2 * r));

%% Calculate pressure difference required for failure

% Calculate alpha for values of z for which we have HL and L 
alpha = (2*pi*x) ./ L;

% Calculate A & B
A=sin(alpha) - (alpha .* cos(alpha)); 
B=(sin(theta) - theta .* cos(theta) ) ./ sin(theta).^3; 

% Finally, let's calculate the pressure differential
delta_p = (k * gamma * L * 4 * pi^2) .* (x ./ L).^3 ...
          .* (B - beta / (3 * k)) ./ A;


%% Quick plot

plotyy(y, -z, y, delta_p);
%show()





	