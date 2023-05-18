function [ H ] = w_height( T,d,H0 )
%function [ H ] = w_height( T,d,H0 )
%   Calculates wave height based on period, water depth, and deep water
%   wave height, using dispersion relationship for wavelenght in function
%  "wavelength"
cd=9.81*T/(2*pi);%deep water phase velocity, p. 153 in Kinsman
cgd=cd/2; %deep water group velocity, p. 153 in Kinsman
L=wavelength(T,d);
c=L/T;
k=2*pi./L;
cg=(c/2).*(1+((2.*k.*d)./sinh(2.*k.*d))); %intermediate water group velocity, 3.4:22.2 in Kinsman, p. 153
H=H0.*(cgd./cg).^.5; %intermediate water wave height for refraction coefficient=1, 3.5:1.3 in Kinsman p. 153

end

