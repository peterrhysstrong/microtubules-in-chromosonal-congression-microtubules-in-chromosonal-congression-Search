function [x,y,z] = sph2cart2(az,elev,r)
%function that works out change from cartesian to spherical coordinates in the standard way 
z = r .* cos(elev);
rsinelev = r .* sin(elev);
x = r.*sin(elev) .* cos(az);
y = r.*sin(elev) .* sin(az);
 