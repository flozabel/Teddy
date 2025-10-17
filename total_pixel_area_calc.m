% By Florian Zabel, University of Basel, Switzerland (2024)
% Contact: florian.zabel@unibas.ch

clc;
clear all;

resolution=0.5; %
for row=1:360
  latitude(row)=90-(row*resolution)+(resolution/2);
  lat=latitude(row);
  radius_a=6378.1370000;  %GRS80, WGS84 Äquatorradius [km]
  radius_p=6356.7523140;  %GRS80, WGS84 Polradius [km]
  radius=(radius_a*cosd(lat)+radius_p*(1-cosd(lat)));
  pixel_area_km(row)=((2*pi*radius)/(360/resolution)*cosd(lat))*((2*pi*radius)/(360/resolution));
end
pixel_area_km=pixel_area_km';
save ('total_pixel_area_0p5degree.mat','pixel_area_km','-v7.3');