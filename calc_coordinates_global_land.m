%calc lat/lon for all samples at 0.5° resolution
function [lat,lon,lat_all,lon_all]=calc_coordinates_global_land(uppery,lowery,leftx,rightx,res,mask)

lat_all=uppery:-res:lowery;  %extend y
lat_all=lat_all-(res/2);
lat_all(end)=[];
lon_all=leftx:res:rightx; %extend x
lon_all=lon_all+(res/2);
lon_all(end)=[];

lon_global=repmat(lon_all,1,length(lat_all));

n=1;
for y=1:length(lat_all)
  for x=1:length(lon_all)
    lat_global(1,n)=lat_all(y);
    n=n+1;
  end
end

yrow=floor((90-lat_global)/res)+1;
xcol=floor((180+lon_global)/res)+1;

n=1;
for m=1:length(lat_global)
  if(mask(yrow(m),xcol(m))==0)
    continue
  else
    lat(n)=lat_global(m);
    lon(n)=lon_global(m);
    n=n+1;
  end
end

if ~exist('lat','var') %FZ added 07.10.2025
  lat=[];
  lon=[];
end

return
%save('coordinates_global_land.mat','lat','lon','-v7.3');

end