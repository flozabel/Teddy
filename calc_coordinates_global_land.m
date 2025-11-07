%calc lat/lon for all samples at a specific resolution
function [lat,lon,yrow,xcol,lat_,lon_,yrow_global,xcol_global,yrow_global_,xcol_global_]=calc_coordinates_global_land(upperyy,loweryy,leftx,rightx,resolution,mask_global)

%[lat,lon,yrow,xcol,lat_,lon_,yrow_global,xcol_global,yrow_global_,xcol_global_]

lat_=upperyy:-resolution:loweryy;
lat_(end)=[];
lon_=leftx:resolution:rightx;
lon_(end)=[];

%convert lat/lon to row,col in a global file without mask
yrow_global_=floor((90-lat_)/resolution)+1;
xcol_global_=floor((180+lon_)/resolution)+1;

%cut global mask according to selected tile
mask=mask_global(yrow_global_,xcol_global_);

lon_temp=repmat(lon_,1,length(lat_));
n=1;
m=1;
for y=1:length(lat_)
  for x=1:length(lon_)
    if(mask(y,x)==0)
      m=m+1;
      continue
    else
      lat(1,n)=lat_(y);
      lon(1,n)=lon_temp(m);
      n=n+1;
      m=m+1;
    end
  end
end

if ~exist('lat','var') %FZ added 07.10.2025
  lat=[];
  lon=[];
end

%convert lat/lon to row,col for tile
yrow=floor((upperyy-lat)/resolution)+1;
xcol=floor((-leftx+lon)/resolution)+1;

%convert lat/lon to row,col in a global file with mask (1d)
yrow_global=floor((90-lat)/resolution)+1;
xcol_global=floor((180+lon)/resolution)+1;

%convert to central coordinates
lat=lat-resolution/2;
lon=lon+resolution/2;
lat_=lat_-resolution/2;
lon_=lon_+resolution/2;

return
%save('coordinates_global_land.mat','lat','lon','-v7.3');

end