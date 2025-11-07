function write_netcdf_dailyval(data,lat,lon,timestamp,var_name,var_name_long,unit,time_unit,filename,comment)

filename_nc=[filename,'.nc4'];

rows=size(data,1);
cols=size(data,2);
years=size(data,3);
doys=size(data,4);

%change dimension order and orientatin for netcdf
data=single(data);
data=permute(data,[2 1 3 4]);

%additional information for header
institution=['University of Basel, Department of Environmental Sciences'];
contact=['Florian Zabel, florian.zabel@unibas.ch'];

creation_time=date;
timesteps=timestamp;

%write netcdf
nccreate(filename_nc,'lon','Datatype','double','Dimensions',{'lon' cols},'ChunkSize',[cols],'Format','netcdf4','DeflateLevel',9);
nccreate(filename_nc,'lat','Datatype','double','Dimensions',{'lat' rows},'ChunkSize',[rows],'Format','netcdf4','DeflateLevel',9);
nccreate(filename_nc,'doy','Datatype','int16','Dimensions',{'doy' doys},'ChunkSize',[doys],'Format','netcdf4','DeflateLevel',9);
nccreate(filename_nc,'year','Datatype','int16','Dimensions',{'year' years},'ChunkSize',[years],'Format','netcdf4','DeflateLevel',9);
nccreate(filename_nc,var_name,'Datatype','single','Dimensions',{'lon' cols 'lat' rows 'year' years 'doy' doys},'ChunkSize',[cols,rows,years,1],'Format','netcdf4','FillValue',single(1.0E20),'DeflateLevel',9);

%write attributes
ncwriteatt(filename_nc,'/','Institution',institution);
ncwriteatt(filename_nc,'/','Contact',contact);
ncwriteatt(filename_nc,'/','Comment',comment);
ncwriteatt(filename_nc,'/','Creation_time',creation_time);

%write data and data attributes
ncwrite(filename_nc,'lon',lon); ncwriteatt(filename_nc,'lon','units','degrees_east'); ncwriteatt(filename_nc,'lon','long name','lon');
ncwrite(filename_nc,'lat',lat); ncwriteatt(filename_nc,'lat','units','degrees_north'); ncwriteatt(filename_nc,'lat','long name','lat');
ncwrite(filename_nc,'doy',int16(timesteps)); ncwriteatt(filename_nc,'doy','units',time_unit);ncwriteatt(filename_nc,'doy','long name','doy');ncwriteatt(filename_nc,'doy','calendar','standard');
ncwrite(filename_nc,'year',int16([0:years-1])); ncwriteatt(filename_nc,'year','units','year');ncwriteatt(filename_nc,'year','long name','year');
ncwrite(filename_nc,var_name,data); ncwriteatt(filename_nc,var_name,'units',unit); ncwriteatt(filename_nc,var_name,'long name',var_name_long);

ncdisp(filename_nc);

end %function