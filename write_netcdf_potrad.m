function write_netcdf_potrad(data,lat,lon,timestamp,var_name,var_name_long,unit,time_unit,filename,comment)

filename_nc=[filename,'.nc4'];

rows=size(data,1);
cols=size(data,2);
doys=size(data,3);
hours=size(data,4);

%change dimension order and orientatin for netcdf
data=permute(data,[2 1 3 4]);

%additional information for header
institution=['University of Basel, Deptartment of Environmental Sciences'];
contact=['Florian Zabel, florian.zabel@unibas.ch'];

creation_time=date;
timesteps=timestamp;

%write netcdf
nccreate(filename_nc,'lon','Datatype','double','Dimensions',{'lon' cols},'ChunkSize',[cols],'Format','netcdf4','DeflateLevel',9);
nccreate(filename_nc,'lat','Datatype','double','Dimensions',{'lat' rows},'ChunkSize',[rows],'Format','netcdf4','DeflateLevel',9);
nccreate(filename_nc,'doy','Datatype','int16','Dimensions',{'doy' doys},'ChunkSize',[doys],'Format','netcdf4','DeflateLevel',9);
nccreate(filename_nc,'hour','Datatype','int16','Dimensions',{'hour' hours},'ChunkSize',[hours],'Format','netcdf4','DeflateLevel',9);
nccreate(filename_nc,var_name,'Datatype','single','Dimensions',{'lon' cols 'lat' rows 'doy' doys 'hour' hours},'ChunkSize',[cols,rows,1,24],'Format','netcdf4','FillValue',single(1.0E20),'DeflateLevel',9);

%write attributes
ncwriteatt(filename_nc,'/','Institution',institution);
ncwriteatt(filename_nc,'/','Contact',contact);
ncwriteatt(filename_nc,'/','Comment',comment);
ncwriteatt(filename_nc,'/','Creation_time',creation_time);

%write data and data attributes
ncwrite(filename_nc,'lon',lon); ncwriteatt(filename_nc,'lon','units','degrees_east'); ncwriteatt(filename_nc,'lon','long name','lon');
ncwrite(filename_nc,'lat',lat); ncwriteatt(filename_nc,'lat','units','degrees_north'); ncwriteatt(filename_nc,'lat','long name','lat');
ncwrite(filename_nc,'doy',int16(timesteps)); ncwriteatt(filename_nc,'doy','units',time_unit);ncwriteatt(filename_nc,'doy','long name','doy');ncwriteatt(filename_nc,'doy','calendar','standard');
ncwrite(filename_nc,'hour',int16([0:hours-1])); ncwriteatt(filename_nc,'hour','units','hh');ncwriteatt(filename_nc,'hour','long name','hour');
ncwrite(filename_nc,var_name,data); ncwriteatt(filename_nc,var_name,'units',unit); ncwriteatt(filename_nc,var_name,'long name',var_name_long);

ncdisp(filename_nc);

end %function