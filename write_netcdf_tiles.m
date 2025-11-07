function write_netcdf_tiles(data,lat,lon,timestamp,var_name,var_name_long,unit,time_unit,filename,comment,total_rows)

filename_nc=[filename,'.nc4'];

rows=size(data,1);
cols=size(data,2);
bands=size(data,3);

%change dimension order and orientatin for netcdf
data=single(data);
data=permute(data,[2 1 3]);

%additional information for header
institution=['University of Basel, Department of Environmental Sciences'];
contact=['Florian Zabel, florian.zabel@unibas.ch'];

creation_time=date;
timesteps=int32(timestamp);

if(~exist(filename_nc,'file')) %if file does not yet exist, write new netcdf file with attributes

  %create variables with dimensions
  nccreate(filename_nc,'lon','Datatype','double','Dimensions',{'lon' cols},'ChunkSize',[cols],'Format','netcdf4','DeflateLevel',9);
  nccreate(filename_nc,'lat','Datatype','double','Dimensions',{'lat' inf},'ChunkSize',[total_rows],'Format','netcdf4','DeflateLevel',9);
  nccreate(filename_nc,'time','Datatype','int32','Dimensions',{'time' bands},'ChunkSize',[bands],'Format','netcdf4','DeflateLevel',9);
  nccreate(filename_nc,var_name,'Datatype','single','Dimensions',{'lon' cols 'lat' inf 'time' bands},'ChunkSize',[cols,total_rows,1],'Format','netcdf4','FillValue',single(1.0E20),'DeflateLevel',9);

  %write attributes
  ncwriteatt(filename_nc,'/','Institution',institution);
  ncwriteatt(filename_nc,'/','Contact',contact);
  ncwriteatt(filename_nc,'/','Comment',comment);
  ncwriteatt(filename_nc,'/','Creation_time',creation_time);

  %write data and data attributes
  ncwrite(filename_nc,'lon',lon); ncwriteatt(filename_nc,'lon','units','degrees_east'); ncwriteatt(filename_nc,'lon','long name','lon');
  ncwrite(filename_nc,'lat',lat); ncwriteatt(filename_nc,'lat','units','degrees_north'); ncwriteatt(filename_nc,'lat','long name','lat');
  ncwrite(filename_nc,'time',timesteps); ncwriteatt(filename_nc,'time','units',time_unit);ncwriteatt(filename_nc,'time','long name','time');ncwriteatt(filename_nc,'time','calendar','standard');
  ncwrite(filename_nc,var_name,data); ncwriteatt(filename_nc,var_name,'units',unit); ncwriteatt(filename_nc,var_name,'long name',var_name_long);

  ncdisp(filename_nc);

elseif(exist(filename_nc,'file')) %if file already exists, append data

  validate_lat=ncread([filename_nc],'lat');
  ncwrite(filename_nc,'lat',lat,numel(validate_lat)+1);
  ncwrite(filename_nc,var_name,data,[1,numel(validate_lat)+1,1]);
  ncdisp(filename_nc);

end

end %function