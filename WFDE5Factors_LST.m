% By Florian Zabel, University of Basel, Switzerland (2024)
% Contact: florian.zabel@unibas.ch

function WFDE5Factors_LST(parameters_out,path_wfde5,outdir,startyear,endyear,lat_all_global,lon_all_global,res,global_mask,parallel_cpus)

units_out=[];
parameters_wfde5=[];

delete(gcp('nocreate'));
parpool('Processes',parallel_cpus);

if(find(contains(parameters_out,'tas','IgnoreCase',true)>=1))
  units_out=[units_out,{'K'}];
  parameters_wfde5=[parameters_wfde5,{'Tair'}];
end
if(find(contains(parameters_out,'hurs','IgnoreCase',true)>=1))
  units_out=[units_out,{'%'}];
  parameters_wfde5=[parameters_wfde5,{'Qair'}];
end
if(find(contains(parameters_out,'rsds','IgnoreCase',true)>=1))
  units_out=[units_out,{'W m-2'}];
  parameters_wfde5=[parameters_wfde5,{'SWdown'}];
end
if(find(contains(parameters_out,'rlds','IgnoreCase',true)>=1))
  units_out=[units_out,{'W m-2'}];
  parameters_wfde5=[parameters_wfde5,{'LWdown'}];
end
if(find(contains(parameters_out,'ps','IgnoreCase',true)>=1))
  units_out=[units_out,{'hPa'}];
  parameters_wfde5=[parameters_wfde5,{'PSurf'}];
end
if(find(contains(parameters_out,'sfcwind','IgnoreCase',true)>=1))
  units_out=[units_out,{'m s-1'}];
  parameters_wfde5=[parameters_wfde5,{'Wind'}];
end
if(find(contains(parameters_out,'pr','IgnoreCase',true)>=1))
  units_out=[units_out,{'mm/h'}];
  parameters_wfde5=[parameters_wfde5,{'Rainf'}];
end

nyears=endyear-startyear+1;

yy=size(global_mask,1);
xx=size(global_mask,2);

%local solar time
lst_offset=local_solar_time(res);
lst_offset=repmat(lst_offset,[yy,1]);

%calculate local midnight in utc
midnight=lst_offset(1:300,:)*(-1)+1;
midnight(midnight<=0)=midnight(midnight<=0)+24;
midnight_lst=ones(yy,xx); %set 0am in local solar time

for parameterloop=1:length(parameters_out)
  parameter_wfde5=parameters_wfde5{parameterloop};
  parameter_out=parameters_out{parameterloop};
  unit_out=units_out{parameterloop};
  
  factor=NaN(yy,xx,24,366);
  amp_prctle=NaN(yy,xx,4,366,'single');
  dailyval=NaN(yy,xx,nyears,366,'single');
  if(strcmpi(parameter_out,'tas')==1)
    dailyval_min=NaN(yy,xx,nyears,366,'single');
    dailyval_max=NaN(yy,xx,nyears,366,'single');
  end
  
  for doy=1:366
    disp(['DOY: ',num2str(doy)]);
    data_wfde5=NaN(yy,xx,24,nyears);
    if(strcmpi(parameter_out,'pr')==1)
      data_wfde5_snowf=NaN(yy,xx,24);
    end
    if(strcmpi(parameter_out,'hurs')==1)
      data_wfde5_t=NaN(yy,xx,24,nyears);
      data_wfde5_p=NaN(yy,xx,24,nyears);
    end
    
    n=1;
    for year=startyear:endyear
      
      Month=month(doy);  %month
      DayMonth=day(doy); %day of month
      band=DayMonth*24-23;
      if(leap_year(year)==0 && Month==2 && DayMonth==29)
        n=n+1;
        continue
      end
      
      data=NaN(yy,xx,24);
      data_compound=NaN(yy,xx,48);
      
      cdir=[path_wfde5,filesep,parameter_wfde5,filesep];
      filelist = dir([cdir,'*_',num2str(year),num2str(Month,'%.2i'),'*.nc']);
      filename=[filelist.name];
      filename=filename(1:length(filename)-3);
      
      disp(['read WFDE5 file ',filename,' for DOY ',num2str(doy)]);
      finfo = ncinfo([cdir,filename,'.nc']);
      varname=finfo.Variables(1,4).Name;
      x=finfo.Dimensions(1,find(strcmpi({finfo.Dimensions.Name},'lon')==1)).Length;
      y=finfo.Dimensions(1,find(strcmpi({finfo.Dimensions.Name},'lat')==1)).Length;
      z=finfo.Dimensions(1,find(strcmpi({finfo.Dimensions.Name},'time')==1)).Length;
      vartime='time';
      times = ncread([cdir,filename,'.nc'],vartime);
      times=times/24;
      timeoffset=datenum(['01.01.',num2str(1900),' 00:00:00'],'dd.mm.yyyy HH:MM:SS');
      formatOut='dd.mm.yyyy HH:MM:SS';
      times=double(times)+timeoffset;
      year_nc=str2double(datestr(times(band),'yyyy'));
      month_nc=str2double(datestr(times(band),'mm'));
      day_nc=str2double(datestr(times(band),'dd'));
      
      %read netcdf file for doy+-12h
      previousmonth=Month-1;
      if(previousmonth<1)
        previousmonth=12;
        previousyear=year-1;
      else
        previousyear=year;
      end
      nextmonth=Month+1;
      if(nextmonth>12)
        nextmonth=1;
        nextyear=year+1;
      else
        nextyear=year;
      end
      filelist_previous = dir([cdir,'*_',num2str(previousyear),num2str(previousmonth,'%.2i'),'*.nc']);
      if(isempty(filelist_previous))
        filelist_previous = dir([cdir,'*_',num2str(year),num2str(previousmonth,'%.2i'),'*.nc']);
      end
      filename_p=[filelist_previous.name];
      filelist_next = dir([cdir,'*_',num2str(nextyear),num2str(nextmonth,'%.2i'),'*.nc']);
      if(isempty(filelist_next))
        filelist_next = dir([cdir,'*_',num2str(year),num2str(nextmonth,'%.2i'),'*.nc']);
      end
      filename_n=[filelist_next.name];
      finfo_p = ncinfo([cdir,filename_p]);
      finfo_n = ncinfo([cdir,filename_n]);
      z_p=finfo_p.Dimensions(1,find(strcmpi({finfo.Dimensions.Name},'time')==1)).Length;
      z_n=finfo_n.Dimensions(1,find(strcmpi({finfo.Dimensions.Name},'time')==1)).Length;
      
      if(band==1)
        data_compound(:,:,1:12)=rot90(ncread([cdir,filename_p],varname,[1 1 z_p-12+1],[x y 12]));
      else
        data_compound(:,:,1:12)=rot90(ncread([cdir,filename,'.nc'],varname,[1 1 band-12],[x y 12]));
      end
      data_compound(:,:,13:36)=rot90(ncread([cdir,filename,'.nc'],varname,[1 1 band],[x y 24]));
      if(band+24>z)
        data_compound(:,:,37:48)=rot90(ncread([cdir,filename_n],varname,[1 1 1],[x y 12]));
      else
        data_compound(:,:,37:48)=rot90(ncread([cdir,filename,'.nc'],varname,[1 1 band+24],[x y 12]));
      end
      for y=1:yy
        for x=1:xx
          for t=1:24
            data(y,x,t)=data_compound(y,x,t+12-lst_offset(y,x));
          end
        end
      end
      data_wfde5(:,:,:,n)=data;
      
      if(strcmpi(parameter_out,'pr')==1)
        filename=strrep(filename,'Rainf','Snowf');
        filename_p=strrep(filename_p,'Rainf','Snowf');
        filename_n=strrep(filename_n,'Rainf','Snowf');
        cdir=strrep(cdir,'Rainf','Snowf');
        disp(['read WFDE5 file ',filename,' for DOY ',num2str(doy)]);
        finfo = ncinfo([cdir,filename,'.nc']);
        varname=finfo.Variables(1,4).Name;
        if(band==1)
          data_compound(:,:,1:12)=rot90(ncread([cdir,filename_p],varname,[1 1 z_p-12+1],[x y 12]));
        else
          data_compound(:,:,1:12)=rot90(ncread([cdir,filename,'.nc'],varname,[1 1 band-12],[x y 12]));
        end
        data_compound(:,:,13:36)=rot90(ncread([cdir,filename,'.nc'],varname,[1 1 band],[x y 24]));
        if(band+24>z)
          data_compound(:,:,37:48)=rot90(ncread([cdir,filename_n],varname,[1 1 1],[x y 12]));
        else
          data_compound(:,:,37:48)=rot90(ncread([cdir,filename,'.nc'],varname,[1 1 band+24],[x y 12]));
        end
        for y=1:yy
          for x=1:xx
            for t=1:24
              data(y,x,t)=data_compound(y,x,t+12-lst_offset(y,x));
            end
          end
        end
        data_wfde5(:,:,:,n)=data_wfde5(:,:,:,n)+data;
      end
      
      if(strcmpi(parameter_out,'hurs')==1)
        filename=strrep(filename,'Qair','Tair');
        filename_p=strrep(filename_p,'Qair','Tair');
        filename_n=strrep(filename_n,'Qair','Tair');
        cdir=strrep(cdir,'Qair','Tair');
        disp(['read WFDE5 file ',filename,' for DOY ',num2str(doy)]);
        finfo = ncinfo([cdir,filename,'.nc']);
        varname=finfo.Variables(1,4).Name;
        if(band==1)
          data_compound(:,:,1:12)=rot90(ncread([cdir,filename_p],varname,[1 1 z_p-12+1],[x y 12]));
        else
          data_compound(:,:,1:12)=rot90(ncread([cdir,filename,'.nc'],varname,[1 1 band-12],[x y 12]));
        end
        data_compound(:,:,13:36)=rot90(ncread([cdir,filename,'.nc'],varname,[1 1 band],[x y 24]));
        if(band+24>z)
          data_compound(:,:,37:48)=rot90(ncread([cdir,filename_n],varname,[1 1 1],[x y 12]));
        else
          data_compound(:,:,37:48)=rot90(ncread([cdir,filename,'.nc'],varname,[1 1 band+24],[x y 12]));
        end
        for y=1:yy
          for x=1:xx
            for t=1:24
              data(y,x,t)=data_compound(y,x,t+12-lst_offset(y,x));
            end
          end
        end
        data_wfde5_t(:,:,:,n)=data;
        
        filename=strrep(filename,'Tair','PSurf');
        filename_p=strrep(filename_p,'Tair','PSurf');
        filename_n=strrep(filename_n,'Tair','PSurf');
        cdir=strrep(cdir,'Tair','PSurf');
        disp(['read WFDE5 file ',filename,' for DOY ',num2str(doy)]);
        finfo = ncinfo([cdir,filename,'.nc']);
        varname=finfo.Variables(1,4).Name;
        if(band==1)
          data_compound(:,:,1:12)=rot90(ncread([cdir,filename_p],varname,[1 1 z_p-12+1],[x y 12]));
        else
          data_compound(:,:,1:12)=rot90(ncread([cdir,filename,'.nc'],varname,[1 1 band-12],[x y 12]));
        end
        data_compound(:,:,13:36)=rot90(ncread([cdir,filename,'.nc'],varname,[1 1 band],[x y 24]));
        if(band+24>z)
          data_compound(:,:,37:48)=rot90(ncread([cdir,filename_n],varname,[1 1 1],[x y 12]));
        else
          data_compound(:,:,37:48)=rot90(ncread([cdir,filename,'.nc'],varname,[1 1 band+24],[x y 12]));
        end
        for y=1:yy
          for x=1:xx
            for t=1:24
              data(y,x,t)=data_compound(y,x,t+12-lst_offset(y,x));
            end
          end
        end
        data_wfde5_p(:,:,:,n)=data;
      end
      
      n=n+1;
    end %year
    
    %convert units
    if strcmpi(parameter_out,'ps')==1
      data_wfde5=data_wfde5/100;   %convert surface pressure [Pa] to [hPa]
    elseif strcmpi(parameter_out,'pr')==1
      data_wfde5=data_wfde5*3600;  %convert mass flux density to mm/h
    elseif(strcmpi(parameter_out,'hurs')==1)
      %convert specific humidity to relative humidity [%]
      data_wfde5_t=single(data_wfde5_t-273.16); %[K] to [°C]
      data_wfde5_p=single(data_wfde5_p/100); %[Pa] to [hPa=mbar]
      %es according to Buck equation; t in °C; es in hPa
      a=data_wfde5_t;
      b=data_wfde5_t;
      c=data_wfde5_t;
      d=data_wfde5_t;
      %over liquid water, T>0°C
      a(data_wfde5_t>=0)=6.1121;
      b(data_wfde5_t>=0)=18.678;
      c(data_wfde5_t>=0)=234.5;
      d(data_wfde5_t>=0)=257.14;
      %over ice water, T<0°C
      a(data_wfde5_t<0)=6.1115;
      b(data_wfde5_t<0)=23.036;
      c(data_wfde5_t<0)=333.7;
      d(data_wfde5_t<0)=279.82;
      es_buck=a.*exp((b-(data_wfde5_t./c)).*(data_wfde5_t./(d+data_wfde5_t))); %saturation vapour pressure [hPa]
      e = data_wfde5 .* data_wfde5_p ./ (0.378 .* data_wfde5 + 0.622);  %vapour pressure [hPa]
      rh = e ./ es_buck *100; %relative humidity [%]
      rh(rh>100)=100;
      rh(rh<0)=0;
      data_wfde5=rh;
    end
    
    %reduce noice for solar radiation
    if(strcmpi(parameter_out,'rsds')==1)
      data_wfde5(data_wfde5<0.1)=0;
    end
    
    %calculate daily mean/sum values for climate analogue fingerprint method
    if(strcmpi(parameter_out,'pr')==1)
      dailyval(:,:,:,doy)=squeeze(sum(data_wfde5(:,:,:,:),3,'omitnan'));
    elseif(strcmpi(parameter_out,'tas')==1)
      dailyval(:,:,:,doy)=squeeze(mean(data_wfde5(:,:,:,:),3,'omitnan'));
      dailyval_min(:,:,:,doy)=squeeze(min(data_wfde5(:,:,:,:),[],3,'omitnan'));
      dailyval_max(:,:,:,doy)=squeeze(max(data_wfde5(:,:,:,:),[],3,'omitnan'));
    else
      dailyval(:,:,:,doy)=squeeze(mean(data_wfde5(:,:,:,:),3,'omitnan'));
    end
    
    clear('data_wfde5','data','data_compound','data_wfde5_snowf','data_wfde5_t','data_wfde5_p');
    
  end %doy
  
  if(~exist(outdir))
    mkdir(outdir)
  end
  
  %write daily sum/mean
  load('lat');
  lat_all_global(301:360)=[];
  load('lon');
  timestamps=1:366;
  var_name='mean';
  var_name_long='daily mean value';
  unit_ncdf=[unit_out];
  time_unit=['day of year (doy) [1-366], with leap day at doy 60'];
  comment=[num2str(startyear),'-',num2str(endyear)];
  filename=[outdir,filesep,parameter_out,'_WFDE5_mean_doy_',num2str(startyear),'-',num2str(endyear)];
  if(strcmpi(parameter_out,'pr')==1)
    var_name='sum';
    var_name_long='daily sum value';
    unit_ncdf=[unit_out];
    filename=[outdir,filesep,parameter_out,'_WFDE5_sum_doy_',num2str(startyear),'-',num2str(endyear)];
  end
  write_netcdf_dailyval(dailyval,lat_all_global,lon_all_global,timestamps,var_name,var_name_long,unit_ncdf,time_unit,filename,comment); %save as netcdf
  clear('dailyval');
  
  %write daily max/min for tas
  if(strcmpi(parameter_out,'tas')==1)
    var_name='min';
    var_name_long='daily min value';
    unit_ncdf=[unit_out];
    filename=[outdir,filesep,parameter_out,'_WFDE5_min_doy_',num2str(startyear),'-',num2str(endyear)];
    write_netcdf_dailyval(dailyval_min,lat_all_global,lon_all_global,timestamps,var_name,var_name_long,unit_ncdf,time_unit,filename,comment); %save as netcdf
    var_name='max';
    var_name_long='daily max value';
    unit_ncdf=[unit_out];
    filename=[outdir,filesep,parameter_out,'_WFDE5_max_doy_',num2str(startyear),'-',num2str(endyear)];
    write_netcdf_dailyval(dailyval_max,lat_all_global,lon_all_global,timestamps,var_name,var_name_long,unit_ncdf,time_unit,filename,comment); %save as netcdf
    clear('dailyval_min','dailyval_max');
  end
  
end %parameter

end
