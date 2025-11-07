function wfde5_all=wfde5_tiles(upperyy,loweryy,leftxx,rightxx,years_wfde5,startyear_wfde5,wfde5dir,parameters_wfde5,xcol_global,yrow_global,parallel_cpus)

%open parallel pool
delete(gcp('nocreate'));
parpool('Processes',min(parallel_cpus,12));

disp(['processing tile ',num2str(upperyy),'° - ',num2str(loweryy),'°']);

%read hourly wfde5 reanalysis data from netcdf file
wfde5_all=zeros(366*24,years_wfde5,length(xcol_global),length(parameters_wfde5),'single');

for varloop=1:length(parameters_wfde5)
  parameter_wfde5=parameters_wfde5{varloop};
  path=[wfde5dir,filesep,parameter_wfde5];

  for yy=1:years_wfde5
    year=yy+startyear_wfde5-1;
    data_all=zeros(744,length(xcol_global),12,'single');
    parfor mm=1:12
      filelist=dir([path,filesep,parameter_wfde5,'*_',num2str(year),num2str(mm,'%.2i'),'*.nc']);
      filename=[filelist.name];
      cdir=[filelist.folder];
      data_all(:,:,mm)=read_wfde5(cdir,filename,xcol_global,yrow_global,parameter_wfde5);
    end%months
    for mm=1:12
      if(leap_year(year)==1)
        dayspermonth = [31 29 31 30 31 30 31 31 30 31 30 31];
      else
        dayspermonth = [31 28 31 30 31 30 31 31 30 31 30 31];
      end
      b=sum(dayspermonth(1:mm))*24;
      a=b-dayspermonth(mm)*24+1;
      wfde5_all(a:b,yy,:,varloop)=data_all(1:dayspermonth(mm)*24,:,mm);
    end
  end%years
end%varloop

if(~exist([wfde5dir,filesep,'wfde5_tiles'],'dir'))
  mkdir([wfde5dir,filesep,'wfde5_tiles'])
end
disp(['saving WFDE5 tile to file ',wfde5dir,filesep,'wfde5_tiles',filesep,'wfde5_all_',num2str(upperyy),'_',num2str(loweryy),'_',num2str(leftxx),'_',num2str(rightxx),'.mat']);
save([wfde5dir,filesep,'wfde5_tiles',filesep,'wfde5_all_',num2str(upperyy),'_',num2str(loweryy),'_',num2str(leftxx),'_',num2str(rightxx),'.mat'],'wfde5_all','-v7.3');

end



function wfde5=read_wfde5(cdir,filename,xcol_global,yrow_global,parameter_wfde5)
finfo = ncinfo([cdir,filesep,filename]);
varname=finfo.Variables(1,4).Name;
x_nc=finfo.Dimensions(1,find(strcmpi({finfo.Dimensions.Name},'lon')==1)).Length;
y_nc=finfo.Dimensions(1,find(strcmpi({finfo.Dimensions.Name},'lat')==1)).Length;
z_nc=finfo.Dimensions(1,find(strcmpi({finfo.Dimensions.Name},'time')==1)).Length;
data=zeros(y_nc,x_nc,z_nc);
data2d=zeros(length(xcol_global),z_nc);
wfde5=zeros(744,length(xcol_global));

disp(['read hourly WFDE5 file ',filename]);
data=ncread([cdir,filesep,filename],varname,[1 1 1],[720 360 z_nc]);
%data=rot90(data,1);
for xloop=1:length(xcol_global)
  x=xcol_global(xloop);
  y=yrow_global(xloop);
  x_rot90=360-y+1;
  y_rot90=x;
  data2d(xloop,:)=squeeze(data(y_rot90,x_rot90,:));
end
wfde5(1:z_nc,:)=permute(data2d(:,1:z_nc),[2 1]);

%add snowfall to rainfall
if strcmpi(parameter_wfde5,'Rainf')==1
  wfde5_snow=zeros(744,length(xcol_global));
  filename=strrep(filename,'Rainf','Snowf');
  cdir=strrep(cdir,'Rainf','Snowf');
  finfo = ncinfo([cdir,filesep,filename]);
  varname=finfo.Variables(1,4).Name;
  disp(['read hourly WFDE5 file ',filename]);
  data=ncread([cdir,filesep,filename],varname,[1 1 1],[720 360 z_nc]);
  %data=rot90(data,1);
  for xloop=1:length(xcol_global)
    x=xcol_global(xloop);
    y=yrow_global(xloop);
    x_rot90=360-y+1;
    y_rot90=x;
    data2d(xloop,:)=data(y_rot90,x_rot90,:);
  end
  wfde5_snow(1:z_nc,:)=permute(data2d(:,1:z_nc),[2 1]);
  wfde5(:,:)=wfde5(:,:)+wfde5_snow(:,:);
end
end
