% By Florian Zabel, University of Basel, Switzerland (2025)
% Contact: florian.zabel@unibas.ch

clc
clear

%calculate aspect and slope from terrain
%read wfde5 terrain
filename=['..',filesep,'WFDE5_v2p1',filesep,'dgm',filesep,'ASurf_WFDE5_CRU_v2.1.nc'];
varname='ASurf';
dgm=rot90(ncread(filename,varname,[1 1],[720 360]));

mask=dgm;
mask(isnan(mask)==0)=1;
mask(isnan(mask)==1)=0;

resolution=0.5;

lons=-179.75:resolution:179.75;
lats=89.75:-resolution:-89.75;

cols=length(lons);
rows=length(lats);

load('total_pixel_area_0p5degree.mat');
dgm(isnan(dgm)==1)=0;
for y=1:rows
  for x=1:cols
    lat=lats(y);
    lon=lons(x);
    xl=x;
    yl=y;
    if(x==1)
      xl=x+1;
    end
    if(x==cols)
      xl=x-1;
    end
    if(y==1)
      yl=y+1;
    end
    if(y==rows)
      yl=y-1;
    end
    dgm_matrix(1,1)=dgm(yl-1,xl-1);
    dgm_matrix(1,2)=dgm(yl-1,xl);
    dgm_matrix(1,3)=dgm(yl-1,xl+1);
    dgm_matrix(2,1)=dgm(yl,xl-1);
    dgm_matrix(2,2)=dgm(yl,xl);
    dgm_matrix(2,3)=dgm(yl,xl+1);
    dgm_matrix(3,1)=dgm(yl+1,xl-1);
    dgm_matrix(3,2)=dgm(yl+1,xl);
    dgm_matrix(3,3)=dgm(yl+1,xl+1);
    delta_x=(dgm_matrix(1,3)+2*dgm_matrix(2,3)+dgm_matrix(3,3)-dgm_matrix(1,1)-2*dgm_matrix(2,1)-dgm_matrix(3,1))/8;    %delta höhe in west-ost-richtung
    delta_y=(dgm_matrix(1,1)+2*dgm_matrix(1,2)+dgm_matrix(1,3)-dgm_matrix(3,1)-2*dgm_matrix(3,2)-dgm_matrix(3,3))/8;    %delta höhe in nord-süd-richtung
    pixel_length=sqrt(pixel_area_km(y));
    pixel_length=pixel_length*1000; %in meter
    %slope in % (0:1)
    slope(y,x)=sqrt(delta_x^2+delta_y^2)/pixel_length;  %nach pythagoras die höhe aus verschiedenen himmelsrichtungen bestimmen -> kartesisches system (a²+b²=c² -> a=sqrt(c²+b²)). resultierende höhe / länge (1000m) = steigung
    %slope in ° (0:90)
    slope(y,x)=180.0/pi*atan(slope(y,x));  %atan gibt das ergebnis in bogenmaß (rad) aus -> Grad = Bogenmaß * 180/pi (Bogenmaß = Grad * pi/180)
    aspect(y,x)=180.0/pi*atan(delta_x/delta_y);
    if(delta_y>0.)
      aspect(y,x)=aspect(y,x)+180;
    end
    if(delta_y<0. && delta_x>0.)
      aspect(y,x)=aspect(y,x)+360;
    end
    %Exposition ändern in 0°=Süd, 90°=West,180°=Nord,270°=Ost)
    if(aspect(y,x)>=180)
      aspect(y,x)=aspect(y,x)-180;
    else
      aspect(y,x)=aspect(y,x)+180;
    end
  end
end

year=2000;
months=[31 29 31 30 31 30 31 31 30 31 30 31];
potRadiation_1=zeros(rows,cols,'single');
potRadiation_2=zeros(rows,cols,'single');
potRadiation_3=zeros(rows,cols,'single');
potRadAll=zeros(rows,cols,366,24,'single');

doy=1;
for monat=1:12
  for tag=1:months(monat)
    h=1;
    disp(['DOY: ',num2str(doy)]);
    for stunde=0:23
      parfor x=1:cols
        lons=-179.75:0.5:179.75;
        lats=89.75:-0.5:-89.75;        
        for y=1:rows
          if(mask(y,x)==0)
            continue
          end          
          %disp(['compiling potential radiation for ',num2str(tag),'.',num2str(monat),'.',num2str(jahr)]);          
          potRadiation_1(y,x)=GetPotRad(year,monat,tag,stunde,20,0,lats(y),lons(x),slope(y,x),aspect(y,x));
          potRadiation_2(y,x)=GetPotRad(year,monat,tag,stunde,40,0,lats(y),lons(x),slope(y,x),aspect(y,x));
          potRadiation_3(y,x)=GetPotRad(year,monat,tag,stunde,60,0,lats(y),lons(x),slope(y,x),aspect(y,x));
        end
      end      
      potRadAll(:,:,doy,h)=(potRadiation_1+potRadiation_2+potRadiation_3)/3;
      h=h+1;
    end
    doy=doy+1;    
  end
end

%save as mat file
save('potrad.mat','potRadAll','-v7.3');

%write as netcdf file
load('lat');
load('lon');
timestamps=1:366;
var_name='prad';
var_name_long='potential hourly incoming shortwave radiation';
unit_ncdf=['W m-2'];
time_unit=['day of year (doy) [1-366], with leap day at doy 60'];
comment=['Calculated by Florian Zabel, University of Basel, Switzerland'];
filename=['potential_radiation'];
write_netcdf_potrad(potRadAll,lat,lon,timestamps,var_name,var_name_long,unit_ncdf,time_unit,filename,comment); %save as netcdf
