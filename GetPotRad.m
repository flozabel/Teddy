% By Florian Zabel, University of Basel, Switzerland (2025)
% Contact: florian.zabel@unibas.ch

function potRadiation=GetPotRad(jahr,monat,tag,stunde,minute,sekunde,lat,long,slope,aspect)

solconst=1367; %[W/m2]

time=datenum([num2str(tag),'.',num2str(monat),'.',num2str(jahr),' ',num2str(stunde),':',num2str(minute),':',num2str(sekunde)],'dd.mm.yyyy HH:MM:SS');
doy=day(datetime(time,'ConvertFrom','datenum'),'dayofyear');

if(monat>2)
  jj=jahr;
  mj=monat;
else
  jj=jahr-1;
  mj=monat+12;
end

dj=tag;
hj=stunde/24+minute/1440+sekunde/86400;
aj=floor(jj/100);
bj=2-aj+floor(aj/4);
jd_null=floor(365.25*(jj+4716))+floor(30.6001*(mj+1))+dj+bj-1524.5;
jd=floor(365.25*(jj+4716))+floor(30.6001*(mj+1))+dj+hj+bj-1524.5;
n=jd-2451545.0;

L=280.460*pi/180+0.9856474*pi/180*n; %[rad]
g=357.528*pi/180+0.9856003*pi/180*n; %[rad]
ekliptikale_laenge=L+1.915*pi/180*sin(g)+0.01997*pi/180*sin(2*g); %[rad]
ekliptik=23.439*pi/180-0.0000004*pi/180*n; %[rad]
rektaszension=atan(cos(ekliptik)*sin(ekliptikale_laenge)/cos(ekliptikale_laenge)); %[rad]
rektaszension=rektaszension*180/pi; %[grad]
if(cos(ekliptikale_laenge)<0)
  rektaszension=rektaszension+180;
end
rektaszension=rektaszension*pi/180; %[rad]
deklination=asin(sin(ekliptik)*sin(ekliptikale_laenge)); %[rad]

djul=(jd_null-2451545.0)/36525;
ut=6.697376+2400.05134*djul+1.002738*(stunde+(minute/60));
stundenwinkel_ut=ut*15*pi/180; %[rad] 
stundenwinkel_ort=stundenwinkel_ut+long*pi/180; %[rad] 
stundenwinkel=stundenwinkel_ort-rektaszension; %[rad] 

nenner=cos(stundenwinkel)*sin(lat*pi/180)-tan(deklination)*cos(lat*pi/180); %[rad]
azimut=atan(sin(stundenwinkel)/nenner); %[rad]
azimut=azimut*180/pi; %[grad]
if(nenner<0.)
  azimut=azimut+180;
end

h=asin(cos(deklination)*cos(stundenwinkel)*cos(lat*pi/180)+sin(deklination)*sin(lat*pi/180)); %[rad]
h=h*180/pi; %[grad]

day_ang=2.*pi*doy/365.25; %[rad]
eccentric=1+0.03344*cos(day_ang-0.048869); %[rad]
solconst= solconst *eccentric; %[W/m2]

potRadiation=solconst*sin(h*pi/180); %[W/m2]

hh=asin(sin(h*pi/180)*cos(slope*pi/180)+cos(h*pi/180)*cos(aspect*pi/180-azimut*pi/180)*sin(slope*pi/180)); %[rad]
hh=hh*180/pi; %[grad]

potRadiation=potRadiation*(sin(hh*pi/180)/sin(h*pi/180)); %[W/m2s]
if(potRadiation<0)
  potRadiation=0;
end

end