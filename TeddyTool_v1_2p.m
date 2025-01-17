% Teddy-Tool v1.2p: Temporal Disaggregation of Daily Climate Model Data (parallelized version)
% The Teddy-Tool disaggregates daily climate model data to hourly values for climate impact and attribution studies
% Changes in v1.1p: Parallel processing of samples
% Changes in v1.2p: Optimization, partitioning into tiles, spatial output
% By Florian Zabel, University of Basel, Switzerland (2024)
% Contact: florian.zabel@unibas.ch
% Tested with Matlab2023b

%compile TeddyTool_v1_2p
%mcc -m TeddyTool_v1_2p.m -a WFDE5Factors.m WFDE5Factors_LST.m wfde5_tiles.m leap_year.m local_solar_time.m calc_coordinates_global_land.m write_netcdf_tiles.m write_netcdf_dailyval.m parsave.m


clc;
clear;
warning off;
%dbstop if error
%dbclear all


disp(newline)
disp('------------------------------------------------')
disp('TEDDY')
disp('Temporal Disaggregation of Daily Climate Model Data')
disp('Version 1.2p')
disp(newline)
disp('by Florian Zabel (2024)')
disp('University of Basel')
disp('contact: florian.zabel@unibas.ch')
disp('------------------------------------------------')
disp(newline)
pause(1)


disp('Starting Teddy-Tool v1.2p for temporal empirical disaggregation of daily climate model data');
disp('Reading ini file');
%read variables from TeddyTool.ini
filename=['TeddyTool.ini'];
fid = fopen(filename);
line=fgetl(fid);
lat=str2num(strtok(line,'!'));        %latitude
line=fgetl(fid);
lon=str2num(strtok(line,'!'));        %longitude
line=fgetl(fid);
htimestep=str2num(strtok(line,'!'));  %select hourly time step for temporal disaggregation (1=1-hourly,2=2-hourly,etc.)
line=fgetl(fid);
outdir=strtrim(strtok(line,'!'));     %output directory
line=fgetl(fid);
factordir=strtrim(strtok(line,'!'));  %directory containing preprocessed subdaily factors derived from wfde5 data
line=fgetl(fid);
wfde5dir=strtrim(strtok(line,'!'));   %directory containing hourly WFDE5 reference data
line=fgetl(fid);
climdir=strtrim(strtok(line,'!'));    %directory containing ISIMIP climate model data
line=fgetl(fid);
models=strtrim(strtok(line,'!'));     %select ISIMIP climate model(s): GFDL-ESM4,IPSL-CM6A-LR,MPI-ESM1-2-HR,MRI-ESM2-0,UKESM1-0-LL
models=deblank(models);
models=strsplit(models,',');
line=fgetl(fid);
scenarios=strtrim(strtok(line,'!'));  %select ISIMIP scenario(s): ssp126,ssp585,ssp370,historical,picontrol
scenarios=deblank(scenarios);
scenarios=strsplit(scenarios,',');
line=fgetl(fid);
period=strtrim(strtok(line,'!'));
period=strsplit(period,'-');
startyear=str2num(cell2mat(period(1)));       %select startyear for ISIMIP climate data
endyear=str2num(cell2mat(period(2)));         %select endyear for ISIMIP climate data
line=fgetl(fid);
parameters_isimip=strtrim(strtok(line,'!'));  %select variable(s): tas,hurs,rsds,rlds,ps,sfcwind,pr
parameters_isimip=deblank(parameters_isimip);
parameters_isimip=strsplit(parameters_isimip,',');
line=fgetl(fid);
period=strtrim(strtok(line,'!'));
period=strsplit(period,'-');
startyear_wfde5=str2num(cell2mat(period(1))); %select startyear for hourly wfde5 data that is used to determine the diurnal profile
endyear_wfde5=str2num(cell2mat(period(2)));   %select endyear for hourly wfde5 data that is used to determine the diurnal profile
years_wfde5=endyear_wfde5-startyear_wfde5+1;
line=fgetl(fid);
timewindow=str2num(strtok(line,'!'));         %select time window of +- n days around DOY to search for similar conditions in the past (default=11)
line=fgetl(fid);
lstflag=str2num(strtok(line,'!'));            %use LST (local solar time) (1=default) or UTC (0) for data processing?
line=fgetl(fid);
rainfallflag=str2num(strtok(line,'!'));       %consider precipitation on consecutive days: 1=yes, 0=no; doy window flag (ln 13) must be >= 1
line=fgetl(fid);
precipnanflag=str2num(strtok(line,'!'));      %write NaN values for precipitation (0=no=default) in case that no precipitation event can be found in the historical hourly reference dataset (1=yes: in this case mass and energy are not preserved!)
line=fgetl(fid);
parallel_cpus=str2num(strtok(line,'!'));      %number of parallel workers
line=fgetl(fid);
tilesize=str2num(strtok(line,'!'));           %tile flag and size: value>0:use tiles and set tile size; 0:use coordinates in line 1,2
fclose(fid);

if(floor(24/htimestep)~=24/htimestep)
  error('timestep cannot be divided by 24');
end

if(tilesize>0)
  boundaries_lat=lat;
  boundaries_lon=lon;
end

res=0.5; %spatial resolution
parameters_isimip_read={'tas','hurs','rsds','rlds','ps','sfcwind','pr','tasmax','tasmin'}; %for method climate analogies, all variables are required
parameters_wfde5_factors={'tas','hurs','rsds','rlds','ps','sfcwind','pr'};
parameters_wfde5={'Tair','Qair','SWdown','LWdown','PSurf','Wind','Rainf'};

var_names_long={'air temperature','humidity','showrtwave radiation','longwave radiation','air pressure','wind speed','precipitation'};
units={'K','%','W m-2','W m-2','hPa','ms-1',['mm per ',num2str(htimestep),' hour']};

for scenarioloop=1:length(scenarios)
  scenario=scenarios{scenarioloop};
  disp(['Processing scenario ',scenario]);

  for modelloop=1:length(models)
    model=models{modelloop};
    disp(['Processing model ',model]);

    if(tilesize>0) %flag for tile-mode
      %calc global coordinates for tiles
      uppery=boundaries_lat(1);
      lowery=boundaries_lat(2);
      leftx=boundaries_lon(1);
      rightx=boundaries_lon(2);

      for upperyy=uppery:-tilesize:lowery+tilesize
        loweryy=upperyy-tilesize;

        crosscheck=0;
        aoutdir=[outdir,filesep,model,filesep];

        if(~exist([aoutdir],'dir'))
          mkdir([aoutdir]);
        end

        %read land-sea-mask (WFDE5 Land-Sea-Mask)
        disp('Reading land-sea mask');
        mask=load('wfde5_mask.mat').mask;
        %mask(mask==0)=1; %without mask

        [lat,lon,lat_all,lon_all]=calc_coordinates_global_land(upperyy,loweryy,leftx,rightx,res,mask); %added Florian Zabel, 27.03.2024

        %convert lat/lon to global row,col
        yrow_global=floor((90-lat)/res)+1;
        xcol_global=floor((180+lon)/res)+1;

        %convert lat/lon to row,col for tile
        res=0.5;
        yrow=floor((upperyy-lat)/res)+1;
        xcol=floor((-leftx+lon)/res)+1;

        %check if tile has already been calculated
        crosscheck=zeros(length(parameters_wfde5_factors),1);
        for varloop=1:length(parameters_wfde5_factors)
          varname=parameters_isimip{varloop};
          filename=[aoutdir,varname,'_',model,'_',scenario,'_',num2str(startyear),'-',num2str(endyear),'.nc4'];
          %read lat for nc4 file
          if(exist(filename,'file'))
            validate_lat=ncread([filename],'lat');
            if(sum(validate_lat==lat_all,'all') == numel(lat_all))
              crosscheck(varloop)=1;
            end
          end
        end
        if(sum(crosscheck)==length(parameters_wfde5_factors))
          disp(['tile ',model,'_',num2str(upperyy),'_',num2str(loweryy),'_',num2str(leftx),'_',num2str(rightx),' already exists and is skipped']);
          continue
        end

        %calculate local solar time offset to UTC
        lst_offsets=local_solar_time(res);

        %calculate local midnight in UTC
        midnights=lst_offsets*(-1)+1;
        midnights(midnights<=0)=midnights(midnights<=0)+24;
        numdays=daysact(['1-jan-',num2str(startyear)],['31-dec-',num2str(endyear)])+1;

        %for older Matlab Versions
        % start_date = datetime(['1-Jan-', num2str(startyear)], 'InputFormat', 'dd-MMM-yyyy');
        % end_date = datetime(['31-Dec-', num2str(endyear)], 'InputFormat', 'dd-MMM-yyyy');
        % numdays = days(end_date - start_date) + 1;

        potRadAll=zeros(366,24,length(xcol),'single');
        %read potential radiation
        if(find(contains(parameters_isimip,'rsds','IgnoreCase',true)>=1))
          disp('Reading potential radiation');
          filename=['potential_radiation.nc4'];
          data=ncread([filename],'prad');
          for xloop=1:length(xcol)
            x=xcol_global(xloop);
            y=yrow_global(xloop);
            if(mask(y,x)==0)
              continue
            end
            potRadAll(:,:,xloop)=data(x,y,:,:);
          end
          clear data
        end

        %check if precalculated factors exist
        for varname=parameters_isimip
          varname=char(varname);
          %read precalculated multi-year mean subdaily factors
          if(strcmpi(varname,'pr')==1)
            filename_wfde5_fact=[factordir,filesep,varname,'_WFDE5_sum_doy_',num2str(startyear_wfde5),'-',num2str(endyear_wfde5),'.nc4'];
          else
            filename_wfde5_fact=[factordir,filesep,varname,'_WFDE5_mean_doy_',num2str(startyear_wfde5),'-',num2str(endyear_wfde5),'.nc4'];
          end
          if(~exist(filename_wfde5_fact,'file'))
            prompt = ['Precompiled daily mean/sum does not exist. Compile from WFDE5 data? (y|n)?'];
            CalcWFDE5FactorsFlag = input(prompt,'s');
            if isempty(CalcWFDE5FactorsFlag)
              CalcWFDE5FactorsFlag = 'y';
            end
            if(strcmpi(CalcWFDE5FactorsFlag,'y')==1)
              if(lstflag==1)
                WFDE5Factors_LST(parameters_isimip,wfde5dir,factordir,startyear_wfde5,endyear_wfde5);
              else
                WFDE5Factors(parameters_isimip,wfde5dir,factordir,startyear_wfde5,endyear_wfde5);
              end
            else
              disp('Stopped processing due to missing subdaily factors');
              return
            end
          end
        end

        %read daily ISIMIP climate model data
        vloop=1;
        data_cm=zeros(numdays,length(xcol),length(parameters_isimip_read),'single');
        time_cm=zeros(numdays,1,'double');
        for varname=parameters_isimip_read
          varname=char(varname);

          path_cm=[climdir,filesep,scenario,filesep,model];
          filelist = dir([path_cm,filesep,model,'*',varname,'_global_daily*.nc']);
          nofileflag=0;
          if(isempty(filelist)==1)
            disp(['ERROR: No files in folder ',path_cm]);
            nofileflag=1;
            break
          end
          n=1;
          for fileloop=1:length(filelist)
            %open ISIMIP climate data
            filename_cm=[path_cm,filesep,filelist(fileloop).name];
            ncid = netcdf.open([filename_cm],'NOWRITE');
            finfo = ncinfo([filename_cm]);
            xcm=finfo.Dimensions(1,find(strcmpi({finfo.Dimensions.Name},'lon')==1)).Length;
            ycm=finfo.Dimensions(1,find(strcmpi({finfo.Dimensions.Name},'lat')==1)).Length;
            zcm=finfo.Dimensions(1,find(strcmpi({finfo.Dimensions.Name},'time')==1)).Length;
            yearoffset=filename_cm(length(filename_cm)-11:length(filename_cm)-8);
            if(str2double(yearoffset)>endyear)
              continue
            end
            %Get timeoffset from netcdf file (different climate models use different time offsets)
            time_varid=netcdf.inqVarID(ncid,'time');
            time_units=netcdf.getAtt(ncid,time_varid,'units');
            parts=strsplit(time_units,' ');
            reference_date=parts{3};
            reference_datetime=strcat(reference_date,{' 00:00:00'});
            ref_datetime=datetime(reference_datetime,'InputFormat','yyyy-MM-dd HH:mm:ss');
            timeoffset=datenum(ref_datetime);
            % if (strcmpi(model,'GSWP3-W5E5'))
            %   if (startyear>=1981)
            %     timeoffset=datenum(['01.01.',num2str(1900),' 00:00:00'],'dd.mm.yyyy HH:MM:SS');
            %   else
            %     timeoffset=datenum(['01.01.',num2str(1860),' 00:00:00'],'dd.mm.yyyy HH:MM:SS');
            %   end
            % else
            %   timeoffset=datenum(['01.01.',num2str(yearoffset),' 00:00:00'],'dd.mm.yyyy HH:MM:SS');
            % end
            timeread=ncread([filename_cm],'time');
            timeread=double(timeread);
            if(timeread(1)~=0)
              modelspec_offset=timeread(1);
            else
              modelspec_offset=0;
            end
            timeoffset=timeoffset+(timeread(1)-modelspec_offset);
            time_cm_read_all=timeread+timeoffset;
            time_cm_read=timeread+timeoffset;
            time_cm_read(time_cm_read>datenum(endyear,12,31))=[];
            time_cm_read(time_cm_read<datenum(startyear,1,1))=[];
            index=find(ismember(time_cm_read_all,time_cm_read)==1);
            if(isempty(index))
              disp(['skip ',filename_cm])
              continue
            end
            acm=index(1);
            bcm=index(end);
            zcm=bcm-acm+1;
            time_cm_read_all=timeread+timeoffset;
            time_cm(n:n+zcm-1,1)=time_cm_read_all(acm:bcm);
            disp(['reading ',filename_cm]);
            %read ISIMIP climate data for lat/lon
            data=ncread([filename_cm],varname,[1 1 acm],[720 360 zcm]);
            netcdf.close(ncid)
            for xloop=1:length(xcol)
              x=xcol_global(xloop);
              y=yrow_global(xloop);
              if(mask(y,x)==0)
                continue
              end
              data_cm(n:n+zcm-1,xloop,vloop)=data(x,y,:);
            end
            clear data
            n=n+zcm;
          end %fileloop
          %convert units
          if(strcmpi(varname,'pr')==1)
            data_cm(:,:,vloop)=data_cm(:,:,vloop)*86400; %convert mass flux density to mm/day
          elseif(strcmpi(varname,'ps')==1)
            data_cm(:,:,vloop)=data_cm(:,:,vloop)/100;   %Pa to hPa
          elseif(strcmpi(varname,'clt')==1)
            data_cm(:,:,vloop)=data_cm(:,:,vloop)/100;   %percent to 0-1
          end
          vloop=vloop+1;
        end %varloop

        if(nofileflag==1)
          continue
        end

        %read preprocessed wfde5 daily mean values for all variables over the climate analogy time period (standard=1980-2018)
        wfde5_mean=zeros(years_wfde5,366,length(xcol),length(parameters_isimip_read),'single');
        vloop=1;
        for varname=parameters_isimip_read
          varname=char(varname);
          if(strcmpi(varname,'pr')==1)
            filename_wfde5_mean=[factordir,filesep,varname,'_WFDE5_sum_doy_',num2str(startyear_wfde5),'-',num2str(endyear_wfde5),'.nc4'];
            vname='sum';
          elseif(strcmpi(varname,'tasmax')==1)
            filename_wfde5_mean=[factordir,filesep,'tas','_WFDE5_max_doy_',num2str(startyear_wfde5),'-',num2str(endyear_wfde5),'.nc4'];
            vname='max';
          elseif(strcmpi(varname,'tasmin')==1)
            filename_wfde5_mean=[factordir,filesep,'tas','_WFDE5_min_doy_',num2str(startyear_wfde5),'-',num2str(endyear_wfde5),'.nc4'];
            vname='min';
          else
            filename_wfde5_mean=[factordir,filesep,varname,'_WFDE5_mean_doy_',num2str(startyear_wfde5),'-',num2str(endyear_wfde5),'.nc4'];
            vname='mean';
          end
          disp(['reading ',filename_wfde5_mean])
          ncid = netcdf.open([filename_wfde5_mean],'NOWRITE');
          data=ncread([filename_wfde5_mean],vname,[1 1 1 1],[720 300 years_wfde5 366]);
          netcdf.close(ncid);
          for xloop=1:length(xcol)
            x=xcol_global(xloop);
            y=yrow_global(xloop);
            if(mask(y,x)==0)
              continue
            end
            wfde5_mean(:,:,xloop,vloop)=data(x,y,:,:);
          end
          clear data
          vloop=vloop+1;
        end

        filename=['bestrank',filesep,'bestrank','_',model,'_',scenario,'_',num2str(startyear),'-',num2str(endyear),'_',num2str(upperyy),'_',num2str(loweryy),'_',num2str(leftx),'_',num2str(rightx),'.mat'];
        if(exist(filename,'file'))
          load(filename);
        else
          %open parallel pool
          delete(gcp('nocreate'));
          parpool('Threads',parallel_cpus);
          %parpool('Threads'); %for older Matlab versions

          bestrank=zeros(numdays,2,length(xcol),'int16');
          ttimesteps=size(data_cm,1);

          for xloop=1:length(xcol)
            if(mask(yrow_global(xloop),xcol_global(xloop))==0)
              continue
            end
            disp(['Processing sample ',num2str(xloop),' of ',num2str(length(xcol)),' to find most similar day']);

            parfor t=1:ttimesteps
              time=time_cm(t);
              doy=day(datetime(time,'ConvertFrom','datenum'),'dayofyear');
              year=str2double(datestr(time,'yyyy'));
              if(doy>=60 && leap_year(year)==0) %doy 60 is skipped if no leap year
                doy=doy+1;
              end

              %find most similar day in past (best fit/fingerprint)
              abserr=NaN(years_wfde5,timewindow*2+1,length(parameters_isimip_read),'single');
              abserr_reshaped=NaN(years_wfde5*(timewindow*2+1),length(parameters_isimip_read),'single');
              timewindowsteps_all = zeros(years_wfde5, timewindow*2+1,'int16');
              %valuesstd = NaN(years_wfde5, timewindow*2+1, 'single');
              %stdthreshold = NaN(length(parameters_isimip_read), 'single');

              for vloop=1:length(parameters_isimip_read)
                wfde5_years=startyear_wfde5:endyear_wfde5;
                for ywfde5=1:years_wfde5
                  timewindowsteps=(doy-timewindow:doy+timewindow);
                  if(leap_year(wfde5_years(ywfde5))==0) %wfde5 leap day is at position 60 (29.02)
                    timewindowsteps(timewindowsteps>=60)=timewindowsteps(timewindowsteps>=60)+1;
                  end
                  timewindowsteps(timewindowsteps<1)=timewindowsteps(timewindowsteps<1)+366;
                  timewindowsteps(timewindowsteps>366)=timewindowsteps(timewindowsteps>366)-366;
                  timewindowsteps_all(ywfde5,:)=timewindowsteps;
                  abserr(ywfde5,:,vloop)=abs(squeeze(wfde5_mean(ywfde5,timewindowsteps,xloop,vloop))-data_cm(t,xloop,vloop)');
                  %valuesstd(ywfde5,:)=squeeze(wfde5_mean(ywfde5,timewindowsteps,xloop,vloop));
                end
                %stdthreshold(vloop)=std(valuesstd,0,'all'); %get standard deviation

                %consider rainfall classes, added 22.10.2022
                pr_clim=NaN(1,3);
                pr_ref=NaN(1,3);
                consecutive_pr_states=NaN(years_wfde5,3,timewindow*2+1);
                if(strcmpi(parameters_isimip_read(vloop),'pr')==1)
                  abserror_precip=abserr(:,:,vloop);
                  if(t==1)
                    pr_clim(1,1)=data_cm(t,xloop,vloop);       %climate model precipitation
                    pr_clim(1,2:3)=data_cm(1:t+1,xloop,vloop); %climate model precipitation
                  elseif(t==size(data_cm,1))
                    pr_clim(1,1:2)=data_cm(t-1:t,xloop,vloop); %climate model precipitation
                    pr_clim(1,3)=data_cm(t,xloop,vloop);       %climate model precipitation
                  else
                    pr_clim(1,:)=data_cm(t-1:t+1,xloop,vloop); %climate model precipitation
                  end
                  pr_clim(pr_clim<1)=0; %reduce drizzle effect for precipitation in clim
                  pr_clim(pr_clim>0)=1; %convert clim to boolean
                  timewindowsteps=(doy-timewindow:doy+timewindow);
                  for mw=1:length(timewindowsteps)    %moving window +- 1 day over timewindow
                    for ywfde5=1:years_wfde5
                      mwsteps=[timewindowsteps(mw)-1,timewindowsteps(mw),timewindowsteps(mw)+1];
                      if(leap_year(wfde5_years(ywfde5))==0)
                        mwsteps(mwsteps>=60)=mwsteps(mwsteps>=60)+1;
                      end
                      mwsteps(mwsteps<1)=mwsteps(mwsteps<1)+366;
                      mwsteps(mwsteps>366)=mwsteps(mwsteps>366)-366;
                      pr_ref(ywfde5,:)=wfde5_mean(ywfde5,mwsteps,xloop,vloop); %reference wfde5 precipitation
                    end
                    pr_ref(pr_ref<1)=0;   %reduce drizzle effect for precipitation in ref
                    pr_ref(pr_ref>0)=1;   %convert ref to boolean
                    consecutive_pr_states(:,:,mw)=pr_ref(:,:)==pr_clim(1,:); %check if consecutive days with or without precipitation is the same in both datasets
                  end
                  %same day must have the same precipitation status
                  proof_dayspecific_pr_state=squeeze(consecutive_pr_states(:,2,:));
                  if(sum(proof_dayspecific_pr_state,'all')>=1)
                    abserror_precip(proof_dayspecific_pr_state==0)=NaN; %mark all days that don't have the same precipitation status on same day
                  end
                  if(rainfallflag==1)
                    %consecutive days (yesterday, today, tomorrow) must have same precipitation status
                    proof_consecutive_pr_states=squeeze(sum(consecutive_pr_states,2));
                    if(sum(proof_consecutive_pr_states==3,'all')>=1)
                      abserror_precip(proof_consecutive_pr_states<3)=NaN;  %mark uneven sequence of precipitation days
                    elseif(sum(proof_consecutive_pr_states==2,'all')>=1)
                      abserror_precip(proof_consecutive_pr_states<2)=NaN;  %if three consecutive days do not exist, take at least two secutive days
                    end
                  end
                  abserr(:,:,vloop)=abserror_precip;
                end
              end
              abserr_reshaped(:,:)=reshape(abserr,(timewindow*2+1)*years_wfde5,length(parameters_isimip_read));
              yearn=repmat((1:years_wfde5)',timewindow*2+1,1);
              doyn=reshape(timewindowsteps_all,[(timewindow*2+1)*years_wfde5,1]);

              %remove nan values for sorting
              masknan=sum(isnan(abserr_reshaped),2);
              test=sum(masknan==0);
              if(test>=1) %only if at least one sample remains after exclusion
                %remove samples from population that do not fit the criterias (e.g. precipitation or standard deviation)
                yearn(masknan>0)=[];
                doyn(masknan>0)=[];
                abserr_reshaped(masknan>0,:)=[];
              end
              [srt,index]=sort(abserr_reshaped(:,:),1);
              rank=zeros(length(yearn),length(parameters_isimip_read),'single');
              for vloop=1:length(parameters_isimip_read)
                pos=index(:,vloop);
                idxrepeat=[false diff(srt(:,vloop))'==0];
                rnkNoSkip=cumsum(~idxrepeat);
                rnk=1:numel(pos);
                rnk(idxrepeat)=rnkNoSkip(idxrepeat);
                rank(pos,vloop)=rnk;
              end
              rank(isnan(abserr_reshaped)==1)=NaN; %in case no matching precipitation day is found, ignore precipitation for finding the most similar day
              ranksum=sum(rank,2,'omitnan');
              [mn, idx]=min(ranksum);
              bestrank(t,:,xloop)=[yearn(idx) doyn(idx)]; %year of best fit, doy of best fit
            end
          end
          if(~exist('bestrank','dir'))
            mkdir('bestrank')
          end
          save(filename,'bestrank','-v7.3');
        end
        clear wfde5_mean

        %read hourly wfde5 reanalysis data from netcdf file
        wfde5_all=zeros(366*24,years_wfde5,length(xcol),length(parameters_wfde5),'single');
        %if exist mat file for tile then load
        filename=[wfde5dir,filesep,'wfde5_tiles',filesep,'wfde5_all_',num2str(upperyy),'_',num2str(loweryy),'_',num2str(leftx),'_',num2str(rightx),'.mat'];
        if(exist(filename,'file'))
          disp(['loading WFDE5 tile from file ',filename]);
          load(filename);
        else %call function to calculate the preprocessed wfde5 tiles
          wfde5_all=wfde5_tiles(upperyy,loweryy,leftx,rightx,years_wfde5,startyear_wfde5,wfde5dir,parallel_cpus,mask);
        end

        %convert wfde5 units
        for xloop=1:length(xcol)
          wfde5_all(:,:,xloop,5)=wfde5_all(:,:,xloop,5)/100; %convert surface pressure [Pa] to [hPa]
          wfde5_all(:,:,xloop,7)=wfde5_all(:,:,xloop,7)*3600; %convert mass flux density to [mm h-1]
          %convert specific humidity to relative humidity [%]
          wfde5_t=wfde5_all(:,:,xloop,1)-273.16;%[K] to [°C]
          %es according to Buck equation; t in °C; es in hPa
          a=wfde5_t;
          b=wfde5_t;
          c=wfde5_t;
          d=wfde5_t;
          %over liquid water, T>0°C
          a(wfde5_t>=0)=6.1121;
          b(wfde5_t>=0)=18.678;
          c(wfde5_t>=0)=234.5;
          d(wfde5_t>=0)=257.14;
          %over ice water, T<0°C
          a(wfde5_t<0)=6.1115;
          b(wfde5_t<0)=23.036;
          c(wfde5_t<0)=333.7;
          d(wfde5_t<0)=279.82;
          es_buck=a.*exp((b-(wfde5_t./c)).*(wfde5_t./(d+wfde5_t))); %saturation vapour pressure [hPa]
          e = wfde5_all(:,:,xloop,2) .* wfde5_all(:,:,xloop,5) ./ (0.378 .* wfde5_all(:,:,xloop,2) + 0.622);  %vapour pressure [hPa]
          rh = e ./ es_buck *100; %relative humidity [%]
          rh(rh>100)=100;
          rh(rh<0)=0;
          wfde5_all(:,:,xloop,2)=rh;
        end
        clear wfde5_t rh a b c d es_buck e

        for vloop=1:length(parameters_wfde5_factors)
          varname=parameters_wfde5_factors{vloop};

          if(crosscheck(vloop)==1)
            disp(['skipping ',varname,', since it already exists'])
            continue
          end

          timestepsh=zeros(length(time_cm)*24,1,'double');
          for tloop=1:length(time_cm)
            aa=tloop*24-23;
            bb=aa+23;
            time=time_cm(tloop);
            timestepsh(aa:bb,1)=time:1/24:time+1-(1/24);
          end

          output_temp=zeros(length(xcol),length(timestepsh),'single');
          parfor xloop=1:length(xcol)

            if(mask(yrow_global(xloop),xcol_global(xloop))==0)
              continue
            end

            disp(['Assign hourly climate for ', varname,', processing sample ',num2str(xloop),' of ',num2str(length(xcol))]);

            lst_offset=lst_offsets(xcol_global(xloop));
            if(lstflag==1)
              midnight=24; %set to 11 pm
            else
              midnight=midnights(xcol_global(xloop));
            end

            wfde5=wfde5_all(:,:,xloop,vloop);
            valh=zeros(length(timestepsh),1,'single');

            for tloop=1:length(time_cm)
              aa=tloop*24-23;
              bb=aa+23;
              %disp(['processing subdaily climate ',num2str(tloop),' of ',num2str(length(time_cm)),' for ',varloop])
              time=time_cm(tloop);
              doy=day(datetime(time,'ConvertFrom','datenum'),'dayofyear');
              year=str2double(datestr(time,'yyyy'));
              yy=bestrank(tloop,1,xloop);
              year_wfde5=yy+startyear_wfde5-1;
              doy_wfde5=bestrank(tloop,2,xloop);
              if(leap_year(year_wfde5)==0 && doy_wfde5>=60)  %leap day is at position 366 (31.12)
                doy_wfde5=doy_wfde5-1;
              end

              if(lstflag==1) %in case of using local solar time
                a=(doy_wfde5-1)*24+1-lst_offset;
                b=a+24-1;
                if(leap_year(year_wfde5)==1)
                  leap_offset=0;
                else
                  leap_offset=24;
                end
                if(leap_year(year_wfde5-1)==1)
                  leap_offset_m=0;
                else
                  leap_offset_m=24;
                end
                if(a<1)
                  wfde5_day=wfde5(1:b,yy);
                  if(yy-1>=1)
                    wfde5_day=[wfde5(size(wfde5,1)-leap_offset_m+a:size(wfde5,1)-leap_offset_m,yy-1);wfde5_day;];
                  else %if next year does not exist, take values from next year
                    wfde5_day=[wfde5(size(wfde5,1)-leap_offset+a:size(wfde5,1)-leap_offset,yy);wfde5_day;];
                  end
                elseif(b>size(wfde5,1)-leap_offset)
                  wfde5_day=wfde5(a:size(wfde5,1)-leap_offset,yy);
                  if(yy+1<=size(wfde5,2))
                    wfde5_day=[wfde5_day;wfde5(1:24-length(wfde5_day),yy+1)];
                  else %if next year does not exist, take values from previous year
                    wfde5_day=[wfde5_day;wfde5(1:24-length(wfde5_day),yy)];
                  end
                else
                  wfde5_day=wfde5(a:b,yy); %hourly data of best fit
                end
                %adjust potential radiation to local solar time
                a=1-lst_offset;
                b=24-lst_offset;
                if(a<1)
                  potRad=squeeze(potRadAll(max(doy-1,1),24+a:24,xloop))';
                  potRad=[potRad;potRadAll(doy,1:b,xloop)'];
                elseif(b>24)
                  potRad=squeeze(potRadAll(doy,a:24,xloop))';
                  potRad=[potRad;potRadAll(min(doy+1,366),1:b-24,xloop)'];
                else
                  potRad=squeeze(potRadAll(doy,1:24,xloop))';
                end
              else %otherwise process data in utc
                a=(doy_wfde5-1)*24+1;
                b=a+24-1;
                wfde5_day=wfde5(a:b,yy); %hourly data of best fit
                potRad=squeeze(potRadAll(doy,1:24,xloop))'; %set potential radiation for doy
              end

              if(strcmpi(varname,'pr'))
                factor=wfde5_day./sum(wfde5_day);
                if(sum(wfde5_day)==0 && data_cm(tloop,xloop,vloop)>1)
                  %if no precipitation in wfde5 data, distribute precipitation > 1mm/day
                  %this can be the case under very rare environments with no historical precipitation in the doy-window, e.g. in deserts
                  %expand doy window to +-50 days around doy
                  doysequence_precip=(doy-50)*24:(doy+50)*24;
                  doysequence_precip(doysequence_precip>8784)=doysequence_precip(doysequence_precip>8784)-8784;
                  doysequence_precip(doysequence_precip<1)=doysequence_precip(doysequence_precip<1)+8784;
                  if(sum(wfde5(doysequence_precip,:),'all')>0)
                    wfde5_pr=wfde5(doysequence_precip,:);
                    wfde5_pr=wfde5_pr(:);
                  else
                    wfde5_pr=wfde5(:);
                  end
                  if(sum(wfde5_pr)>0) %if historical precipitation exists
                    %linear regression between precipitation amount and duration
                    wfde5_pr_ds=zeros(length(1:24:length(wfde5_pr)-24),1);
                    wfde5_pr_dh=zeros(length(1:24:length(wfde5_pr)-24),1);
                    n=1;
                    for m=1:24:length(wfde5_pr)-24
                      wfde5_pr_ds(n,1)=sum(wfde5_pr(m:m+24-1));     %precipitation sum per day
                      wfde5_pr_dh(n,1)=sum((wfde5_pr(m:m+24-1)>0)); %hours of precipitation per day
                      n=n+1;
                    end
                    %f = fittype('poly1');
                    %[myfit gof] = fit(wfde5_pr_ds,wfde5_pr_dh,f);
                    %coefficients=coeffvalues(myfit);
                    %plot(myfit,wfde5_pr_ds,wfde5_pr_dh);
                    coefficients=polyfit(wfde5_pr_ds, wfde5_pr_dh,1); %changed 19.09.2024
                    precip_duration=min(24,round(coefficients(1)*data_cm(tloop,xloop,vloop)+coefficients(2)));
                    nighttime=find(potRad==0);
                    if(precip_duration<length(nighttime))
                      random_nighttime=sort(randsample(nighttime,precip_duration));  %randomly select hours
                    else
                      random_nighttime=sort(randsample(1:24,precip_duration));
                    end
                    intensities=rand(length(random_nighttime),1);
                    intensities=intensities/sum(intensities);
                    wfde5_day(random_nighttime)=intensities;
                    factor=wfde5_day./sum(wfde5_day);
                  else %if no historical precipitation exists
                    if (precipnanflag==1)
                      factor=wfde5_day./sum(wfde5_day);  %if precipitation nan flag == 1, write NaN: note that mass- and energy in not preserved in this case!
                    else
                      nighttime=find(potRad==0);
                      precip_duration=randsample(1:length(nighttime),1); %randomly select duration
                      random_nighttime=sort(randsample(nighttime,precip_duration));  %randomly select hours
                      intensities=rand(length(random_nighttime),1);
                      intensities=intensities/sum(intensities);
                      wfde5_day(random_nighttime)=intensities;
                      factor=wfde5_day./sum(wfde5_day);
                    end
                  end
                elseif(sum(wfde5_day)==0)
                  %if precipitation nan flag == 1, write NaN: note that mass- and energy in not preserved in this case!
                  %if no precipitation in wfde5 data, distribute precipitation <1mm/day randomly to hour without radiation.
                  %if no hour without radiation occurs, distribute precipitation to midnight
                  %this should only be the case for drizzling in climate model data vs. no precipitation in wfde5
                  nighttime=find(potRad==0); %all hours without sunlight
                  if(isempty(nighttime))
                    wfde5_day(midnight)=1; %set factor to 1 at local midnight
                  else
                    random_nighttime=randsample(nighttime,1);  %randomly select one hour
                    wfde5_day(random_nighttime)=1;
                  end
                  factor=wfde5_day./sum(wfde5_day);
                end
                wfde5_pr=[];
                nighttime=[];
                random_nighttime=[];
                intensities=[];
                %clear wfde5_pr nighttime random_nighttime intensities
              else
                factor=wfde5_day./mean(wfde5_day);
              end

              %check for NaN-values
              if(strcmpi(varname,'pr')==1 && precipnanflag==1)
                if(tloop==1)
                  disp('crosscheck skipped for precipitation due to NaN flag'); %skip
                end
              elseif(strcmpi(varname,'rsds')==1)
                if(tloop==1)
                  disp('crosscheck skipped for radiation'); %skip
                end
              else
                if (sum(isnan(factor))>0)
                  disp(['crosscheck failed, NaN values occured for ',varname,' on ', datestr(time,'dd.mm.yyyy')]);
                  %disp('execution paused, click return to continue');
                  %pause
                end
              end

              %apply hourly fraction
              valh(aa:bb,1)=data_cm(tloop,xloop,vloop).*factor;

              %for temperature: scale between min and max
              if(strcmpi(varname,'tas')==1)
                data_max=data_cm(tloop,xloop,length(parameters_isimip_read)-1);
                data_min=data_cm(tloop,xloop,length(parameters_isimip_read));
                data=valh(aa:bb,1);
                valh(aa:bb,1)=((data-min(data))/(max(data)-min(data)))*(data_max-data_min)+data_min;
              end

              %for temperature: correct diurnal data to conserve climate model mean while maintaining min and max
              if(strcmpi(varname,'tas')==1)
                deltamean=data_cm(tloop,xloop,vloop)-mean(valh(aa:bb));
                dist=(min(valh(aa:bb)-data_min,abs(valh(aa:bb)-data_max))./((data_max-data_min)/2));
                valhcor=min(valh(aa:bb)+(dist./sum(dist).*deltamean*24),data_max);
                n=0;
                while (round(data_cm(tloop,xloop,vloop)-mean(valhcor),4)>0) %round to 4 decimals
                  deltamean=data_cm(tloop,xloop,vloop)-mean(valhcor);
                  dist=(min(valhcor-data_min,abs(valhcor-data_max))./((data_max-data_min)/2));
                  valhcor=min(valhcor+(dist./sum(dist).*deltamean*24),data_max);
                  n=n+1;
                  if(n>=10)
                    break
                  end
                end
                valh(aa:bb)=valhcor;
              end

              %for relative humidity and cloud cover: limit to max value of 100 and 1 respectively and correct mean value correspondigly
              if(strcmpi(varname,'clt')==1 || strcmpi(varname,'hurs')==1)
                if(strcmpi(varname,'clt')==1)
                  maxthreshold=1;
                elseif(strcmpi(varname,'hurs')==1)
                  maxthreshold=100;
                end
                n=0;
                if(max(valh(aa:bb))>maxthreshold)
                  valhcor=min(valh(aa:bb),maxthreshold);
                  while(round(data_cm(tloop,xloop,vloop)-mean(valhcor),4)>0) %round to 4 decimals
                    deltamean=data_cm(tloop,xloop,vloop)-mean(valhcor);
                    valhcor=min(valhcor+(deltamean*24./sum(valhcor<maxthreshold)),maxthreshold);
                    n=n+1;
                    if(n>=10)
                      break
                    end
                  end
                  valh(aa:bb)=valhcor;
                end
              end

              %for solar radiation, check if values exceed potential maximum radiation for each hour
              if(strcmpi(varname,'rsds')==1)
                valhor=valh(aa:bb);
                valhcor=valh(aa:bb);
                valhcor=min(valhcor,potRad);
                valhcor(valhcor<120)=valhor(valhcor<120); %WMO defined sunshine duration as direct solar irradiance > 120 W/m2. Equivalent to level of solar irradiance shortly after sunrise or before sunset in cloud-free conditions
                valhcor(valhcor<0)=0;
                n=0;
                while(round(data_cm(tloop,xloop,vloop)-mean(valhcor),4)>0.1)
                  deltamean=data_cm(tloop,xloop,vloop)-mean(valhcor);
                  t=find(valhcor>0);
                  valhcor(t)=max(valhcor(t)+deltamean,0);
                  valhcor=min(valhcor,potRad);
                  n=n+1;
                  if(n>=10)
                    break
                  end
                end
                valh(aa:bb)=valhcor;
              end
            end

            valh(valh<0)=0; %set negative values to 0

            %check for NaN values in the dataset
            if(strcmpi(varname,'pr')==1 && precipnanflag==1)
              %skip
            elseif(strcmpi(varname,'rsds')==1)
              %set NaN values for high latitudes during winter without radiation to 0
              valh(isnan(valh)==1)=0;
            else
              if(sum(isnan(valh)==1)>0)
                test = find(isnan(valh));
                teststr=datestr(double(timestepsh(test(:))));
                if(size(teststr,1)>1)
                  disp(['NaN occured in dataset ',varname,' at timestep ',teststr(1,:), ' - ',teststr(end,:)]);
                else
                  disp(['NaN occured in dataset ',varname,' at timestep ',teststr(1,:)]);
                end
                %disp('execution paused, click return to continue');
                %pause
              end
            end

            output_temp(xloop,:)=valh;
          end %xloop

          %aggregate hourly to n-hourly
          if(htimestep>1)
            output_temp_agg=zeros(length(xcol),size(output_temp,2)/htimestep,'single');
            for xloop=1:length(xcol)
              valh=output_temp(xloop,:);
              val=zeros(length(valh)/htimestep,1,'single');
              timesteps=zeros(length(valh)/htimestep,1,'double');
              counter=1;
              for nstep=1:htimestep:length(valh)
                if(strcmpi(varname,'pr')==1)
                  val(counter,1)=sum(valh(nstep:nstep+htimestep-1)); %sum for precipitation
                else
                  val(counter,1)=mean(valh(nstep:nstep+htimestep-1)); %mean
                end
                timesteps(counter,1)=timestepsh(nstep);
                counter=counter+1;
              end
              output_temp_agg(xloop,:)=val;
            end
            output_temp=output_temp_agg;
            clear output_temp_agg
          else
            timesteps=timestepsh;
          end

          %save data as NetCDF file and append missing tiles to existing NetCDF file, if file already exists
          output_spatial=NaN(length(lat_all),length(lon_all),length(timesteps),'single');
          for xloop=1:length(xcol)
            output_spatial(yrow(xloop),xcol(xloop),:)=output_temp(xloop,:);
          end
          disp(['writing NetCDF file for ',varname,', tile ',num2str(upperyy),'-',num2str(loweryy)]);
          timestamp_n(:,1)=0:htimestep:size(output_spatial,3)*htimestep-1;
          var_name_long=var_names_long{vloop};
          unit=units{vloop};
          time_unit=['hours since ',num2str(startyear),'-01-01, 00:00:00'];
          filename=[aoutdir,varname,'_',model,'_',scenario,'_',num2str(startyear),'-',num2str(endyear)];
          comment='Temporal disaggregation by using the Teddy Tool v1.2p, Zabel and Poschlod (2023): Temporal disaggregation of daily climate model data for climate impact analysis. https://doi.org/10.5194/gmd-16-5383-2023';
          total_rows=(uppery-lowery)/res;
          write_netcdf_tiles(output_spatial,lat_all,lon_all,timestamp_n,varname,var_name_long,unit,time_unit,filename,comment,total_rows);

          clear output_spatial output_temp timesteps timestepsh

        end %variable
      end %tiles
    end %if tile-mode
  end %model
end %scenario
