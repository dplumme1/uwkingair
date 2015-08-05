% Example provided by Dr. David Plummer
%
% Here’s my example Matlab code that reads in netcdf KA data 
% (either 1- or 25-Hz standard flight data, as well as the derived 
% cloud physics data). 

% I put the main user-changeable parameters in the first five lines. 
% The input filenames and directories are in lines 7 and 19, with the
% relevant variables following. Everything that comes after that is
% automated until the figure name and plotting parameters, which 
% start in line 90. 

% It should work with any standard files given the correct flight date and time.

sps = 1; % Choice of 1 Hz or 25 Hz variables 
flight = '0802'; % Reads from file with this date in filename
year = '2013';
tstart = '141300'; % start and end of plot, HHMMSS
tend='141700';

    f=netcdf.open(['./',flight,'/',year,flight,'.c',num2str(sps),'.nc'],'NC_NOWRITE'); % file is ./MMDD/yearMMDD.cX.nc, X is either 1 or 25 HZ

	    varid = netcdf.inqVarID(f,'time'); % general time variable
	    time_sec = netcdf.getVar(f,varid);
    
	    varid = netcdf.inqVarID(f,'trose'); % air temperature
	    trose = netcdf.getVar(f,varid);
	    trose(trose<-32000)=nan;
    
    netcdf.close(f)


    f=netcdf.open(['./',flight,'/',year,flight,'_2d.c1.nc'],'NC_NOWRITE'); % similar name but for cloud physics probe obs

	    varid = netcdf.inqVarID(f,'time'); % time variable for 2D data
	    time_sec_2d = netcdf.getVar(f,varid);
	    
	    varid = netcdf.inqVarID(f,'mass2_2dp_OBL'); % mass derived from 2D-P
	    mass2_2dp_OBL = netcdf.getVar(f,varid);
	    mass2_2dp_OBL(mass2_2dp_OBL<-32000)=nan;

	    varid = netcdf.inqVarID(f,'mass2_cip_IBR'); % mass derived from CIP
	    mass2_cip_IBR = netcdf.getVar(f,varid);
	    mass2_cip_IBR(mass2_cip_IBR<-32000)=nan;

    netcdf.close(f)

% Convert times to Matlab format for plotting

    timevec = nan(length(time_sec),1); 
    for i = 1:length(timevec)
        timevec(i) = datenum([2013 01 01 00 00 time_sec(i)]);
    end
   
    timevec_2d = nan(length(time_sec_2d));
    for i = 1:length(timevec_2d)
        timevec_2d(i) = datenum([2013 01 01 00 00 time_sec_2d(i)]);
    end



%%%%%%%%%%%%%%%%%%% Find the indices to use when plotting the selected time period (ind1:ind2)

start_time = datenum([year,'-',flight(1:2),'-',flight(3:4),' ',tstart(1:2),':',tstart(3:4),':',tstart(5:6)]);
end_time = datenum([year,'-',flight(1:2),'-',flight(3:4),' ',tend(1:2),':',tend(3:4),':',tend(5:6)]);

timediff = abs(timevec-start_time);
ind1 = min(find(timediff==min(timediff)));
timediff = abs(timevec-end_time);
ind2 = min(find(timediff==min(timediff)));




%%%%%%%%%%%%%%%%%%%% 25-Hz data are ordered in a 2D array, 1-Hz are 1D vectors

if(sps<2) % 1-Hz data
  ptime = timevec(ind1:ind2);

% variables that can be 1 Hz or 25 Hz
  ptemp = trose(ind1:ind2);

else % 25-Hz data
  ptime = nan(length(ind1:ind2)*sps,1);
  ptemp = nan(length(ind1:ind2)*sps,1);

  indnew = 1;
  for indorig = ind1:ind2
    ptime(indnew:indnew+sps-1) = timevec(indorig):((1/sps)/86400):timevec(indorig)+((1-(1/sps))/86400); % interpolate time variable to 25-Hz increments
    ptemp(indnew:indnew+sps-1) = trose(:,indorig);
    indnew = indnew+sps;
  end
end

% Variables that are always 1-Hz
  ptime1s = timevec(ind1:ind2);
  pmass2_cip = mass2_cip_IBR(ind1:ind2);
  pmass2_2dp = mass2_2dp_OBL(ind1:ind2);



%%%%%%%%%%%%%%%%%%%% Example figure, air temperature & CIP/2D-P mass

    figname = ['./',flight,'/',year,flight,'_',datestr(start_time,'HHMMSS'),'_',datestr(end_time,'HHMMSS'),'_',num2str(sps),'Hz.png'];

    h = figure('Position',[100 100 1300 1000]);
    set(h,'Renderer','painters')
    subplot(211)
    plot(ptime,ptemp,'linewidth',2)
    xlim([timevec(ind1) timevec(ind2)])
    set(gca,'fontsize',15)
    datetick('x','keeplimits')
    title([flight,': ',num2str(sps),'-Hz Air temperature (C)'])

    subplot(212)
    plot(ptime1s,pmass2_cip,'color','r','linewidth',2)
    xlim([timevec(ind1) timevec(ind2)])
    hold on
    plot(ptime1s,pmass2_2dp,'color','b','linewidth',2)
    set(gca,'fontsize',15)
    datetick('x','keeplimits')
    title([flight,': 1-Hz Mass2 CIP (red), 2DP (blue), (g m^-^3)'])
    
    set(gcf,'color','w')
    I = getframe(h);
    imwrite(I.cdata,figname,'PNG');
  