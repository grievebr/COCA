%% Downscale ROMS to benthic data
load('C:/Users/brian.grieve/Documents/Ocean_models/ROMS_COBALT/downscaling_start2.mat')
% Subset
latsub = latrom(381:600,216:295);
lonsub = lonrom(381:600,216:295);

% total years of ROMS
yearvec = [1980:2014];
for y=yearvec
    % Determine dates of survey    
    y_dateid = alltrawl{find(alltrawl{:,'EST_YEAR'}==y),'date_id'};
    y_season = alltrawl{find(alltrawl{:,'EST_YEAR'}==y),'SEASON'};
    spring_min = min(y_dateid(y_season==2));
    spring_max = max(y_dateid(y_season==2)); 
    fall_min = min(y_dateid(y_season==4));
    fall_max = max(y_dateid(y_season==4));
    % Format dates to Julian
    days = [datetime(string(spring_min),'InputFormat','yyyyMMdd'), datetime(string(spring_max),'InputFormat','yyyyMMdd'), datetime(string(fall_min),'InputFormat','yyyyMMdd'), datetime(string(fall_max),'InputFormat','yyyyMMdd')];
    days = juliandate(days,'modifiedjuliandate')-juliandate(datetime(sprintf('%d-01-01 00:00:00',y)),'modifiedjuliandate')+1;
    days2 = [days(1):days(2), days(3):days(4)];
    
    % Convert model time to Julian
    filename = strcat('C:/Users/brian.grieve/Documents/Ocean_models/bwt_1980-2014/bwt_1980-2014/bwt_',sprintf('%d',y),'.nc'); %disp(filename);
    time = ncread(filename,'ocean_time');
    datevec = datetime('1900-01-01')+seconds(time(:));
    jdays = juliandate(datevec,'modifiedjuliandate')-juliandate(datetime(sprintf('%d-01-01 01:00:00',y)),'modifiedjuliandate')+1;
    dayvec = nan(length(days2),1); % The stupidest for loop I've ever written.
    for d = 1:length(days2)
        dayvec(d) = find(jdays==days2(d));
    end
    
    temp = ncread(filename,'bwt');
    temp = permute(temp, [3 1 2]);
    for d = 1:length(dayvec)  
        % Format dateid
        if day(datevec(dayvec(d)))<10 
            daystr = strcat('0',string(day(datevec(dayvec(d)))));
        else
            daystr = string(day(datevec(dayvec(d))));
        end
        if month(datevec(dayvec(d)))<10 
            monstr = strcat('0',string(month(datevec(dayvec(d)))));
        else
            monstr = string(month(datevec(dayvec(d))));
        end
        fileout = strcat(sprintf('E:/ROMS_COBALT_DOWNSCALE2/%d',y),sprintf('/bwt_%d%s%s.mat',y,monstr,daystr)); disp(fileout)
        slice = squeeze(temp(dayvec(d),:,:));
        slice(slice>1E36) = nan; % nan values
        slicepad = nan(size(latrom)); 
        slicepad(381:600,216:295) = slice;
        % downscale with nearest neighbor
        cobalt = griddata(lonrom,latrom,slicepad,lonmarta,latmarta,'nearest');
        save(fileout,'cobalt');
        clear cobalt
    end
    
    clear temp cobalt filename fileout
end
%% Regrid habitat variables to same, most frequent resolution 

load('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/latlonmats.mat'); 

% Bathymetry
bathylat = csvread('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/bathylat.csv',1,0);
bathylon = csvread('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/bathylon.csv',1,0);
bathydata = csvread('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/bathymat.csv',1,0);
bathydata(bathydata==0) = nan;
bathyres = griddata(bathylon,bathylat,bathydata,lonmarta,latmarta,'nearest'); 
save('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/bathyres.mat','bathyres');
csvwrite('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/bathyres.csv',bathyres)

% Hard Bottom
HBlat = csvread('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/HBlat.csv',1,0);
HBlon = csvread('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/HBlon.csv',1,0);
HBdata = csvread('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/HBmat.csv',1,0);
HBdata(HBdata==-999) = nan;
HBres = griddata(HBlon,HBlat,HBdata,lonmarta,latmarta,'nearest'); 
save('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/HBres.mat','HBres');
csvwrite('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/HBres.csv',HBres)

% BPI Broad
bpilat = csvread('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/bpilat.csv',1,0);
bpilon = csvread('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/bpilon.csv',1,0);
bpidata = csvread('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/bpimat.csv',1,0);
bpidata(bpidata==-999) = nan;
bpires = griddata(bpilon,bpilat,bpidata,lonmarta,latmarta,'nearest'); 
save('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/bpires.mat','bpires');
csvwrite('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/bpires.csv',bpires)

% BPI fine. Uses same resolution as BPI Broad
bpidata = csvread('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/bpifinemat.csv',1,0);
bpidata(bpidata==-999) = nan;
bpifineres = griddata(bpilon,bpilat,bpidata,lonmarta,latmarta,'nearest'); 
save('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/bpifineres.mat','bpifineres');
csvwrite('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/bpifineres.csv',bpifineres)

% Sediment classes
sedlat = csvread('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/sedlat.csv',1,0);
sedlon = csvread('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/sedlon.csv',1,0);
seddata = csvread('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/sedmat.csv',1,0);
sedres = griddata(sedlon,sedlat,seddata,lonmarta,latmarta,'nearest'); 
save('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/sedres.mat','sedres');
csvwrite('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/sedres.csv',sedres)

% Patch Dynamics
compactmat = csvread('C:/Users/brian.grieve/Documents/coca/data/compactmat.csv',1,0); 
save('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/compact.mat','compactmat');
areamat = csvread('C:/Users/brian.grieve/Documents/coca/data/areamat.csv',1,0); 
save('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/area.mat','areamat');
perimmat = csvread('C:/Users/brian.grieve/Documents/coca/data/perimmat.csv',1,0); 
save('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/perim.mat','perimmat');



load('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/bathyres.mat');
load('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/bpires.mat');
load('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/bpifineres.mat');
load('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/HBres.mat');
load('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/sedres.mat');
load('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/HBres.mat');
load('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/area.mat');
load('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/compact.mat');
load('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/perim.mat');


habvarnames = ['lat';'lon';'bat';'bpi';'sed';'hbo';'ara';'per';'com'];
habvars = [latmarta(:),lonmarta(:)];
habvars(:,3) = bathyres(:);
habvars(:,4) = bpires(:);
habvars(:,5) = sedres(:);
habvars(:,6) = HBres(:);
habvars(:,7) = areamat(:);
habvars(:,8) = patchmat(:);
habvars(:,9) = compactmat(:);

csvwrite('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/habvars.csv',habvars);

%% Regrid temperature bias
tempbias_spring = readmatrix('C:/Users/brian.grieve/Documents/coca/Data/tempbiasspring_newmod.txt');
tempbias_fall = readmatrix('C:/Users/brian.grieve/Documents/coca/Data/tempbiasfall_newmod.txt');
lonsub = readmatrix('C:/Users/brian.grieve/Documents/coca/Data/lonsub.txt');
latsub = readmatrix('C:/Users/brian.grieve/Documents/coca/Data/latsub.txt');
lonsub = lonsub-360;

springbias = griddata(lonsub,latsub,tempbias_spring,lonmarta,latmarta,'nearest'); 
fallbias = griddata(lonsub,latsub,tempbias_fall,lonmarta,latmarta,'nearest'); 

save('C:/Users/brian.grieve/Documents/coca/Data/ROMS_bias.mat','fallbias', 'springbias')

