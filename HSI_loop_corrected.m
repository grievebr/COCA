%% Script to acquire Habitat Suitability Index for species ensemble
% Set Season and years to calculate

yearvec = [1980:2014];
season = 2;

if  season==2
    bias = newspring;
    seasname = 'spring';
elseif season==4
    bias = newfall;
    seasname = 'fall';
end

addpath C:/Users/brian.grieve/Documents/m_map/m_map
addpath C:/Users/brian.grieve/Documents/scripts/matlab/functions
specieslist = readtable('C:/Users/brian.grieve/Documents/coca/BellR_species_listConversions.xlsx');
specieslist{18,'CommonName'} = {'Windowpane flounder'};

load('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/latlonmats.mat');
load('C:/Users/brian.grieve/Documents/Marta_benthicdata/latlon/tempbiases.mat');
cruiseinfo = readtable('C:/Users/brian.grieve/Documents/coca/Data/cruiseinfo.csv');
load('alltrawl_mat.mat');
u_svspp = specieslist{:,'SVSPP'};

% Model Domain
domlat = [35 45];
domlon = [-76 -65];
moddomlat = find(latgrid(1,:)>=domlat(1) & latgrid(1,:)<=domlat(2));
moddomlon = find(longrid(:,1)>=domlon(1) & longrid(:,1)<=domlon(2));





spp_catchability = nan([25 length(yearvec)]);
%%
for ind_vec = 1:length(u_svspp)
    spp = u_svspp(ind_vec);
    spp_name = string(specieslist{ind_vec,1});
    disp({spp,spp_name})
    
    % Load benthic maps
    fname = sprintf('C:/Users/brian.grieve/Documents/coca/Benthic_maps/s%d/binomial/pred%d.txt',season,spp);
    predmat = readmatrix(fname);
    predmat(find(predmat==-999)) = nan;
    predmat = reshape(predmat,size(longrid));
    
    % Normalize to 2SD
    PAsd = nanmean(predmat(:)) + 2*nanstd(predmat(:));
    predmat_std = predmat;
    predmat_std(find(predmat_std>PAsd)) = PAsd;
    predmat_std = predmat_std/PAsd;
    
    % Pull partition
    fname = sprintf('C:/Users/brian.grieve/Documents/coca/Data/cval%d/habmodallHB%d.csv',season,season);
    cvaldf = readtable(fname);
    partition = cvaldf{ind_vec,'partition'};
    fishtrawl = alltrawl(find(alltrawl{:,'SVSPP'}==spp & alltrawl{:,'SEASON'}==season),:);
    
    % Set Johnson-Lewin coefficients
    coefs = csvread(sprintf('C:/Users/brian.grieve/Documents/coca/Rich_GAMS/tmb/%scoefs/tmbcoefs%d.csv',seasname,spp),1,1);
    Topt = coefs(partition,2);
    Er = coefs(partition,3);
    Ed = coefs(partition,4);
    cc = coefs(partition,5);
    ccexp = coefs(partition,6);
    maxabund = BAfunction(Topt,Topt,Er,Ed,cc,ccexp); % Abundance at optimal temperature
    
    y_availability = nan([length(yearvec) 2]);
    y_availability(:,1) = yearvec;
    
    for year=yearvec
        disp(year)
        
        % Isolate cruises per season and year
        ic = find(cruiseinfo{:,'EST_YEAR'}==year & cruiseinfo{:,'SEASON'}==season);
        cruiseno = cruiseinfo{ic,'CRUISE6'};
        ic2 = find(strcmp(string(fishtrawl{:,'CRUISE6'}),cruiseno)==1);
        
        % Days of cruise
        cruisedates = unique(fishtrawl{ic2,'date_id'});
        
        d_availability = nan([length(cruisedates),1]);
        
        for d = 1:length(cruisedates)
            dateid = char(string(cruisedates(d)));
            filename = sprintf('E:/ROMS_COBALT_DOWNSCALE2/%d/bwt_%s.mat',year,dateid);
            
            % Import and debias temperature models
            temp = importdata(filename,'cobalt');
            temp = temp+bias;
            
            % Predict abundance
            temppred = BAfunction(temp,Topt,Er,Ed,cc,ccexp);
            %ombine models
            hsi = (sqrt(temppred/maxabund).*predmat_std);
            % sum HSI across domain for the day
            hsi_domain = hsi(moddomlon,moddomlat);
            sumhsi = nansum(hsi_domain(:))*.0169; % .0169 is original value, .0121 is modified
            
            % pull individual samples for each day
            i_samples = find(strcmp(dateid,string(fishtrawl{:,'date_id'}))==1 & strcmp('NOAA',string(fishtrawl{:,'Department'}))==1);
            daysamples = fishtrawl(i_samples,:);
            daysamples = daysamples(find(~isnan(daysamples{:,'stratarea'})),:);
            awsi = nan([height(daysamples) 1]); % Area Weighted Suitability Index
            if length(awsi)==0
                continue
            else
                for k = 1:height(daysamples)
                    % Find HSI for survey point
                    i_x = find(abs(daysamples{k,'LON'}-longrid(:,1))==nanmin(abs(daysamples{k,'LON'}-longrid(:,1))));
                    i_y = find(abs(daysamples{k,'LAT'}-latgrid(1,:))==nanmin(abs(daysamples{k,'LAT'}-latgrid(1,:))));
                    hsi_ijk = hsi(i_x,i_y);
                    
                    % Pull strata area and number of times it was sampled per cruise
                    stratarea = daysamples{k,'stratarea'};
                    stratpercruise = daysamples{k,'stratpercruise'};
                    awsi(k) = hsi_ijk*(stratarea/stratpercruise); % Area weighted suitability index, Kohout et al 2015
                end
            end
            % sum total area surveyed/potential HSI that day
            d_availability(d) = nansum(awsi/sumhsi);
        end
        % yearly availability
        y_availability(year-min(yearvec)+1,2) = nansum(d_availability);
    end
    % save for species
    spp_catchability(ind_vec,:) = y_availability(:,2);
end
