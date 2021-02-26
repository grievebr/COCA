library(mgcv); library(chron); library(fishmethods); library(readxl); library(stringr); library(tidyr); library(mefa)

# Need a mode function
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Calibration factors between surveys
species_conversions = read_excel("C:/Users/brian.grieve/Documents/coca/BellR_species_listConversions.xlsx",sheet = "Sheet1")
species_conversions[,15:18] = 0;
names(species_conversions) = c('CommonName','ScientificName',names(species_conversions[3:14]),'HB_pw_SPRING','HB_pw_FALL','HB_pW_SPRING','HB_pW_FALL')

bigelow_fall = read.csv("C:/Users/brian.grieve/Downloads/bigelow_fall_calibration.csv")
bigelow_fall = bigelow_fall[is.element(bigelow_fall$svspp,u_svspp),1:5]
bigelow_spring = read.csv("C:/Users/brian.grieve/Downloads/bigelow_spring_calibration.csv")
bigelow_spring = bigelow_spring[is.element(bigelow_spring$svspp,u_svspp),1:5]
#species_conversions = species_conversions[-which(is.na(species_conversions$SVSPP)),];
species_conversions[(is.na(species_conversions)==1 | species_conversions==0)] = 1; 
for (spp in u_svspp){
  i_hb = which(bigelow_spring$svspp==spp)
  i_spp = which(species_conversions$SVSPP==spp)
  species_conversions[i_spp,'HB_pw_SPRING'] = bigelow_spring[i_hb,'pw'];
  species_conversions[i_spp,'HB_pW_SPRING'] = bigelow_spring[i_hb,'pW'];
  
  i_hb = which(bigelow_fall$svspp==spp)
  i_spp = which(species_conversions$SVSPP==spp)
  species_conversions[i_spp,'HB_pw_FALL'] = bigelow_fall[i_hb,'pw'];
  species_conversions[i_spp,'HB_pW_FALL'] = bigelow_fall[i_hb,'pW'];
}


# Format trawl survey data
# Prior to 2018 NOAA Bottom Trawl
noaatrawl = read.csv("C:/Users/brian.grieve/Documents/coca/Bell_Data_Request_2018_bis_20180727.txt")
varnames = names(noaatrawl);

# 2018 NOAA Trawl data
noaatrawl$AREA = NA; 
Bell_Data_Request_2019 <- read.csv("C:/Users/brian.grieve/Documents/coca/Data/Bell_Data_Request_2019.txt")
newi = vector(mode='numeric',length=NCOL(Bell_Data_Request_2019))
for (i in 1:length(newi)){newi[i] = which(names(noaatrawl)==names(Bell_Data_Request_2019)[i])}
noaa19 = data.frame(matrix(data=NA,nrow=nrow(Bell_Data_Request_2019),ncol=ncol(noaatrawl))); names(noaa19) = names(noaatrawl)
noaa19[,newi] = Bell_Data_Request_2019;
noaatrawl = rbind(noaatrawl,noaa19[which(noaa19$CRUISE6==201804),]);
noaatrawl = noaatrawl[,varnames];

# Information about trawl types to calculate area swept
geardata = read_excel("C:/Users/brian.grieve/Documents/coca/Trawl_gear.xlsx",col_types = c("numeric", "text", "date","date", "text", "text", "text", "text"))

# Code Seasons to numbers
n = nrow(noaatrawl); 
seas = vector(mode='double',length=n)
seas[which(noaatrawl$SEASON=='FALL')]=4; 
seas[which(noaatrawl$SEASON=='SPRING')]=2;
seas[which(noaatrawl$SEASON=='SUMMER')]=3;
seas[which(noaatrawl$SEASON=='WINTER')]=1;
seas[which(noaatrawl$EST_MONTH==4 & seas==1)]=2;
noaatrawl$SEASON=seas;

seas = vector(mode='double',length=nrow(geardata))
seas[which(geardata$Season=='Fall')]=4; 
seas[which(geardata$Season=='Spring')]=2;
seas[which(geardata$Season=='Summer')]=3;
seas[which(geardata$Season=='Winter')]=1;
geardata$Season=seas;

# Reformat abundance, biomass, and PA vectors into dataframe with SVSPP as a single column
biovec = seq(32,106,3); #seq(18,92,3);
abunvec = seq(33,106,3); #seq(19,92,3); 
pavec = seq(34,106,3); #seq(20,92,3); 
biovec1 = seq(biovec[1],length(names(noaatrawl)),2)
abunvec1 = seq(abunvec[1],length(names(noaatrawl)),2)

varnames = names(noaatrawl); 
newnames = vector(mode='character',length(names(biovec))); 
newnames[1:(biovec[1]-1)] = varnames[1:(biovec[1]-1)]; 
newnames[biovec] = varnames[biovec1]; 
newnames[abunvec] = varnames[abunvec1]; 
for (i in pavec){
  newnames[i] = sub('_WT','_PA',newnames[i-2])}


#create new dataframe with P/A
trawl = data.frame(matrix(data=NA,nrow=n,ncol=max(pavec)))
trawl[,1:(biovec[1]-1)] = noaatrawl[,1:(biovec[1]-1)]
trawl[,biovec] = noaatrawl[,biovec1]; 
trawl[,abunvec] = noaatrawl[,abunvec1]; 
for (i in biovec){
  trawl[which(is.na(trawl[,i])),i] = 0;
  i1 = i+1
  trawl[which(is.na(trawl[,i1])),i1] = 0;
  i2 = i+2; 
  trawl[,i2] = trawl[,i1]>0
}
names(trawl) = newnames; 

# Fill in missing trawl information from sources
trawl[(which(is.na(trawl$TOWDUR) & trawl$SV_VESSEL=='DE'|trawl$SV_VESSEL=='AL'|trawl$SV_VESSEL=='AT')),'TOWDUR']=30;
trawl[(which(is.na(trawl$TOWDUR) & trawl$SV_VESSEL=='HB')),'TOWDUR']=20; # Politis2014
trawl$MEAN_WING_SPRD_METERS[which(is.na(trawl$MEAN_WING_SPRD_METERS) & trawl$SVGEAR==10)] = mean(trawl$MEAN_WING_SPRD_METERS[which(trawl$SVGEAR==10)],na.rm=T)
trawl[which(trawl$BOTSPEED==0 & trawl$SV_VESSEL!= 'HB'),'BOTSPEED']=3.8 # 3.8 Standard (Stauffer 2003, Jech 2014);
trawl[which(trawl$BOTSPEED==0 & trawl$SV_VESSEL== 'HB'),'BOTSPEED']=3; # Politis 2014
# Net sizes (Azarovitz 1981); 
trawl$MEAN_WING_SPRD_METERS[which(trawl$SVGEAR==11)]=10.4; 
trawl$MEAN_WING_SPRD_METERS[which(trawl$SVGEAR==12)]=9.2;
trawl$MEAN_WING_SPRD_METERS[which(trawl$SVGEAR==41)]=11.8;
trawl$MEAN_WING_SPRD_METERS[which(trawl$SVGEAR==17)]=10.4; # Modified 36 Yankee trawl with chain cookie/rubber sweep

calc_sa_wing = ((trawl$TOWDUR*60)*trawl$BOTSPEED*.5144*trawl$MEAN_WING_SPRD_METERS)/1E6; # (Seconds)*(knots)to(m/s)*(m)to(km) = Swept area in km.


# Unique values
u_cruiseID = na.omit(unique(trawl$CRUISE6))
u_years = na.omit(unique(trawl$EST_YEAR)) 
u_svspp = na.omit(unique(species_conversions$SVSPP))

# species-specific conversion factors based on gear used in specific cruises
cruise_gear = expand.grid(u_svspp,u_cruiseID); 
cruise_gear = data.frame(cbind(cruise_gear$Var2,cruise_gear$Var1,matrix(data=NA,nrow=nrow(cruise_gear),ncol = 15)))
names(cruise_gear) = c('cruiseID','SVSPP','season','year','vessel','vessel_conversion_NUM','vessel_conversion_WT','trawltype','svgear','gear_conversion_NUM','gear_conversion_WT','doors','door_conversion_NUM','door_conversion_WT','swept_area','total_conversion_NUM','total_conversion_WT')
for (i in u_cruiseID){
  ind_trawl = which(trawl$CRUISE6==i)
  ind_cruisegear = which(cruise_gear$cruiseID==i)
  if (length(unique(cruise_gear[ind_cruisegear,'EST_YEAR']))>1){warning(sprintf('Cruise spanned multiple years %d',i))}
  if (length(unique(cruise_gear[ind_cruisegear,'SEASON']))>1){warning(sprintf('Cruise spanned multiple seasons %d',i))}
  if (length(unique(cruise_gear[ind_cruisegear,'SV_VESSEL']))>1){warning(sprintf('Cruise spanned multiple vessels %d',i))}
  if (length(unique(cruise_gear[ind_cruisegear,'SVGEAR']))>1){warning(sprintf('Cruise spanned multiple gears %d',i))}
  cruise_gear[ind_cruisegear,'year'] = getmode(trawl[ind_trawl,'EST_YEAR']); 
  cruise_gear[ind_cruisegear,'season'] = getmode(trawl[ind_trawl,'SEASON']);
  cruise_gear[ind_cruisegear,'vessel'] = getmode(trawl[ind_trawl,'SV_VESSEL']); 
  cruise_gear[ind_cruisegear,'svgear'] = getmode(trawl[ind_trawl,'SVGEAR']);
  cruise_gear[ind_cruisegear,'swept_area'] = mean(calc_sa_wing[ind_trawl],na.rm=T)
}
cruise_gear$trawltype[which(cruise_gear$svgear==10)] = '4 Seam, 3 Bridle'; #(Politis 2014)
cruise_gear$trawltype[which(cruise_gear$svgear==11)] = '36 Yankee';
cruise_gear$trawltype[which(cruise_gear$svgear==12)] = '3/4 36 Yankee';
cruise_gear$trawltype[which(cruise_gear$svgear==17)] = 'Mod. 36 Yankee chainsweep';
cruise_gear$trawltype[which(cruise_gear$svgear==41)] = '41 Yankee';


# Conversion Factors to NEFSC Survey standard configuration: Reid 1999 (species_conversions table)
# Trawls: #41 Yankee (svgear==41) to #36 Yankee (svgear==11) - Spring 1973-1981 only
# Doors: BMV to Polyvalent - Spring 1985 to present
# Vessels: Delaware II to Albatross IV - Various, some during same survey
cruise_gear$doors[which(cruise_gear$year>=1985)] = 'polyvalent'
cruise_gear$doors[c(which(cruise_gear$year<=1984),which(cruise_gear$year==1985 & cruise_gear$season==1))] = 'BMV';
cruise_gear[which(cruise_gear$svgear==11),c('gear_conversion_NUM','gear_conversion_WT')]=1
cruise_gear[which(cruise_gear$doors=='polyvalent'),c('door_conversion_NUM','door_conversion_WT')]=1
cruise_gear[which(cruise_gear$vessel=='AL'),c('vessel_conversion_NUM','vessel_conversion_WT')]=1 # Setting to Albatross units
cruise_gear[which(cruise_gear$vessel=='AT'),c('vessel_conversion_NUM','vessel_conversion_WT')]=1 # No conversion (reid 1999)
cruise_gear[which(cruise_gear$vessel=='PC'),c('vessel_conversion_NUM','vessel_conversion_WT')]=1 # No conversion (used fall 2017 only, similar to HB)
cruise_gear[which(cruise_gear$svgear==12 | cruise_gear$svgear==17 | cruise_gear$svgear==10),c('gear_conversion_NUM','gear_conversion_WT')]=1 # No conversion (reid 1999); HB gear accounted for in vessel conversion

# apply for species
for (i in u_svspp){
  ind_spp_conv = which(species_conversions$SVSPP==i);
  ind_spp = which(cruise_gear$SVSPP==i)
  ind_cruisegear = which(cruise_gear$SVSPP==i & cruise_gear$svgear==41)
  cruise_gear[ind_cruisegear,'gear_conversion_NUM'] = species_conversions[ind_spp_conv,'GCF_NUM'];
  cruise_gear[ind_cruisegear,'gear_conversion_WT'] = species_conversions[ind_spp_conv,'GCF_WT'];
  ind_cruisegear = which(cruise_gear$SVSPP==i & cruise_gear$doors=='BMV' & cruise_gear$vessel!='HB')
  cruise_gear[ind_cruisegear,'door_conversion_NUM'] = species_conversions[ind_spp_conv,'DCF_NUM'];
  cruise_gear[ind_cruisegear,'door_conversion_WT'] = species_conversions[ind_spp_conv,'DCF_WT'];
  ind_cruisegear = which(cruise_gear$SVSPP==i & cruise_gear$vessel=='DE')
  cruise_gear[ind_cruisegear,'vessel_conversion_NUM'] = species_conversions[ind_spp_conv,'VCF_NUM'];
  cruise_gear[ind_cruisegear,'vessel_conversion_WT'] = species_conversions[ind_spp_conv,'VCF_WT'];
  ind_cruisegear = which(cruise_gear$SVSPP==i & cruise_gear$vessel=='HB' & cruise_gear$season==2)
  cruise_gear[ind_cruisegear,'vessel_conversion_NUM'] = 1/species_conversions[ind_spp_conv,'HB_NUM_SPRING'];
  cruise_gear[ind_cruisegear,'vessel_conversion_WT'] = 1/species_conversions[ind_spp_conv,'HB_pW_SPRING'];
  cruise_gear[ind_cruisegear,'door_conversion_NUM'] = 1;
  cruise_gear[ind_cruisegear,'door_conversion_WT'] = 1;
  ind_cruisegear = which(cruise_gear$SVSPP==i & cruise_gear$vessel=='HB' & cruise_gear$season==4)
  cruise_gear[ind_cruisegear,'vessel_conversion_NUM'] = 1/species_conversions[ind_spp_conv,'HB_NUM_FALL'];
  cruise_gear[ind_cruisegear,'vessel_conversion_WT'] = 1/species_conversions[ind_spp_conv,'HB_pW_FALL'];
  cruise_gear[ind_cruisegear,'door_conversion_NUM'] = 1;
  cruise_gear[ind_cruisegear,'door_conversion_WT'] = 1;
}
  
# Get final conversion. Divide by gear, multiply by door & vessel (http://globec.whoi.edu/globec-dir/data_doc/NMFS_trawl_factors.html)
cruise_gear[,'total_conversion_NUM'] = (cruise_gear[,'door_conversion_NUM']*cruise_gear[,'vessel_conversion_NUM'])/cruise_gear[,'gear_conversion_NUM'];
cruise_gear[,'total_conversion_WT'] = (cruise_gear[,'door_conversion_WT']*cruise_gear[,'vessel_conversion_WT'])/cruise_gear[,'gear_conversion_WT'];


trawl$ind = 1:nrow(trawl)
trawl.corrected = rep(trawl[,c(1:31,107)],times=length(u_svspp));
trawl.corrected$SVSPP = rep(u_svspp,each=nrow(trawl))
trawl.corrected$num = as.vector(unlist(trawl[,abunvec]));
trawl.corrected$wt = as.vector(unlist(trawl[,biovec]));
trawl.corrected$presence = as.vector(unlist(trawl[,pavec]));

for (i in u_cruiseID){
  for (u in u_svspp){
    tc.ind = which(trawl.corrected$CRUISE6==i & trawl.corrected$SVSPP==u)
    cg.ind = which(cruise_gear$cruiseID==i & cruise_gear$SVSPP==u)
    trawl.corrected[tc.ind,'AREA_SWEPT_WINGS_MEAN_KM2'] = cruise_gear[cg.ind,'swept_area']
    trawl.corrected[tc.ind,'num_raw'] = trawl.corrected[tc.ind,'num'];
    trawl.corrected[tc.ind,'wt_raw'] = trawl.corrected[tc.ind,'wt']
    trawl.corrected[tc.ind,'num'] = trawl.corrected[tc.ind,'num']*cruise_gear[cg.ind,'total_conversion_NUM'];
    trawl.corrected[tc.ind,'wt'] = trawl.corrected[tc.ind,'wt']*cruise_gear[cg.ind,'total_conversion_WT'];
  }}

# Datetime formatting and night determination
esttime = chron(times.=trawl.corrected$EST_TIME,format='h:m:s');
trawl.corrected$time.num = hours(esttime)+minutes(esttime)/60
zenith.df = astrocalc4r(day=trawl.corrected$EST_DAY,month=trawl.corrected$EST_MONTH,year=trawl.corrected$EST_YEAR,hour=trawl.corrected$time.num,timezone=c(rep(-5,times=nrow(trawl.corrected))),lat=trawl.corrected$LAT,lon=trawl.corrected$LON)
trawl.corrected$zenith = zenith.df$zenith
trawl.corrected$night = trawl.corrected$time.num>zenith.df$sunset | trawl.corrected$time.num<zenith.df$sunrise
    
# Now add state surveys. Each were a different format

# Massachusetts DMF
sheet_ids = c('015','023','026','032','033','034','036','072','073','074','075','076','101','102','103','105','106','107','108','139','141','155','177','193','301')
for (i in 1:length(sheet_ids)){
sheet_fall = read_excel("C:/Users/brian.grieve/Documents/coca/statedata/Mass/MDMF_FALL_RB.xlsx",sheet=paste('F_',sheet_ids[i],sep=''))
sheet_spring = read_excel("C:/Users/brian.grieve/Documents/coca/statedata/Mass/MDMF_SPRING_RB.xlsx",sheet=paste('S_',sheet_ids[i],sep=''))
if (i==1){
  MDMF = rbind(sheet_fall,sheet_spring)}
if (i!=1){
MDMF = rbind(MDMF,sheet_fall,sheet_spring);}} 
names(MDMF) = c('LON','LAT','YEAR','SEASON','STRATUM','STATION','AVGDEPTH','BTEMP','EST_YEAR','EST_MONTH','EST_DAY','EST_TIME','SVSPP','num','wt')

seas = vector(mode='double',length=nrow(MDMF))
seas[which(MDMF$SEASON=='FALL')]=4; 
seas[which(MDMF$SEASON=='SPRING')]=2;
MDMF$SEASON=seas;
MDMF$EST_TIME=chron(times. = str_remove(as.character(MDMF$EST_TIME),'1899-12-31 '),format = c(times='h:m:s'));
#

MDMF$Department = 'Massachusetts'
MDMF$MEAN_WING_SPRD_METERS = 9.2 # Arovitz 1982. Slightly different head/footrope lengths, so may not be perfect
MDMF$BOTSPEED = 2.5 #USFWS 2016 Performance report
MDMF$TOWDUR = 20
MDMF$AREA_SWEPT_WINGS_MEAN_KM2 = (MDMF$MEAN_WING_SPRD_METERS*MDMF$BOTSPEED*.5144*MDMF$TOWDUR*60)/1E6

# Connecticut 
sheet_ids = c('Cod','Haddock','Yellowtail flounder','Winter flounder','Summer flounder','American Plaice','Atlantic herring','Pollock','Menhaden','Striped bass','Spiny dogfish','Black sea bass','Winter skate','Little skate','Windowpane flounder','Ocean pout','Silver hake','Tautog','Lobster','Alewife','Blueback herring')
conn_svspp = c(73,74,105,106,103,102,32,75,36,139,15,141,23,26,108,193,72,177,301,33,34)
for (i in 1:length(sheet_ids)){
  sheet = read_excel("C:/Users/brian.grieve/Documents/coca/statedata/RichBell_LISTS_Spp1984-2017.xlsx",sheet=sheet_ids[i],skip=2,col_types = c("text", "date", "numeric", 
                     "text", "numeric", "text", "text","numeric", "numeric", "text", "date","numeric", "numeric", "numeric","numeric","text", "numeric", "numeric", 
                     "numeric", "numeric", "numeric","text", "text", "text", "text", "text","text", "text", "text", "numeric","numeric", "text", "numeric", "text", 
                     "text", "numeric", "numeric", "numeric","numeric", "numeric", "numeric","numeric"))
  sheet$SVSPP = conn_svspp[i]
  if (i==1){
    conn = sheet;} 
  if (i!=1){
    conn = rbind(conn,sheet);}
  }
seas = vector(mode='double',length=nrow(conn))
seas[which(conn$Season=='FA')]=4; 
seas[which(conn$Season=='SP')]=2;
conn$Season=seas;
conn$AREA_SWEPT_WINGS_MEAN_KM2 = conn$`A_swept(sqnmi)`*3.4345; # nmi2 to km2
conn[which(is.na(conn$AREA_SWEPT_WINGS_MEAN_KM2)),'AREA_SWEPT_WINGS_MEAN_KM2'] = (30*60*3.5*.5144*9.3477)/1E6; #Duration(min) * seconds * speed(knots) * meters/sec * netwidth / km2. Does not account for tow length standardization because abundance already is. 

newnames = names(conn);
newnames[c(3,4,5,6,12,13,14,15,20,33,34,40,41)] = c('EST_YEAR','SEASON','STRATUM','CRUISE6','TOWDUR','LAT','LON','AVGDEPTH','BTEMP','num','wt','BOTSPEED','DOPDISTB')
names(conn) = newnames;
conn[which(conn$wt=='.' & conn$num==0),'wt'] = 0;
conn$Department = 'Connecticut'
conn$Date = chron(as.character(conn$Date),format = c(dates = 'y-m-d'))
conn$EST_MONTH = as.numeric(months(conn$Date))
conn$EST_DAY = as.numeric(days(conn$Date))
conn$EST_TIME = chron(times. = str_remove(as.character(conn$`Time Start`),'1899-12-31 '),format = c(times = 'h:m:s'))
esttime = chron(times.=conn$EST_TIME,format='h:m:s');
conn$time.num = hours(esttime)+minutes(esttime)/60
zenith.df = astrocalc4r(day=conn$EST_DAY,month=conn$EST_MONTH,year=conn$EST_YEAR,hour=conn$time.num,timezone=c(rep(-5,times=nrow(conn))),lat=conn$LAT,lon=conn$LON)
conn$zenith = zenith.df$zenith
conn$night = conn$time.num>zenith.df$sunset | conn$time.num<zenith.df$sunrise

# New Jersey
nj = read_excel("C:/Users/brian.grieve/Documents/coca/statedata/BellRich_TNC_Feb2019/HEAD+CATCH.xlsx", col_types = c("date", "numeric", "numeric", 
      "numeric", "numeric", "numeric", "numeric", "text", "numeric"))
nj$SPECIES[which(nj$SPECIES=='Black Seabass')] = 'Black Sea Bass';
nj$svspp = 0; nj$lon = 0; nj$lat = 0;
nj$STATION = floor(nj$STATION)
coords = read_excel("C:/Users/brian.grieve/Documents/coca/statedata/BellRich_TNC_Feb2019/Station_Coordinates.xlsx", col_types = c("text", 
          "numeric", "numeric", "numeric", "numeric"))
coords[which(coords$Station=='23A'),'Station'] = '23';
coords[which(coords$Station=='29A'),'Station'] = '29';
coords[which(coords$Station=='47A'),'Station'] = '47';
coords$Station = as.numeric(coords$Station)
coords$lat = (coords$Start_Latitude+coords$End_Latitude)/2
coords$lon = (coords$Start_Longitude+coords$End_Longitude)/2
njspec = read_excel("C:/Users/brian.grieve/Documents/coca/statedata/BellRich_TNC_Feb2019/SPECIES - DelBay Trawl.xlsx");
u = unique(njspec$svspp)[-which(is.na(unique(njspec$svspp)))]
for (i in u){
  svspp_spec = njspec[which(njspec$svspp==i),'SPECIES'];
  nj[which(nj$SPECIES==svspp_spec$SPECIES),'svspp'] = i;
}
nj = nj[-which(nj$svspp==0),]

u = unique(nj$STATION)
for (i in u){
  i_station = which(coords$Station==i) 
  nj[which(nj$STATION==i),'lat'] = coords[i_station,'lat']
  nj[which(nj$STATION==i),'lon'] = coords[i_station,'lon']
}

nj2 = expand(nj,crossing(nesting(DATE,STATION),svspp))
u_stationdate = nesting(nj$DATE,nj$STATION);
nj2[,c(4:12)]=NA;
names(nj2) = c('DATE','STATION','SVSPP','EST_MONTH','AVGDEPTH','SALINITY','BTEMP','DO','SPECIES','num','LON','LAT')
for (i in 1:nrow(u_stationdate)){
  ind1 = which(nj$DATE==as.numeric(u_stationdate[i,1]) & nj$STATION==as.numeric(u_stationdate[i,2]))
  ind2 = which(nj2$DATE==as.numeric(u_stationdate[i,1]) & nj2$STATION==as.numeric(u_stationdate[i,2]))
  for (i2 in ind2){
    nj2[i2,c(4:9,11,12)] = nj[ind1[1],c(2,4:8,11,12)]}
  for (i2 in ind1){
    nj2[ind2[which(nj2[ind2,'SVSPP']==as.numeric(nj[i2,'svspp']))],'num'] = as.numeric(nj[i2,'NUMBER'])}
}
nj2[which(is.na(nj2$num)),'num'] = 0;
nj2$EST_YEAR = as.numeric(as.character(years(nj2$DATE)))
nj2$EST_DAY = as.numeric(as.character(days(nj2$DATE)))
nj2$SEASON = NA;

# NJ offshore
nj.o = read_excel("C:/Users/brian.grieve/Documents/coca/statedata/njoceantrawlsurveydata/RBellSpeciesCatch.xlsx")
nj.o = nj.o[-which(nj.o$SPP==999),]
nj.o$EST_YEAR = as.numeric(str_sub(nj.o$YRMODA,1,4));
nj.o$EST_MONTH = as.numeric(str_sub(nj.o$YRMODA,5,6));
nj.o$EST_DAY = as.numeric(str_sub(nj.o$YRMODA,7,8));
nj.o$LAT = rowMeans(cbind(nj.o$SLAT,nj.o$ELAT),na.rm=T)/100
nj.o$LON = rowMeans(cbind(nj.o$SLONG,nj.o$ELONG),na.rm=T)/100*-1
nj.o$AVGDEPTH = rowMeans(cbind(nj.o$MINDEPTH,nj.o$MAXDEPTH),na.rm=T);
nj.o$AVGDEPTH_BEGEND = rowMeans(cbind(nj.o$STARTDEPTH,nj.o$ENDDEPTH),na.rm=T); 
nj.o$SEASON = NA;

newnames = names(nj.o)
newnames[c(2,8,11,12,13,24,29,32,33)] = c('CRUISE6','SV_VESSEL','EST_TIME','TOWDUR','DOPDISTB','BTEMP','SVSPP','num','wt')
names(nj.o) = newnames;
nj.o[which(is.na(nj.o$num)),'num'] = 0;
nj.o[which(is.na(nj.o$wt)),'wt'] = 0;

nj2$Department = 'NJ'
nj2$AREA_SWEPT_WINGS_MEAN_KM2 = 2.6*.5144*(10*60)*16*.55/1E6 # 2.6 knots to m/s * (16'otter*.55constant) to km2. Alternatively, use .4NM*net. Email from Andrew Hassell
nj.o$AREA_SWEPT_WINGS_MEAN_KM2 = 2.6*.5144*(20*60)*30.5*.55/1E6
nj.o$Department = 'NJ.o'
nj.o = nj.o[-which(nj.o$LON>-70),]
i = which(str_length(nj.o$EST_TIME)==4); 
nj.o$EST_TIME[i] = paste(str_sub(nj.o$EST_TIME[i],1,2),':',str_sub(nj.o$EST_TIME[i],3,4),':00',sep='')
i = which(str_length(nj.o$EST_TIME)==3)
nj.o$EST_TIME[i] = paste('0',str_sub(nj.o$EST_TIME[i],1,1),':',str_sub(nj.o$EST_TIME[i],2,3),':00',sep='')
nj.o$EST_TIME = chron(times.= nj.o$EST_TIME, format = 'h:m:s')
esttime = chron(times.=nj.o$EST_TIME,format='h:m:s');
nj.o$time.num = hours(esttime)+minutes(esttime)/60
zenith.df = astrocalc4r(day=nj.o$EST_DAY,month=nj.o$EST_MONTH,year=nj.o$EST_YEAR,hour=nj.o$time.num,timezone=c(rep(-5,times=nrow(nj.o))),lat=nj.o$LAT,lon=nj.o$LON)
nj.o$zenith = zenith.df$zenith
nj.o$night = nj.o$time.num>zenith.df$sunset | nj.o$time.num<zenith.df$sunrise

# MENH
spdf = data.frame(matrix(data=NA,nrow=24,ncol=2))
names(spdf) = c('common','svspp')
spdf$common = c('AcadianRedfish','Alewife','AmericanLobster','AmericanPlaice','AtlanticCod','AtlanticHalibut','AtlanticHerring','AtlanticMenhaden',
                'Blackseabass','BluebackHerring','Haddock','LittleSkate','Oceanpout','Pollock','SilverHake','SpinyDogfish','StripedBass','SummerFlounder',
                'WhiteHake','WindowpaneFlounder','WinterFlounder','WinterSkate','WitchFlounder','YellowtailFlounder')
spdf$svspp = c(155,33,301,102,73,101,32,36,141,34,74,26,193,75,72,15,139,103,76,108,106,23,107,105)
for (i in 1:nrow(spdf)){
  fish = read_excel(paste("C:/Users/brian.grieve/Documents/coca/statedata/MeNH/Data Files/",spdf[i,1],".xlsx",sep=''), 
         col_types = c("text", "numeric", "date","numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                       "numeric", "numeric", "numeric", "text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))
  fish$SVSPP = spdf[i,'svspp']
  if (i==1){
    menh = fish;}
  if (i!=1){
    menh = rbind(menh,fish)}}
  menh$LAT = (menh$start_lat+menh$end_lat)/2;
  menh$LON = (menh$start_lon+menh$end_lon)/2;
  menh$AVGDEPTH = rowMeans(cbind(menh$start_depth,menh$end_depth),na.rm=TRUE);
  
  # Depth Strata and SA
  menh[which(menh$AVGDEPTH<=40),'depthstrata'] = 1;
  menh[which(menh$AVGDEPTH>40 & menh$AVGDEPTH<=70),'depthstrata'] = 2;
  menh[which(menh$AVGDEPTH>70 & menh$AVGDEPTH<=100),'depthstrata'] = 3;
  menh[which(menh$AVGDEPTH>100),'depthstrata'] = 4;
  
  # This SA is the average at given depth zones
  menh[which(menh$depthstrata==1),'AREA_SWEPT_WINGS_MEAN_KM2'] = 0.014017
  menh[which(menh$depthstrata==2),'AREA_SWEPT_WINGS_MEAN_KM2'] = 0.015551
  menh[which(menh$depthstrata==3),'AREA_SWEPT_WINGS_MEAN_KM2'] = 0.016052
  menh[which(menh$depthstrata==4),'AREA_SWEPT_WINGS_MEAN_KM2'] = 0.01628
  
  # Average Wing Spread at given depth zones
  menh[which(menh$depthstrata==1),'MEAN_WING_SPRD_METERS'] = 9.23
  menh[which(menh$depthstrata==2),'MEAN_WING_SPRD_METERS'] = 10.24
  menh[which(menh$depthstrata==3),'MEAN_WING_SPRD_METERS'] = 10.57
  menh[which(menh$depthstrata==4),'MEAN_WING_SPRD_METERS'] = 10.72
  
  menh$Start_Date = chron(dates.=str_sub(as.character(menh$Start_Date),1,10),format = c(dates = 'y-m-d'))
  menh$EST_DAY = as.numeric(days(menh$Start_Date))
  menh$EST_MONTH = as.numeric(months(menh$Start_Date));
  menh$EST_YEAR = years(menh$Start_Date)
  menh$Department = 'menh'
  menh$SEASON = NA;
  menh[which(is.na(menh$EXP_Catch)),'EXP_Catch'] = 0;
  menh[which(is.na(menh$EXP_Weight)),'EXP_Weight'] = 0;
  
  newnames = names(menh);
  newnames[c(1,2,5,12,16,17,18)] = c('CRUISE6','TOW','STRATUM','DOPDISTB','num','wt','BTEMP')
  names(menh) = newnames;
  
  
# Delaware
sheet_ids = c('Black Seabass','Tautog','Summer Flounder','Atl. Herring','Atl. Menhaden','Striped Bass','Spiny Dogfish','Winter Skate','Little Skate',
              'Windowpane','Silver Hake','Alewife','Blueback Herring','Winter Flounder')
del_svspp = c(141,177,103,32,36,139,15,23,26,108,72,33,34,106)
for (i in 1:length(sheet_ids)){
  fish = read_excel("C:/Users/brian.grieve/Documents/coca/statedata/Delaware_all.xlsx", sheet = sheet_ids[i], 
                    col_types = c("numeric", "numeric", "numeric", "numeric", "numeric", "text", "numeric", "numeric", "numeric", "numeric",
                    "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))
  fish$SVSPP = del_svspp[i];
  if (i==1){
    delaware = fish;}
  if (i!=1){
    delaware = rbind(delaware,fish)}
}

delaware$LAT = rowMeans(cbind(delaware$BegLat,delaware$EndLat),na.rm=TRUE);
delaware$LON = rowMeans(cbind(delaware$BegLon,delaware$EndLon),na.rm=TRUE)*-1;
delaware = delaware[-which(is.na(delaware$LAT+delaware$LON)),]
delaware[which(is.na(delaware$`Tow Length (NM)`)),'Tow Length (NM)'] = 1
delaware$AREA_SWEPT_WINGS_MEAN_KM2 = (42*.55*12/39)*delaware$`Tow Length (NM)`*1852/1E6 # Footrope feet to meters * constant * tow length NM to km2 
names(delaware)[c(2:6,11,13,18,20,21)] = c('EST_MONTH','EST_DAY','EST_YEAR','STATION','EST_TIME','DOPDISTB','BTEMP','AVGDEPTH','num','wt')
delaware$Department = 'Delaware'
i = which(str_count(delaware$EST_TIME)==3);
delaware$EST_TIME[i] = paste('0',delaware$EST_TIME[i],sep='');
delaware$EST_TIME[(which(str_count(delaware$EST_TIME)!=4))] = '1200';
delaware$EST_TIME = paste(str_sub(delaware$EST_TIME,1,2),':',str_sub(delaware$EST_TIME,3,4),':00',sep='');
delaware$EST_TIME = chron(times. = delaware$EST_TIME,format=c(times='h:m:s'))
delaware$SEASON = NA;
delaware$date = NA;
datestring = paste(delaware$EST_YEAR,'-',delaware$EST_MONTH,'-',delaware$EST_DAY,sep='')
i = which(delaware$EST_DAY<28)
delaware$date[i] = chron(dates.=datestring[i],format='y-m-d')
i = which(delaware$EST_DAY>=28)
delaware$date[i] = chron(dates.=paste(delaware$EST_YEAR[i],'-',delaware$EST_MONTH[i],'-',delaware$EST_DAY[i]-15,sep=''),format='y-m-d')
delaware$date[i] = delaware$date[i]+15;
delaware$EST_DAY = as.numeric(days(delaware$date))
esttime = chron(times.=delaware$EST_TIME,format='h:m:s');
delaware$time.num = hours(esttime)+minutes(esttime)/60
zenith.df = astrocalc4r(day=delaware$EST_DAY,month=delaware$EST_MONTH,year=delaware$EST_YEAR,hour=delaware$time.num,timezone=c(rep(-5,times=nrow(delaware))),lat=delaware$LAT,lon=delaware$LON)
delaware$zenith = zenith.df$zenith
delaware$night = delaware$time.num>zenith.df$sunset | delaware$time.num<zenith.df$sunrise


# Maryland
spdf = data.frame(matrix(data=NA,nrow=12,ncol=2))
names(spdf) = c('common','svspp')
spdf$common = c('Alewife','Atl_herring','Atl_menhaden',
                'Black_sea_bass','Blueback_herring','Pollock','Silver_hake','Striped_bass','Summer_flounder',
                'Tautog','Windowpane','Winter_flounder')
spdf$svspp = c(33,32,36,141,34,75,72,139,103,177,108,106)
for (i in 1:nrow(spdf)){
  fish = read.csv(paste("C:/Users/brian.grieve/Documents/coca/statedata/Maryland/time-mod/",spdf[i,'common'],".csv",sep=''),header=T)
  fish$SVSPP = spdf[i,'svspp']
  if (i==1){
    maryland = fish}
  if (i!=1){
    maryland = rbind(maryland,fish)}
}
maryland[which(is.na(maryland$SpeciesCount)),'SpeciesCount']=0;  
i = which(maryland$ProblematicSampleYN=='Yes')
maryland = maryland[-i,];
maryland$AreaCovered = maryland$AreaCovered*144/39.37^2/1000000 # feet2 to in2 to m2 to km2
maryland$AVGDEPTH = rowMeans(cbind(maryland$Depthstart,maryland$DepthStop),na.rm=TRUE);
# THIS CHANGES STEMP TO BTEMP. Shallow waters, .984 correlation in common samples. Put in paper. 
names(maryland)[c(2,3,4,6,13,14,27,38)] = c('EST_YEAR','EST_MONTH','EST_DAY','EST_TIME','TOWDUR','AREA_SWEPT_WINGS_MEAN_KM2','BTEMP','num')
maryland$Department = 'Maryland'
maryland$EST_TIME = chron(times.=maryland$EST_TIME,format = c(times='h:m:s'))

# Lat/Lon from station data
station_data = read_excel("C:/Users/brian.grieve/Documents/coca/statedata/Maryland/station_data.xlsx", sheet = "Sheet2")
maryland[,c('LON','LAT')] = NA
u = unique(maryland$SiteNumber)
for (i in 1:length(u)){
  site = which(maryland$SiteNumber==u[i]);
  maryland[site,'LON'] = station_data[[i,'LON']]
  maryland[site,'LAT'] = station_data[[i,'LAT']]
}
maryland = maryland[-which(is.na(maryland$LON)),]
maryland$SEASON = NA
esttime = chron(times.=maryland$EST_TIME,format='h:m:s');
maryland$time.num = hours(esttime)+minutes(esttime)/60
maryland[which(maryland$time.num==0),'time.num'] = 12;
zenith.df = astrocalc4r(day=maryland$EST_DAY,month=maryland$EST_MONTH,year=maryland$EST_YEAR,hour=maryland$time.num,timezone=c(rep(-5,times=nrow(maryland))),lat=maryland$LAT,lon=maryland$LON)
maryland$zenith = zenith.df$zenith
maryland$night = maryland$time.num>zenith.df$sunset | maryland$time.num<zenith.df$sunrise


# neamap
neamap = read.csv("C:/Users/brian.grieve/Documents/coca/statedata/NEAMAP Catch Data_Bell.csv")
n = nrow(neamap); 
spdf = data.frame(matrix(NA,15,2)); names(spdf) = c('common','SVSPP'); 
spdf$common = names(neamap)[11:25];
spdf$SVSPP = c(33,36,34,301,26,193,141,103,15,139,177,106,108,23,105)
neamap2 = data.frame(matrix(data=NA,nrow = 15*nrow(neamap),ncol=15));
neamap2[,1:3] = neamap[,1:3];
neamap2[,7:13] = neamap[,4:10]
names(neamap2) = c('CRUISE6','STATION','DATE','EST_YEAR','EST_MONTH','EST_DAY','EST_TIME','SEASON','AREA_SWEPT_WINGS_MEAN_KM2','LAT','LON','AVGDEPTH','BTEMP','SVSPP','num')
neamap2$DATE = chron(dates.=as.character(neamap2$DATE), format = c(dates='dd-mon-yy'))
i = which(str_count(neamap2$EST_TIME)==4);
neamap2$EST_TIME = as.character(neamap2$EST_TIME)
neamap2$EST_TIME[i] = paste('0',neamap2$EST_TIME[i],sep='');
neamap2$EST_TIME = paste(str_sub(neamap2$EST_TIME,1,5),':00',sep='');
neamap2$EST_TIME = chron(times.=neamap2$EST_TIME, format = c(times='h:m:s'))
neamap2$EST_DAY = as.numeric(days(neamap2$DATE)); 
neamap2$EST_YEAR = as.numeric(as.character(years(neamap2$DATE)))
neamap2$EST_MONTH = as.numeric(months(neamap2$DATE))
neamap2$SEASON = 0;
neamap2$SEASON[which(neamap2$EST_MONTH>=3 & neamap2$EST_MONTH<=6)] = 2
neamap2$SEASON[which(neamap2$EST_MONTH>=9 & neamap2$EST_MONTH<=11)] = 4
neamap2$AREA_SWEPT_WINGS_MEAN_KM2 = neamap2$AREA_SWEPT_WINGS_MEAN_KM2/1E6
for (i in 1:nrow(spdf)){
  neamap2[((i-1)*n+1):(i*n),'SVSPP'] = spdf$SVSPP[i];
  neamap2[((i-1)*n+1):(i*n),'num'] = neamap[,spdf$common[i]]}
neamap2$Department = 'neamap'
esttime = chron(times.=neamap2$EST_TIME,format='h:m:s');
neamap2$time.num = hours(esttime)+minutes(esttime)/60
zenith.df = astrocalc4r(day=neamap2$EST_DAY,month=neamap2$EST_MONTH,year=neamap2$EST_YEAR,hour=neamap2$time.num,timezone=c(rep(-5,times=nrow(neamap2))),lat=neamap2$LAT,lon=neamap2$LON)
neamap2$zenith = zenith.df$zenith
neamap2$night = neamap2$time.num>zenith.df$sunset | neamap2$time.num<zenith.df$sunrise

#RI
load("C:/Users/brian.grieve/Documents/coca/statedata/RIDEM Seasonal Trawl_4172019.RData")
i = which(is.na(abundance$nBottomTemp)|abundance$nBottomTemp==0)
abundance[i,'nBottomTemp'] = abundance[i,'nSurfaceTemp']; # shallow, .968 correlation.  
i = which(is.na(abundance$nBottomTemp)|is.na(abundance$nStartDepth)|is.na(abundance$LAT))
abundance = abundance[-i,]
biomass = biomass[-i,]
n = nrow(abundance); 
spdf = data.frame(matrix(NA,23,2)); names(spdf) = c('common','SVSPP'); 
spdf$common = names(abundance)[12:34];
spdf$SVSPP = c(73,74,105,106,103,107,102,32,75,76,36,139,15,141,23,26,108,193,72,177,301,33,34)
ridata = data.frame(matrix(data=NA,nrow = 23*n,ncol=14));
ridata[,1:11] = rep(abundance[,1:11],times=23);
names(ridata) = c('CRUISE6','Tow','LAT','LON','EST_YEAR','EST_MONTH','EST_DAY','AVGDEPTH','STEMP','BTEMP','AREA_SWEPT_WINGS_MEAN_KM2','SVSPP','num','wt')
for (i in 1:nrow(spdf)){
  ridata[((i-1)*n+1):(i*n),'SVSPP'] = spdf$SVSPP[i];
  ridata[((i-1)*n+1):(i*n),'num'] = abundance[,spdf$common[i]]
  ridata[((i-1)*n+1):(i*n),'wt'] = biomass[,spdf$common[i]]}
ridata[which(is.na(ridata$num)),'num'] = 0;
ridata[which(is.na(ridata$wt)),'wt'] = 0;
ridata$Department = 'RI'
ridata$AREA_SWEPT_WINGS_MEAN_KM2 = ridata$AREA_SWEPT_WINGS_MEAN_KM2/1E6;
ridata$SEASON = NA;


# Combine all data sources
alltrawl = trawl.corrected; 
alltrawl$Department = 'NOAA'
alltrawl$EST_TIME = chron(times.=alltrawl$EST_TIME,format=c(times.='h:m:s'))
names_noaa = names(alltrawl);
trawlsurveys = c('MDMF','conn','nj2','nj.o','menh','delaware','maryland','neamap2','ridata')
for (tss in trawlsurveys){
  ts_data = eval(parse(text=tss));
  names_ts = names(ts_data);
  ts_data_rearranged = data.frame(matrix(data=NA,nrow=nrow(ts_data),ncol=ncol(alltrawl))); names(ts_data_rearranged)=names_noaa;
  ind_ts = which(is.element(names_ts,names_noaa))
  ind_noaa = which(is.element(names_noaa,names_ts))
  ts_data_rearranged[,sort(names_noaa[ind_noaa])] = ts_data[,sort(names_ts[ind_ts])];
  alltrawl = rbind(alltrawl,ts_data_rearranged)
}
alltrawl[which(alltrawl$EST_MONTH>=1 & alltrawl$EST_MONTH<=3 & is.na(alltrawl$SEASON)),'SEASON']=1
alltrawl[which(alltrawl$EST_MONTH==12 & is.na(alltrawl$SEASON)),'SEASON']=1
alltrawl[which(alltrawl$EST_MONTH>=4 & alltrawl$EST_MONTH<=6 & is.na(alltrawl$SEASON)),'SEASON']=2
alltrawl[which(alltrawl$EST_MONTH>=7 & alltrawl$EST_MONTH<=8 & is.na(alltrawl$SEASON)),'SEASON']=3
alltrawl[which(alltrawl$EST_MONTH>=9 & alltrawl$EST_MONTH<=11 & is.na(alltrawl$SEASON)),'SEASON']=4
i = which(alltrawl$EST_MONTH>=10 & alltrawl$EST_DAY>=10);
alltrawl$date_id[i] = as.numeric(paste(alltrawl$EST_YEAR[i],alltrawl$EST_MONTH[i],alltrawl$EST_DAY[i],sep=''))
i = which(alltrawl$EST_MONTH<10 & alltrawl$EST_DAY>=10)
alltrawl$date_id[i] = as.numeric(paste(alltrawl$EST_YEAR[i],0,alltrawl$EST_MONTH[i],alltrawl$EST_DAY[i],sep=''))
i = which(alltrawl$EST_MONTH>=10 & alltrawl$EST_DAY<10)
alltrawl$date_id[i] = as.numeric(paste(alltrawl$EST_YEAR[i],alltrawl$EST_MONTH[i],0,alltrawl$EST_DAY[i],sep=''))
i = which(alltrawl$EST_MONTH<10 & alltrawl$EST_DAY<10)
alltrawl$date_id[i] = as.numeric(paste(alltrawl$EST_YEAR[i],0,alltrawl$EST_MONTH[i],0,alltrawl$EST_DAY[i],sep=''))
alltrawl$DATE = chron(dates.=paste(str_sub(alltrawl$date_id,1,4),'/',str_sub(alltrawl$date_id,5,6),'/',str_sub(alltrawl$date_id,7,8),sep=''),format=c('y/m/d'))
alltrawl$EST_YEAR = as.numeric(alltrawl$EST_YEAR);
alltrawl[which(is.na(alltrawl$night)),'night'] = FALSE;

i = which(alltrawl$Department!='NOAA')
alltrawl[i,'num_raw'] = alltrawl[i,'num']
alltrawl[i,'wt_raw'] = alltrawl[i,'wt']
alltrawl[i,'presence'] = alltrawl[i,'num_raw']>0

# Swept Area Conversion 
SA_alb = getmode(alltrawl$AREA_SWEPT_WINGS_MEAN_KM2[which(alltrawl$SV_VESSEL=='AL')]); # Swept Area to Albatross units
alltrawl$SA_rat = SA_alb/alltrawl$AREA_SWEPT_WINGS_MEAN_KM2; alltrawl$SA_rat[which(alltrawl$SV_VESSEL=='HB')]=1
alltrawl$num_sa = alltrawl$num*SA_rat

# Correct for day/night catches
alltrawl$num_sa_night = alltrawl$num_sa;
daynight = data.frame(matrix(data=NA,nrow=nrow(species_conversions),ncol=6)); names(daynight) = c('CommonName','SVSPP','daymean','nightmean','pvalue','ratio')
daynight$CommonName = species_conversions$CommonName; daynight$SVSPP = species_conversions$SVSPP;
for (spp in u_svspp){
  ind = which(alltrawl$SVSPP==spp & alltrawl$Department=='NOAA')
  fishdata = alltrawl[ind,];
  boxplot(formula(fishdata$num_sa~fishdata$night))
  
  # If significantly different, divide by nighttime catch to correct to daytime
  tt = t.test(formula(fishdata$num_sa~fishdata$night))
  ind_dn = which(daynight$SVSPP==spp)
  daynight[ind_dn,'pvalue'] = tt$p.value;
  daynight[ind_dn,'daymean'] = round(as.numeric(tt$estimate['mean in group FALSE']),2);
  daynight[ind_dn,'nightmean'] = round(as.numeric(tt$estimate['mean in group TRUE']),2);
  daynight[ind_dn,'ratio'] = round(daynight[ind_dn,'daymean']/daynight[ind_dn,'nightmean'],3)
  if (tt$p.value<0.05){
    i_n = which(alltrawl$SVSPP==i & alltrawl$night==1);
    alltrawl[i_n,'num_sa_night'] = alltrawl[i_n,'num_sa']*daynight[ind_dn,'ratio'];}
}


# Get STRATUM area and count how many times each one was sampled per cruise
stratarea = read.csv('C:/Users/brian.grieve/Documents/coca/Data/strata_area.csv')
stratarea$Area = stratarea$Area/.29; # naut.mile2 to km2
u_strat = unique(alltrawl$STRATUM)
for (i in u_strat){
  ind = which(alltrawl$STRATUM==i);
  ind2 = which(stratarea$Strata==i);
  if(length(ind2)!=0){
    alltrawl[ind,'stratarea'] = stratarea[ind2,'Area'];}
}
alltrawl[which(alltrawl$STRATUM==1351),'stratarea'] = stratarea[which(stratarea$Strata==1350),'Area'];

alltrawl$stratpercruise = NA
u_cruiseID = unique(alltrawl[which(alltrawl$Department=='NOAA'),'CRUISE6'])
for (i in u_cruiseID){
  ind = which(alltrawl$CRUISE6==i & alltrawl$SVSPP==106 & alltrawl$Department=='NOAA');
  cruisedata = alltrawl[ind,];
  u_strata = unique(cruisedata[,'STRATUM']); 
  for (strat in u_strata){
    ind_strat = which(alltrawl$STRATUM==strat & alltrawl$CRUISE6==i & alltrawl$Department=='NOAA');
    alltrawl[ind_strat, 'stratpercruise'] = length(which(cruisedata$STRATUM==strat));}}



# Partition Data Set
for(spp in u_svspp){
  for(seas in c(4)){
    i_spp = which(alltrawl$SVSPP==spp & alltrawl$SEASON==seas);
    fishdata = alltrawl[i_spp,];
    
    
    # Partition dataset
    remainder = nrow(fishdata)%%5;
    plength = floor(nrow(fishdata)/5)
    if (remainder!=0){
      set.seed(seed1)
      partition_inds = sample(rep(1:5,plength),replace=FALSE)
      partition_inds[seq(plength*5+1,plength*5+remainder,1)]=1:remainder
    } else {
      set.seed(seed1)
      partition_inds = sample(rep(1:5,plength),replace=FALSE)}
    fishdata$partition = partition_inds;
    alltrawl[i_spp,'partition'] = partition_inds;
  }
}
    
