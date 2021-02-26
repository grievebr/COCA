# Evaluate ROMS COBALT with trawl observations

library(ncdf4); library(stringr); library(RANN); library(rgdal); library(gstat); library(lattice); library(plot3D); library(colorspace)

load("C:/Users/brian.grieve/Documents/coca/Data/alltrawlnn.RData")
surveydata=alltrawl[which(alltrawl$date_id>=19800101 & alltrawl$date_id<= 20111214 & alltrawl$SVSPP==108),] # 108 has most observations. 


dates = seq.Date(from = as.Date('1980/02/01'),to = as.Date('2014/12/08'), by = 'day');
modelinfo = data.frame(matrix(data = NA, nrow = length(dates), ncol = 7)); names(modelinfo) = c('ind','date','dateid','year','month','day','modjday')
modelinfo$ind = 1:length(dates); 
modelinfo$date = dates;
modelinfo$dateid = paste(str_sub(dates,1,4),str_sub(dates,6,7),str_sub(dates,9,10),sep='');
modelinfo$year = as.numeric(str_sub(dates,1,4))
modelinfo$month = as.numeric(str_sub(dates,6,7))
modelinfo$day = as.numeric(str_sub(dates,9,10))
for (i in 1980:2011){
  i_year = which(modelinfo$year==i)
  modelinfo[i_year,'modjday'] = julian(dates[i_year],origin = modelinfo[i_year[1],'date'])+1}

# Compare daily ROMS BWT and trawl BWT of same day
filename = 'C:/Users/brian.grieve/Documents/Ocean_models/NWA_grd.nc';
nc = nc_open(filename); 
lon = ncvar_get(nc,'lon_rho'); 
lat = ncvar_get(nc,'lat_rho'); 
lon = lon[381:600,216:295]; # Subset of NWA ROMS
lat = lat[381:600,216:295];
m = nrow(lon);
n = ncol(lon);

# adjust for negative coordinates
latlon = cbind(as.vector(lat),as.vector(lon-360))

# Determine which surveys and models overlap 
survdays = which(is.element(modelinfo$dateid,surveydata$date_id))
survpoints = which(is.element(surveydata$date_id,modelinfo$dateid))

compmat = data.frame(matrix(data=NA,nrow=length(survpoints),ncol=15)); names(compmat) = c('fishlat','fishlon','BTEMP','modlat','modlon','modtemp','k','department','EST_YEAR','EST_MONTH','EST_DAY','survday','i_k','ind_x','ind_y')
compmat[,c('survday','EST_YEAR','EST_MONTH','EST_DAY','department','fishlat','fishlon','BTEMP')] = surveydata[survpoints,c('date_id','EST_YEAR','EST_MONTH','EST_DAY','Department','LAT','LON','BTEMP')];

for (i in survdays){
  if (modelinfo[i,'day']==1 & modelinfo[i,'month']==1){print(survdays[i,'dateid'])}
  mod.filename = paste('C:/Users/brian.grieve/Documents/Ocean_models/bwt_1980-2014/bwt_1980-2014/bwt_',modelinfo[i,'year'],'.nc',sep='')
  nctemp = nc_open(mod.filename)
  temp = ncvar_get(nctemp,'bwt');
  temp = temp[,,modelinfo[i,'modjday']];
  temp[which(temp>1E30)] = NA
  nc_close(nctemp);
  d = which(compmat$survday==modelinfo[i,'dateid'])
  for (i_d in d){
    # Determine ROMS cells of trawl coordinates
    nndf = nn2(latlon,query=compmat[i_d,c('fishlat','fishlon')],k=1,searchtype='standard');
    i_k = nndf[[1]]
    nnindist = nndf[[2]];
    ind_x = i_k%/%m+1
    ind_y = i_k%%m
    compmat[i_d,'i_k'] = i_k;
    compmat[i_d,'ind_x'] = ind_x;
    compmat[i_d,'ind_y'] = ind_y;
    compmat[i_d,'modlat'] = lat[ind_y,ind_x]
    compmat[i_d,'modlon'] = lon[ind_y,ind_x]-360;
    compmat[i_d,'modtemp'] = temp[ind_y,ind_x]
    compmat[i_d,'k'] = round(nnindist[which(nnindist==min(nnindist))],3)}}
compmat2 = compmat[-which(is.na(compmat$BTEMP)|is.na(compmat$modtemp)|compmat$BTEMP==0),]


# Determine bias by season
# Months
spring = 3:6;
fall = 8:10;

# Spring
ind_spring = is.element(compmat2[,'EST_MONTH'],spring);
trawltemp = matrix(NA,nrow=nrow(lat),ncol=ncol(lat));
romstemp = matrix(NA,nrow=nrow(lat),ncol=ncol(lat));
compmat_spring = compmat2[ind_spring,]; 
u = unique(compmat_spring$i_k); 
for (i in u){
  ind = which(compmat_spring$i_k==i);
  ind_x = compmat_spring[ind[1],'ind_x']
  ind_y = compmat_spring[ind[1],'ind_y']
  trawltemp[ind_y,ind_x] = mean(compmat_spring[ind,'BTEMP'],na.rm=T)
  romstemp[ind_y,ind_x] = mean(compmat_spring[ind,'modtemp'],na.rm=T)}
tempbiasspring = trawltemp-romstemp;

image2D(x=lon-360,y=lat,z=tempbiasspring,zlim=c(-10,10),xlim=c(-76,-65),ylim=c(35,46),main='Spring modelbias',col=palette(diverge_hsv(20)))
#write.table(tempbiasspring,'C:/Users/brian.grieve/Documents/coca/Data/tempbiasspring.txt',row.names = FALSE ,col.names=FALSE)

# Fall
ind_fall = is.element(compmat2[,'EST_MONTH'],fall);
trawltemp = matrix(NA,nrow=nrow(lat),ncol=ncol(lat));
romstemp = matrix(NA,nrow=nrow(lat),ncol=ncol(lat));
compmat_fall = compmat2[ind_fall,]; 
u = unique(compmat_fall$i_k); 
for (i in u){
  ind = which(compmat_fall$i_k==i);
  ind_x = compmat_fall[ind[1],'ind_x']
  ind_y = compmat_fall[ind[1],'ind_y']
  trawltemp[ind_y,ind_x] = mean(compmat_fall[ind,'BTEMP'],na.rm=T)
  romstemp[ind_y,ind_x] = mean(compmat_fall[ind,'modtemp'],na.rm=T)}
tempbiasfall = trawltemp-romstemp;

image2D(x=lon-360,y=lat,z=-1*tempbiasfall,zlim=c(-10,10),xlim=c(-76,-65),ylim=c(35,46),main='Fall modelbias',col=palette(diverge_hsv(20)))
#rite.table(tempbiasfall,'C:/Users/brian.grieve/Documents/coca/Data/tempbiasfall.txt',row.names = FALSE ,col.names=FALSE)


# Map annual bias
for (y in 1980:2011){
  ind_year = which(compmat2$EST_YEAR==y); #print(length(ind_year))
  compmat_yr = compmat2[ind_year,]
  trawltemp = matrix(NA,nrow=nrow(lat),ncol=ncol(lat));
  romstemp = matrix(NA,nrow=nrow(lat),ncol=ncol(lat));
  u = unique(compmat_yr$i_k); 
  for (i in u){
    ind = which(compmat_yr$i_k==i);
    ind_x = compmat_yr[ind[1],'ind_x']
    ind_y = compmat_yr[ind[1],'ind_y']
    trawltemp[ind_y,ind_x] = mean(compmat_yr[ind,'BTEMP'],na.rm=T)
    romstemp[ind_y,ind_x] = mean(compmat_yr[ind,'modtemp'],na.rm=T)}
  tempbias = romstemp-trawltemp;
  image2D(x=lon-360,y=lat,z=1*tempbias,zlim=c(-10,10),xlim=c(-76,-65),ylim=c(35,46),main=y,col=palette(diverge_hsv(20)))}
  

