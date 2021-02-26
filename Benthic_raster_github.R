# Convert benthic rasters to matrices

library(rgdal); library(raster); library(sp); library(RANN); library(sf)
setwd('C:/Users/brian.grieve/Documents/Marta_benthicdata/RASTERS_100_LATLONG')

# Open GeoTiffs and convert to numeric matrices

bathyclasses.sp = readGDAL('C:/Users/brian.grieve/Documents/Marta_benthicdata/SPATIALDATA_UPDATE_03_2019/bathy_classes.tif');
bathyclasses.mat = as.matrix(bathyclasses.sp);
bathyclasses.coord = coordinates(bathyclasses.sp);
bathyclasses.df = data.frame(cbind(bathyclasses.coord,as.vector(bathyclasses.mat))); names(bathyclasses.df) = c('x','y','bathyclasses')
write.csv(bathyclasses.df,'C:/Users/brian.grieve/Documents/Marta_benthicdata/bathyclassesdf.csv',row.names=FALSE,col.names=TRUE);
rm(bathyclasses.sp,bathyclasses.mat,bathyclasses.coord,bathyclasses.df)

bpi.broad.classes.sp = readGDAL('C:/Users/brian.grieve/Documents/Marta_benthicdata/SPATIALDATA_UPDATE_03_2019/bpi_broad_classes.tif'); # 0=flat; 1000=valley; 2000=peaks
bpi.broad.classes.mat = as.matrix(bpi.broad.classes.sp);
bpi.broad.classes.coord = coordinates(bpi.broad.classes.sp);
bpi.broad.classes.df = data.frame(cbind(bpi.broad.classes.coord,as.vector(bpi.broad.classes.mat))); names(bpi.broad.classes.df) = c('x','y','bpi.broad.classes')
write.csv(bpi.broad.classes.df,'C:/Users/brian.grieve/Documents/Marta_benthicdata/bpibroadclassesdf.csv',row.names=FALSE);
rm(bpi.broad.classes.sp,bpi.broad.classes.mat,bpi.broad.classes.coord,bpi.broad.classes.df)

bpi.fine.classes.sp = readGDAL('C:/Users/brian.grieve/Documents/Marta_benthicdata/SPATIALDATA_UPDATE_03_2019/bpi_fine_classes.tif');
bpi.fine.classes.mat = as.matrix(bpi.fine.classes.sp);
bpi.fine.classes.coord = coordinates(bpi.fine.classes.sp);
bpi.fine.classes.df = data.frame(cbind(bpi.fine.classes.coord,as.vector(bpi.fine.classes.mat))); names(bpi.fine.classes.df) = c('x','y','bpi.fine.classes')
write.csv(bpi.fine.classes.df,'C:/Users/brian.grieve/Documents/Marta_benthicdata/bpifineclassesdf.csv',row.names=FALSE,col.names=TRUE);
rm(bpi.fine.classes.sp,bpi.fine.classes.mat,bpi.fine.classes.coord,bpi.fine.classes.df)

sedclasses.sp = readGDAL('C:/Users/brian.grieve/Documents/Marta_benthicdata/SPATIALDATA_UPDATE_03_2019/sediment_classes.tif');
sedclasses.mat = as.matrix(sedclasses.sp);
sedclasses.coord = coordinates(sedclasses.sp);
sedclasses.df = data.frame(cbind(sedclasses.coord,as.vector(sedclasses.mat))); names(sedclasses.df) = c('x','y','sedclasses')
write.csv(sedclasses.df,'C:/Users/brian.grieve/Documents/Marta_benthicdata/sedclassesdf.csv',row.names=FALSE,col.names=TRUE);
rm(sedclasses.sp,sedclasses.mat,sedclasses.coord,sedclasses.df)

HB.sp = readGDAL('C:/Users/brian.grieve/Documents/Marta_benthicdata/DATA_UPDATE_02_2020/HB_possible.tif');
HB.mat = as.matrix(HB.sp);
HB.coord = coordinates(HB.sp);
HB.df = data.frame(cbind(HB.coord,as.vector(HB.mat))); names(HB.df) = c('x','y','HB')
write.csv(HB.df,'C:/Users/brian.grieve/Documents/Marta_benthicdata/HB.csv',row.names=FALSE,col.names=TRUE);
rm(HB.sp,HB.mat,HB.coord,HB.df)

# Determine locations of benthic classifications
setwd('C:/Users/brian.grieve/Documents/Marta_benthicdata')
dflist = c('sedclasses','bpibroadclasses','bathyclasses','sumbpiFbathsed','sumbpiBbathsed','bpifineclasses','HB');
n = ncol(alltrawl)
nandistance.df = data.frame(matrix(data=NA,nrow=length(dflist),ncol=3)); names(nandistance.df) = c('Variable','dist.num','finite.num')
for (i in 1:length(dflist)){
  print(i)
  var.df = read.csv(paste('C:/Users/brian.grieve/Documents/Marta_benthicdata/',dflist[i],'df.csv',sep=''));
  ind.nn = nn2(var.df[,1:2],cbind(alltrawl$LON,alltrawl$LAT),k=1);
  na.dist = which(ind.nn$nn.dist>.25); 
  nandistance.df[i,] = cbind(dflist[i],length(na.dist),length(which(is.finite(var.df[ind.nn$nn.idx[na.dist],3])==1)))
  var.df[ind.nn$nn.idx[na.dist],3] = NA;
  alltrawl[,n+i] = var.df[ind.nn$nn.idx,3];}
names(alltrawl) = c(names(alltrawl)[1:n],dflist)


# Station locations for patch dynamics
stationpoints = st_read('C:/Users/brian.grieve/Documents/Marta_benthicdata/Station_points_with_distance_classes_UPDATED/Station_points_with_distance_classes_UPDATED.shp')
station.df = as.data.frame(stationpoints)
station.geom = st_geometry(stationpoints)
features = st_read('C:/Users/brian.grieve/Documents/Marta_benthicdata/DATA_UPDATE_02_2020/Combined_features_Compactness.shp')
features.geom = st_geometry(features)

u_latlon = (cbind(station.df$V1, station.df$V2))
u_Id = unique(features$Id)
alltrawl.features = data.frame(matrix(data=NA,nrow=nrow(u_latlon),ncol=8)); 
names(alltrawl.features) = c('i','LAT','LON','Id','gridcode','Area','Perim','COMPACT');
alltrawl.features$LAT = u_latlon[,1];
alltrawl.features$LON = u_latlon[,2];

for(i in 1:length(u_Id)){
  if(is.element(i,seq(1,length(u_Id),100))){print(i)}
  points = which(point.in.polygon(u_latlon[,2],u_latlon[,1],features.geom[[i]][[1]][,1],features.geom[[i]][[1]][,2])==1)
  if(length(points)>0){
    alltrawl.features[points,'i'] = i;
    alltrawl.features[points,'Id'] = features[[i,'Id']];
    alltrawl.features[points,'gridcode'] = features[[i,'gridcode']];
    alltrawl.features[points,'Area'] = features[[i,'Area']];
    alltrawl.features[points,'Perim'] = features[[i,'Perim']];
    alltrawl.features[points,'COMPACT'] = features[[i,'COMPACT']];}}
alltrawl.features = data.frame(alltrawl.features,NEAR_FID = station.df[,'NEAR_FID'], NEAR_DIST = station.df[,'NEAR_DIST'], CODE_ON = station.df[,'CODE_ON'], CODE_NEA = station.df[,'CODE_NEA']); 



alltrawl$PatchID = NA; alltrawl$gridcode = NA; alltrawl$PatchArea = NA; alltrawl$PatchPerim = NA; alltrawl$PatchCompact = NA; alltrawl$PatchNearID = NA; alltrawl$PatchNearDist = NA; alltrawl$CODE_ON = NA; alltrawl$CODE_NEA = NA; 
for (i_latlon in 1:nrow(u_latlon)){
  if(is.element(i_latlon,seq(1,nrow(u_latlon),100))){print(i_latlon)}
  xq = u_latlon[i_latlon,2];
  yq = u_latlon[i_latlon,1];
  i = which(alltrawl$LAT==yq & alltrawl$LON==xq);
  if(length(i)>0){
    alltrawl[i,c('PatchID','gridcode','PatchArea','PatchPerim','PatchCompact','PatchNearID','PatchNearDist','CODE_ON','CODE_NEA')] = alltrawl.features[i_latlon,4:12]}}


# Create XY Grids of shapefiles
features = st_read('C:/Users/brian.grieve/Documents/Marta_benthicdata/DATA_UPDATE_02_2020/Combined_features_Compactness.shp')
features.geom = st_geometry(features)
sumbpiB.bath.sed.sp = readGDAL('C:/Users/brian.grieve/Documents/Marta_benthicdata/SPATIALDATA_UPDATE_03_2019/sum_bpbroad_bathy_sedim.tif');
sumbpiB.bath.sed.coord = coordinates(sumbpiB.bath.sed.sp);

lonvec = unique(sumbpiB.bath.sed.coord[,1])
latvec = unique(sumbpiB.bath.sed.coord[,2])
latgrid = matrix(data=latvec,nrow=10704,ncol=9869,byrow=T)
longrid = matrix(data=lonvec,nrow=10704,ncol=9869,byrow=F)

# Obtain patch information for each pixel
patchpixel = expand.grid(lonvec,latvec); names(patchpixel) = c('lon','lat')
patchpixel$ind = 1:nrow(patchpixel)
patchpixel$Id = NA;
patchpixel$perim = NA
patchpixel$area = NA;
patchpixel$compact = NA;
patchpixel = patchpixel[-which(is.na(sumbpiB.bath.sed.df$sumbpiB.bath.sed)),]
                               
u.Id = unique(features$Id)
for (i.id in 1:length(u.Id)){
  if(is.element(i.id,seq(1,length(u.Id),100))){print(i.id)}
  points = which(point.in.polygon(patchpixel[,1],patchpixel[,2],features.geom[[i.id]][[1]][,1],features.geom[[i.id]][[1]][,2])==1)
  if(length(points)>0){
  patchpixel[points,'perim'] = features$Perim[i.id]
  patchpixel[points,'area'] = features$Area[i.id]
  patchpixel[points,'Id'] = features$Id[i.id]
  patchpixel[points,'compact'] = features$COMPACT[i.id]
}}

patchmat = longrid2; patchmat[,] = NA;
patchmat[patchpixel$ind] = patchpixel$compact
compactmat = patchmat; 
patchmat[patchpixel$ind] = patchpixel$perim
perimmat = patchmat; 
patchmat[patchpixel$ind] = patchpixel$area
areamat = patchmat
write.csv(patchmat,file='C:/Users/brian.grieve/Documents/coca/Data/patchmat.csv',row.names = FALSE,col.names=FALSE)
write.csv(areamat,file='C:/Users/brian.grieve/Documents/coca/Data/areamat.csv',row.names = FALSE,col.names=FALSE)
write.csv(compactmat,file='C:/Users/brian.grieve/Documents/coca/Data/compactmat.csv',row.names = FALSE,col.names=FALSE)


