# This script runs through all possible combinations of seven habitat variables to determine the best habitat model


library(dplyr); library(mgcv); library(bbmle); library(Metrics); library(caret); library(readxl); library(rgdal); library(gstat); library(lattice); library(PresenceAbsence)


load("C:/Users/brian.grieve/Documents/coca/Data/alltrawlnn.RData")
load('C:/Users/brian.grieve/Documents/coca/Data/habitatvars2.RData')
species_conversions = read_excel("C:/Users/brian.grieve/Documents/coca/BellR_species_listConversions.xlsx")
u_svspp = unique(alltrawl$SVSPP)

# Functions for the temperature portion. The final model is trained using TMB results, but preliminary results using this 
# method are still used here for some fit statistics/model selection
B_Afunction<-function(T.c, Topt.c, Er, Ed, cc){
  bk<-8.62*10^-5   #k is Boltzmann constant 
  T.k<-T.c+273 # converts celcius to kelvin  
  Topt.k<-Topt.c+273 # converts celcius to kelvin  
  P<-exp(-Er/(bk*T.k))/(1+exp((-1/(bk*T.k)*(Ed-((Ed/Topt.k)+bk*log(Er/(Ed-Er)))*T.k))))*(cc*10^ccexp)
  return(P)
} 
#Log function for maximum likelihood 
B_AfunctionNLL2_nbin<-function(T.c, Topt.c,Er, Ed, cc, k){ 
  meanresponse<-B_Afunction(T.c, Topt.c,Er, Ed, cc)
  -sum(dnbinom(survdat$abundance, mu = meanresponse, size=k, log = TRUE))  
}



# list of all models fit. Model selection was determined by lowest AIC score            
gflist = list()
{
gflist$gf1 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num))
gflist$gf2 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num))
gflist$gf3 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(PatchCompact,k=k.num))
gflist$gf4 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num)) 
gflist$gf5 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num))
gflist$gf6 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses)) 
gflist$gf7 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(log(Patch),k=k.num)) 
gflist$gf8 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(PatchCompact,k=k.num))
gflist$gf9 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + s(log(Patch),k=k.num))
gflist$gf10 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + s(PatchCompact,k=k.num))
gflist$gf11 = formula(presence ~ factor(bathyclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num))
gflist$gf12 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num))
gflist$gf13 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(PatchCompact,k=k.num))
gflist$gf14 = formula(presence ~ factor(bpibroadclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num)) 
gflist$gf15 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses))
gflist$gf16 = formula(presence ~ factor(bathyclasses) + factor(sedclasses))
gflist$gf17 = formula(presence ~ factor(bathyclasses) + s(log(Patch),k=k.num))
gflist$gf18 = formula(presence ~ factor(bathyclasses) + s(PatchCompact,k=k.num))
gflist$gf19 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses)) 
gflist$gf20 = formula(presence ~ factor(bpibroadclasses) + s(log(Patch),k=k.num)) 
gflist$gf21 = formula(presence ~ factor(bpibroadclasses) + s(PatchCompact,k=k.num))
gflist$gf22 = formula(presence ~ factor(sedclasses) + s(log(Patch),k=k.num))
gflist$gf23 = formula(presence ~ factor(sedclasses) + s(PatchCompact,k=k.num))
gflist$gf24 = formula(presence ~ s(log(Patch),k=k.num) + s(PatchCompact,k=k.num))
gflist$gf25 = formula(presence ~ factor(bathyclasses))
gflist$gf26 = formula(presence ~ factor(bpibroadclasses))
gflist$gf27 = formula(presence ~ factor(sedclasses))
gflist$gf28 = formula(presence ~ s(log(Patch),k=k.num))
gflist$gf29 = formula(presence ~ s(PatchCompact,k=k.num))
gflist$gf30 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
gflist$gf31 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses))
gflist$gf32 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
gflist$gf33 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses)) 
gflist$gf34 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
gflist$gf35 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + factor(bpifineclasses)) 
gflist$gf36 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses)) 
gflist$gf37 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
gflist$gf38 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses))
gflist$gf39 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
gflist$gf40 = formula(presence ~ factor(bathyclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
gflist$gf41 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses))
gflist$gf42 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
gflist$gf43 = formula(presence ~ factor(bpibroadclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses)) 
gflist$gf44 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(bpifineclasses))
gflist$gf45 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + factor(bpifineclasses))
gflist$gf46 = formula(presence ~ factor(bathyclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses))
gflist$gf47 = formula(presence ~ factor(bathyclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
gflist$gf48 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + factor(bpifineclasses)) 
gflist$gf49 = formula(presence ~ factor(bpibroadclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses)) 
gflist$gf50 = formula(presence ~ factor(bpibroadclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
gflist$gf51 = formula(presence ~ factor(sedclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses))
gflist$gf52 = formula(presence ~ factor(sedclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
gflist$gf53 = formula(presence ~ s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
gflist$gf54 = formula(presence ~ factor(bathyclasses) + factor(bpifineclasses))
gflist$gf55 = formula(presence ~ factor(bpibroadclasses) + factor(bpifineclasses))
gflist$gf56 = formula(presence ~ factor(sedclasses) + factor(bpifineclasses))
gflist$gf57 = formula(presence ~ s(log(Patch),k=k.num) + factor(bpifineclasses))
gflist$gf58 = formula(presence ~ s(PatchCompact,k=k.num) + factor(bpifineclasses))
gflist$gf59 = formula(presence ~ factor(bpifineclasses))

gflist$gf60 = formula(presence ~ factor(HB))
gflist$gf61 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(HB))
gflist$gf62 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + factor(HB))
gflist$gf63 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(PatchCompact,k=k.num) + factor(HB))
gflist$gf64 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(HB)) 
gflist$gf65 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(HB))
gflist$gf66 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + factor(HB)) 
gflist$gf67 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(log(Patch),k=k.num) + factor(HB)) 
gflist$gf68 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(PatchCompact,k=k.num) + factor(HB))
gflist$gf69 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + factor(HB))
gflist$gf70 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + s(PatchCompact,k=k.num) + factor(HB))
gflist$gf71 = formula(presence ~ factor(bathyclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(HB))
gflist$gf72 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + factor(HB))
gflist$gf73 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(PatchCompact,k=k.num) + factor(HB))
gflist$gf74 = formula(presence ~ factor(bpibroadclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(HB)) 
gflist$gf75 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(HB))
gflist$gf76 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + factor(HB))
gflist$gf77 = formula(presence ~ factor(bathyclasses) + s(log(Patch),k=k.num) + factor(HB))
gflist$gf78 = formula(presence ~ factor(bathyclasses) + s(PatchCompact,k=k.num) + factor(HB))
gflist$gf79 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + factor(HB)) 
gflist$gf80 = formula(presence ~ factor(bpibroadclasses) + s(log(Patch),k=k.num) + factor(HB)) 
gflist$gf81 = formula(presence ~ factor(bpibroadclasses) + s(PatchCompact,k=k.num) + factor(HB))
gflist$gf82 = formula(presence ~ factor(sedclasses) + s(log(Patch),k=k.num) + factor(HB))
gflist$gf83 = formula(presence ~ factor(sedclasses) + s(PatchCompact,k=k.num) + factor(HB))
gflist$gf84 = formula(presence ~ s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(HB))
gflist$gf85 = formula(presence ~ factor(bathyclasses) + factor(HB))
gflist$gf86 = formula(presence ~ factor(bpibroadclasses) + factor(HB))
gflist$gf87 = formula(presence ~ factor(sedclasses) + factor(HB))
gflist$gf88 = formula(presence ~ s(log(Patch),k=k.num) + factor(HB))
gflist$gf89 = formula(presence ~ s(PatchCompact,k=k.num) + factor(HB))
gflist$gf90 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
gflist$gf91 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses) + factor(HB))
gflist$gf92 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
gflist$gf93 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB)) 
gflist$gf94 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
gflist$gf95 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + factor(bpifineclasses) + factor(HB)) 
gflist$gf96 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses) + factor(HB)) 
gflist$gf97 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
gflist$gf98 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses) + factor(HB))
gflist$gf99 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
gflist$gf100 = formula(presence ~ factor(bathyclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
gflist$gf101 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses) + factor(HB))
gflist$gf102 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
gflist$gf103 = formula(presence ~ factor(bpibroadclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB)) 
gflist$gf104 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(bpifineclasses) + factor(HB))
gflist$gf105 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + factor(bpifineclasses) + factor(HB))
gflist$gf106 = formula(presence ~ factor(bathyclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses) + factor(HB))
gflist$gf107 = formula(presence ~ factor(bathyclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
gflist$gf108 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + factor(bpifineclasses) + factor(HB)) 
gflist$gf109 = formula(presence ~ factor(bpibroadclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses) + factor(HB)) 
gflist$gf110 = formula(presence ~ factor(bpibroadclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
gflist$gf111 = formula(presence ~ factor(sedclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses) + factor(HB))
gflist$gf112 = formula(presence ~ factor(sedclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
gflist$gf113 = formula(presence ~ s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
gflist$gf114 = formula(presence ~ factor(bathyclasses) + factor(bpifineclasses) + factor(HB))
gflist$gf115 = formula(presence ~ factor(bpibroadclasses) + factor(bpifineclasses) + factor(HB))
gflist$gf116 = formula(presence ~ factor(sedclasses) + factor(bpifineclasses) + factor(HB))
gflist$gf117 = formula(presence ~ s(log(Patch),k=k.num) + factor(bpifineclasses) + factor(HB))
gflist$gf118 = formula(presence ~ s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
gflist$gf119 = formula(presence ~ factor(bpifineclasses) + factor(HB))


}

for (Season in c(2,4)){

# Correlation between presence and variables. Primarily used to select either Patch Perimeter or Area because
# their impacts are highly correlated. 
patchcor = data.frame(spp=NA,spname=NA,area=NA,perim=NA,usepatch=NA,bpi=NA,bathy=NA,sed=NA,compact=NA);
for (spp in u_svspp[c(1:25)]){
  spname = (species_conversions[which(species_conversions$SVSPP==spp),'Common name']); 
  ind.vec = which(u_svspp==spp)
  fishdata = alltrawl[which(alltrawl$SVSPP==spp),];
  patchcor[ind.vec,'spp'] = spp;
  patchcor[ind.vec,'spname'] = spname;
  patchcor[ind.vec,'area'] = round(cor(fishdata$presence,fishdata$PatchArea,use='pairwise.complete.obs'),4)
  patchcor[ind.vec,'perim'] = round(cor(fishdata$presence,fishdata$PatchPerim,use='pairwise.complete.obs'),4)
  patchcor[ind.vec,'bpi'] = round(cor(fishdata$presence,fishdata$bpibroadclasses,use='pairwise.complete.obs'),4)
  patchcor[ind.vec,'bathy'] = round(cor(fishdata$presence,fishdata$bathyclasses,use='pairwise.complete.obs'),4)
  patchcor[ind.vec,'sed'] = round(cor(fishdata$presence,fishdata$sedclasses,use='pairwise.complete.obs'),4)
  patchcor[ind.vec,'compact'] = round(cor(fishdata$presence,fishdata$PatchCompact,use='pairwise.complete.obs'),4)
  if (abs(patchcor[ind.vec,'perim'])>abs(patchcor[ind.vec,'area'])){
    patchcor[ind.vec,'usepatch'] = 'Perim'
  }else{
    patchcor[ind.vec,'usepatch'] = 'Area'}
}


# Runs through five partitions fitting every combination of habitat variables
for (spp in u_svspp[c(1:25)]){ 
  spname = (species_conversions[which(species_conversions$SVSPP==spp),'Common name']); print(spname)
  ind.vec = which(u_svspp==spp)
  fishdata = alltrawl[which(alltrawl$SVSPP==spp & alltrawl$SEASON==Season),];
  fishdata[which(fishdata$sedclasses==128),'sedclasses'] = NA;
  if(patchcor[ind.vec,'usepatch']=='Perim'){fishdata$Patch = fishdata$PatchPerim; habvars2$Patch = habvars2$PatchPerim}
  if(patchcor[ind.vec,'usepatch']=='Area'){fishdata$Patch = fishdata$PatchArea; habvars2$Patch = habvars2$PatchArea}

  # remove NAs
  i.na = which(is.na(fishdata[,'presence'])|is.na(fishdata$BTEMP)|is.na(fishdata$PatchArea)|is.na(fishdata$PatchPerim)|is.na(fishdata$PatchCompact)|is.na(fishdata$bathyclasses)|is.na(fishdata$bpibroadclasses)|is.na(fishdata$bpifineclasses)|is.na(fishdata$sedclasses)|is.na(fishdata$HB))
  if(length(i.na)>0){
    fishdata = fishdata[-i.na,]}

   
  habcval.df = data.frame(matrix(data=NA,nrow=5,ncol=20)); 
  names(habcval.df) = c('SVSPP','Season','partition','modelnumber','AIC','dev.exp','r.sq','auc','maxtss_thresh','maxtss_sens','maxtss_spec','prev_thresh','prev_sens','prev_spec','mult_cor','maxtss_cor','prev_cor','mult_rmse','maxtss_rmse','prev_rmse')  
  habcval.df$partition = 1:5; 
  habcval.df$SVSPP = spp; 
  habcval.df$Season = Season
  habcval.mod = vector('list',5)
  
  for (part in 1:5){
    print(part)
    fishtrain = fishdata[which(fishdata$partition!=part),]
    fishtest = fishdata[which(fishdata$partition==part),]
    
    # Degrees of freedom of GAM curve. Lower number helps prevent overfitting
    k.num = 5
    
    habcv.part = data.frame(gf=1:119, bathyclasses=NA, bpibroadclasses=NA, sedclasses=NA, bpifineclasses=NA, Patch=NA, PatchCompact=NA, HB=NA, AIC=NA, devexp=NA, R2=NA)
    
    # Fit model, record p value of variables and basic fit statistics
    gf1 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num))
    gammod = gam(gf1, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[1,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[1,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[1,'sedclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[1,'bpifineclasses'] = NA
    habcv.part[1,'HB'] = NA;
    habcv.part[1,'Patch'] = gamanova$s.table[1,4]
    habcv.part[1,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[1,'AIC'] = gammod$aic
    habcv.part[1,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[1,'R2'] = round(gamanova$r.sq,3);
    
    gf2 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num))
    gammod = gam(gf2, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[2,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[2,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[2,'sedclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[2,'bpifineclasses'] = NA
    habcv.part[2,'HB'] = NA;
    habcv.part[2,'Patch'] = gamanova$s.table[1,4]
    habcv.part[2,'PatchCompact'] = NA
    habcv.part[2,'AIC'] = gammod$aic
    habcv.part[2,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[2,'R2'] = round(gamanova$r.sq,3);
    
    gf3 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(PatchCompact,k=k.num))
    gammod = gam(gf3, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[3,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[3,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[3,'sedclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[3,'bpifineclasses'] = NA
    habcv.part[3,'HB'] = NA;
    habcv.part[3,'Patch'] = NA
    habcv.part[3,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[3,'AIC'] = gammod$aic
    habcv.part[3,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[3,'R2'] = round(gamanova$r.sq,3);
    
    gf4 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num)) 
    gammod = gam(gf4, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[4,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[4,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[4,'sedclasses'] = NA
    habcv.part[4,'bpifineclasses'] = NA
    habcv.part[4,'HB'] = NA;
    habcv.part[4,'Patch'] = gamanova$s.table[1,4]
    habcv.part[4,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[4,'AIC'] = gammod$aic
    habcv.part[4,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[4,'R2'] = round(gamanova$r.sq,3);
    
    gf5 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num))
    gammod = gam(gf5, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[5,'bathyclasses'] = NA
    habcv.part[5,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[5,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[5,'bpifineclasses'] = NA
    habcv.part[5,'HB'] = NA;
    habcv.part[5,'Patch'] = gamanova$s.table[1,4]
    habcv.part[5,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[5,'AIC'] = gammod$aic
    habcv.part[5,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[5,'R2'] = round(gamanova$r.sq,3);
    
    gf6 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses)) 
    gammod = gam(gf6, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[6,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[6,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[6,'sedclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[6,'bpifineclasses'] = NA
    habcv.part[6,'HB'] = NA;
    habcv.part[6,'Patch'] = NA
    habcv.part[6,'PatchCompact'] = NA
    habcv.part[6,'AIC'] = gammod$aic
    habcv.part[6,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[6,'R2'] = round(gamanova$r.sq,3);
    
    gf7 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(log(Patch),k=k.num)) 
    gammod = gam(gf7, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[7,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[7,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[7,'sedclasses'] = NA
    habcv.part[7,'bpifineclasses'] = NA
    habcv.part[7,'HB'] = NA;
    habcv.part[7,'Patch'] = gamanova$s.table[1,4]
    habcv.part[7,'PatchCompact'] = NA
    habcv.part[7,'AIC'] = gammod$aic
    habcv.part[7,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[7,'R2'] = round(gamanova$r.sq,3);
    
    gf8 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(PatchCompact,k=k.num))
    gammod = gam(gf8, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[8,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[8,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[8,'sedclasses'] = NA
    habcv.part[8,'bpifineclasses'] = NA
    habcv.part[8,'HB'] = NA;
    habcv.part[8,'Patch'] = NA
    habcv.part[8,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[8,'AIC'] = gammod$aic
    habcv.part[8,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[8,'R2'] = round(gamanova$r.sq,3);
    
    gf9 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + s(log(Patch),k=k.num))
    gammod = gam(gf9, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[9,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[9,'bpibroadclasses'] = NA
    habcv.part[9,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[9,'bpifineclasses'] = NA
    habcv.part[9,'HB'] = NA;
    habcv.part[9,'Patch'] = gamanova$s.table[1,4]
    habcv.part[9,'PatchCompact'] = NA
    habcv.part[9,'AIC'] = gammod$aic
    habcv.part[9,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[9,'R2'] = round(gamanova$r.sq,3);
    
    gf10 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + s(PatchCompact,k=k.num))
    gammod = gam(gf10, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[10,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[10,'bpibroadclasses'] = NA
    habcv.part[10,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[10,'bpifineclasses'] = NA
    habcv.part[10,'HB'] = NA;
    habcv.part[10,'Patch'] = NA
    habcv.part[10,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[10,'AIC'] = gammod$aic
    habcv.part[10,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[10,'R2'] = round(gamanova$r.sq,3);
    
    gf11 = formula(presence ~ factor(bathyclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num))
    gammod = gam(gf11, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[11,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[11,'bpibroadclasses'] = NA
    habcv.part[11,'sedclasses'] = NA
    habcv.part[11,'bpifineclasses'] = NA
    habcv.part[11,'HB'] = NA;
    habcv.part[11,'Patch'] = gamanova$s.table[1,4]
    habcv.part[11,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[11,'AIC'] = gammod$aic
    habcv.part[11,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[11,'R2'] = round(gamanova$r.sq,3);
    
    gf12 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num))
    gammod = gam(gf12, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[12,'bathyclasses'] = NA
    habcv.part[12,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[12,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[12,'bpifineclasses'] = NA
    habcv.part[12,'HB'] = NA;
    habcv.part[12,'Patch'] = gamanova$s.table[1,4]
    habcv.part[12,'PatchCompact'] = NA
    habcv.part[12,'AIC'] = gammod$aic
    habcv.part[12,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[12,'R2'] = round(gamanova$r.sq,3);
    
    gf13 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(PatchCompact,k=k.num))
    gammod = gam(gf13, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[13,'bathyclasses'] = NA
    habcv.part[13,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[13,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[13,'bpifineclasses'] = NA
    habcv.part[13,'HB'] = NA;
    habcv.part[13,'Patch'] = NA
    habcv.part[13,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[13,'AIC'] = gammod$aic
    habcv.part[13,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[13,'R2'] = round(gamanova$r.sq,3);
    
    gf14 = formula(presence ~ factor(bpibroadclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num)) 
    gammod = gam(gf14, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[14,'bathyclasses'] = NA
    habcv.part[14,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[14,'sedclasses'] = NA
    habcv.part[14,'bpifineclasses'] = NA
    habcv.part[14,'HB'] = NA;
    habcv.part[14,'Patch'] = gamanova$s.table[1,4]
    habcv.part[14,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[14,'AIC'] = gammod$aic
    habcv.part[14,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[14,'R2'] = round(gamanova$r.sq,3);
    
    gf15 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses))
    gammod = gam(gf15, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[15,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[15,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[15,'sedclasses'] = NA
    habcv.part[15,'bpifineclasses'] = NA
    habcv.part[15,'HB'] = NA;
    habcv.part[15,'Patch'] = NA
    habcv.part[15,'PatchCompact'] = NA
    habcv.part[15,'AIC'] = gammod$aic
    habcv.part[15,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[15,'R2'] = round(gamanova$r.sq,3);
    
    gf16 = formula(presence ~ factor(bathyclasses) + factor(sedclasses))
    gammod = gam(gf16, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[16,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[16,'bpibroadclasses'] = NA
    habcv.part[16,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[16,'bpifineclasses'] = NA
    habcv.part[16,'HB'] = NA;
    habcv.part[16,'Patch'] =NA
    habcv.part[16,'PatchCompact'] = NA
    habcv.part[16,'AIC'] = gammod$aic
    habcv.part[16,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[16,'R2'] = round(gamanova$r.sq,3);
    
    gf17 = formula(presence ~ factor(bathyclasses) + s(log(Patch),k=k.num))
    gammod = gam(gf17, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[17,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[17,'bpibroadclasses'] = NA
    habcv.part[17,'sedclasses'] = NA
    habcv.part[17,'bpifineclasses'] = NA
    habcv.part[17,'HB'] = NA;
    habcv.part[17,'Patch'] = gamanova$s.table[1,4]
    habcv.part[17,'PatchCompact'] = NA
    habcv.part[17,'AIC'] = gammod$aic
    habcv.part[17,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[17,'R2'] = round(gamanova$r.sq,3);
    
    gf18 = formula(presence ~ factor(bathyclasses) + s(PatchCompact,k=k.num))
    gammod = gam(gf18, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[18,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[18,'bpibroadclasses'] = NA
    habcv.part[18,'sedclasses'] = NA
    habcv.part[18,'bpifineclasses'] = NA
    habcv.part[18,'HB'] = NA;
    habcv.part[18,'Patch'] =NA
    habcv.part[18,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[18,'AIC'] = gammod$aic
    habcv.part[18,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[18,'R2'] = round(gamanova$r.sq,3);
    
    gf19 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses)) 
    gammod = gam(gf19, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[19,'bathyclasses'] = NA
    habcv.part[19,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[19,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[19,'bpifineclasses'] = NA
    habcv.part[19,'HB'] = NA;
    habcv.part[19,'Patch'] = NA
    habcv.part[19,'PatchCompact'] = NA
    habcv.part[19,'AIC'] = gammod$aic
    habcv.part[19,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[19,'R2'] = round(gamanova$r.sq,3);
    
    gf20 = formula(presence ~ factor(bpibroadclasses) + s(log(Patch),k=k.num)) 
    gammod = gam(gf20, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[20,'bathyclasses'] = NA
    habcv.part[20,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[20,'sedclasses'] = NA
    habcv.part[20,'bpifineclasses'] = NA
    habcv.part[20,'HB'] = NA;
    habcv.part[20,'Patch'] = gamanova$s.table[1,4]
    habcv.part[20,'PatchCompact'] = NA
    habcv.part[20,'AIC'] = gammod$aic
    habcv.part[20,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[20,'R2'] = round(gamanova$r.sq,3);
    
    gf21 = formula(presence ~ factor(bpibroadclasses) + s(PatchCompact,k=k.num))
    gammod = gam(gf21, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[21,'bathyclasses'] = NA
    habcv.part[21,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[21,'sedclasses'] = NA
    habcv.part[21,'bpifineclasses'] = NA
    habcv.part[21,'HB'] = NA;
    habcv.part[21,'Patch'] = NA
    habcv.part[21,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[21,'AIC'] = gammod$aic
    habcv.part[21,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[21,'R2'] = round(gamanova$r.sq,3);
    
    gf22 = formula(presence ~ factor(sedclasses) + s(log(Patch),k=k.num))
    gammod = gam(gf22, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[22,'bathyclasses'] = NA
    habcv.part[22,'bpibroadclasses'] = NA
    habcv.part[22,'sedclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[22,'bpifineclasses'] = NA
    habcv.part[22,'HB'] = NA;
    habcv.part[22,'Patch'] = gamanova$s.table[1,4]
    habcv.part[22,'PatchCompact'] = NA
    habcv.part[22,'AIC'] = gammod$aic
    habcv.part[22,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[22,'R2'] = round(gamanova$r.sq,3);
    
    gf23 = formula(presence ~ factor(sedclasses) + s(PatchCompact,k=k.num))
    gammod = gam(gf23, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[23,'bathyclasses'] = NA
    habcv.part[23,'bpibroadclasses'] = NA
    habcv.part[23,'sedclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[23,'bpifineclasses'] = NA
    habcv.part[23,'HB'] = NA;
    habcv.part[23,'Patch'] = NA
    habcv.part[23,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[23,'AIC'] = gammod$aic
    habcv.part[23,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[23,'R2'] = round(gamanova$r.sq,3);
    
    gf24 = formula(presence ~ s(log(Patch),k=k.num) + s(PatchCompact,k=k.num))
    gammod = gam(gf24, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[24,'bathyclasses'] = NA
    habcv.part[24,'bpibroadclasses'] = NA
    habcv.part[24,'sedclasses'] = NA
    habcv.part[24,'bpifineclasses'] = NA
    habcv.part[24,'HB'] = NA;
    habcv.part[24,'Patch'] = gamanova$s.table[1,4]
    habcv.part[24,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[24,'AIC'] = gammod$aic
    habcv.part[24,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[24,'R2'] = round(gamanova$r.sq,3);
    
    gf25 = formula(presence ~ factor(bathyclasses))
    gammod = gam(gf25, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[25,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[25,'bpibroadclasses'] = NA
    habcv.part[25,'sedclasses'] = NA
    habcv.part[25,'bpifineclasses'] = NA
    habcv.part[25,'HB'] = NA;
    habcv.part[25,'Patch'] = NA
    habcv.part[25,'PatchCompact'] = NA
    habcv.part[25,'AIC'] = gammod$aic
    habcv.part[25,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[25,'R2'] = round(gamanova$r.sq,3);
    
    gf26 = formula(presence ~ factor(bpibroadclasses))
    gammod = gam(gf26, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[26,'bathyclasses'] = NA
    habcv.part[26,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[26,'sedclasses'] = NA
    habcv.part[26,'bpifineclasses'] = NA
    habcv.part[26,'HB'] = NA;
    habcv.part[26,'Patch'] = NA
    habcv.part[26,'PatchCompact'] = NA
    habcv.part[26,'AIC'] = gammod$aic
    habcv.part[26,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[26,'R2'] = round(gamanova$r.sq,3);
    
    gf27 = formula(presence ~ factor(sedclasses))
    gammod = gam(gf27, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[27,'bathyclasses'] = NA
    habcv.part[27,'bpibroadclasses'] = NA
    habcv.part[27,'sedclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[27,'bpifineclasses'] = NA
    habcv.part[27,'HB'] = NA;
    habcv.part[27,'Patch'] = NA
    habcv.part[27,'PatchCompact'] = NA
    habcv.part[27,'AIC'] = gammod$aic
    habcv.part[27,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[27,'R2'] = round(gamanova$r.sq,3);
    
    gf28 = formula(presence ~ s(log(Patch),k=k.num))
    gammod = gam(gf28, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[28,'bathyclasses'] = NA
    habcv.part[28,'bpibroadclasses'] = NA
    habcv.part[28,'sedclasses'] = NA
    habcv.part[28,'bpifineclasses'] = NA
    habcv.part[28,'HB'] = NA;
    habcv.part[28,'Patch'] = gamanova$s.table[1,4]
    habcv.part[28,'PatchCompact'] = NA
    habcv.part[28,'AIC'] = gammod$aic
    habcv.part[28,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[28,'R2'] = round(gamanova$r.sq,3);
    
    gf29 = formula(presence ~ s(PatchCompact,k=k.num))
    gammod = gam(gf29, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[29,'bathyclasses'] = NA
    habcv.part[29,'bpibroadclasses'] = NA
    habcv.part[29,'sedclasses'] = NA
    habcv.part[29,'bpifineclasses'] = NA
    habcv.part[29,'HB'] = NA;
    habcv.part[29,'Patch'] = NA
    habcv.part[29,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[29,'AIC'] = gammod$aic
    habcv.part[29,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[29,'R2'] = round(gamanova$r.sq,3);
    
    ### Fine included
    
    gf30 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
    gammod = gam(gf30, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[30,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[30,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[30,'sedclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[30,'bpifineclasses'] = gamanova$pTerms.table[4,3]
    habcv.part[30,'HB'] = NA;
    habcv.part[30,'Patch'] = gamanova$s.table[1,4]
    habcv.part[30,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[30,'AIC'] = gammod$aic
    habcv.part[30,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[30,'R2'] = round(gamanova$r.sq,3);
    
    gf31 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses))
    gammod = gam(gf31, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[31,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[31,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[31,'sedclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[31,'bpifineclasses'] = gamanova$pTerms.table[4,3]
    habcv.part[31,'HB'] = NA;
    habcv.part[31,'Patch'] = gamanova$s.table[1,4]
    habcv.part[31,'PatchCompact'] = NA
    habcv.part[31,'AIC'] = gammod$aic
    habcv.part[31,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[31,'R2'] = round(gamanova$r.sq,3);
    
    gf32 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
    gammod = gam(gf32, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[32,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[32,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[32,'sedclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[32,'bpifineclasses'] = gamanova$pTerms.table[4,3]
    habcv.part[32,'HB'] = NA;
    habcv.part[32,'Patch'] = NA
    habcv.part[32,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[32,'AIC'] = gammod$aic
    habcv.part[32,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[32,'R2'] = round(gamanova$r.sq,3);
    
    gf33 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses)) 
    gammod = gam(gf33, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[33,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[33,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[33,'sedclasses'] = NA
    habcv.part[33,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[33,'HB'] = NA;
    habcv.part[33,'Patch'] = gamanova$s.table[1,4]
    habcv.part[33,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[33,'AIC'] = gammod$aic
    habcv.part[33,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[33,'R2'] = round(gamanova$r.sq,3);
    
    gf34 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
    gammod = gam(gf34, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[34,'bathyclasses'] = NA
    habcv.part[34,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[34,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[34,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[34,'HB'] = NA;
    habcv.part[34,'Patch'] = gamanova$s.table[1,4]
    habcv.part[34,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[34,'AIC'] = gammod$aic
    habcv.part[34,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[34,'R2'] = round(gamanova$r.sq,3);
    
    gf35 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + factor(bpifineclasses)) 
    gammod = gam(gf35, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[35,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[35,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[35,'sedclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[35,'bpifineclasses'] = gamanova$pTerms.table[4,3]
    habcv.part[35,'HB'] = NA;
    habcv.part[35,'Patch'] = NA
    habcv.part[35,'PatchCompact'] = NA
    habcv.part[35,'AIC'] = gammod$aic
    habcv.part[35,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[35,'R2'] = round(gamanova$r.sq,3);
    
    gf36 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses)) 
    gammod = gam(gf36, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[36,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[36,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[36,'sedclasses'] = NA
    habcv.part[36,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[36,'HB'] = NA;
    habcv.part[36,'Patch'] = gamanova$s.table[1,4]
    habcv.part[36,'PatchCompact'] = NA
    habcv.part[36,'AIC'] = gammod$aic
    habcv.part[36,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[36,'R2'] = round(gamanova$r.sq,3);
    
    gf37 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
    gammod = gam(gf37, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[37,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[37,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[37,'sedclasses'] = NA
    habcv.part[37,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[37,'HB'] = NA;
    habcv.part[37,'Patch'] = NA
    habcv.part[37,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[37,'AIC'] = gammod$aic
    habcv.part[37,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[37,'R2'] = round(gamanova$r.sq,3);
    
    gf38 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses))
    gammod = gam(gf38, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[38,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[38,'bpibroadclasses'] = NA
    habcv.part[38,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[38,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[38,'HB'] = NA;
    habcv.part[38,'Patch'] = gamanova$s.table[1,4]
    habcv.part[38,'PatchCompact'] = NA
    habcv.part[38,'AIC'] = gammod$aic
    habcv.part[38,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[38,'R2'] = round(gamanova$r.sq,3);
    
    gf39 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
    gammod = gam(gf39, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[39,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[39,'bpibroadclasses'] = NA
    habcv.part[39,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[39,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[39,'HB'] = NA;
    habcv.part[39,'Patch'] = NA
    habcv.part[39,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[39,'AIC'] = gammod$aic
    habcv.part[39,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[39,'R2'] = round(gamanova$r.sq,3);
    
    gf40 = formula(presence ~ factor(bathyclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
    gammod = gam(gf40, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[40,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[40,'bpibroadclasses'] = NA
    habcv.part[40,'sedclasses'] = NA
    habcv.part[40,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[40,'HB'] = NA;
    habcv.part[40,'Patch'] = gamanova$s.table[1,4]
    habcv.part[40,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[40,'AIC'] = gammod$aic
    habcv.part[40,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[40,'R2'] = round(gamanova$r.sq,3);
    
    gf41 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses))
    gammod = gam(gf41, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[41,'bathyclasses'] = NA
    habcv.part[41,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[41,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[41,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[41,'HB'] = NA;
    habcv.part[41,'Patch'] = gamanova$s.table[1,4]
    habcv.part[41,'PatchCompact'] = NA
    habcv.part[41,'AIC'] = gammod$aic
    habcv.part[41,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[41,'R2'] = round(gamanova$r.sq,3);
    
    gf42 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
    gammod = gam(gf42, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[42,'bathyclasses'] = NA
    habcv.part[42,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[42,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[42,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[42,'HB'] = NA;
    habcv.part[42,'Patch'] = NA
    habcv.part[42,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[42,'AIC'] = gammod$aic
    habcv.part[42,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[42,'R2'] = round(gamanova$r.sq,3);
    
    gf43 = formula(presence ~ factor(bpibroadclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses)) 
    gammod = gam(gf43, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[43,'bathyclasses'] = NA
    habcv.part[43,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[43,'sedclasses'] = NA
    habcv.part[43,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[43,'HB'] = NA;
    habcv.part[43,'Patch'] = gamanova$s.table[1,4]
    habcv.part[43,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[43,'AIC'] = gammod$aic
    habcv.part[43,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[43,'R2'] = round(gamanova$r.sq,3);
    
    gf44 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(bpifineclasses))
    gammod = gam(gf44, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[44,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[44,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[44,'sedclasses'] = NA
    habcv.part[44,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[44,'HB'] = NA;
    habcv.part[44,'Patch'] = NA
    habcv.part[44,'PatchCompact'] = NA
    habcv.part[44,'AIC'] = gammod$aic
    habcv.part[44,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[44,'R2'] = round(gamanova$r.sq,3);
    
    gf45 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + factor(bpifineclasses))
    gammod = gam(gf45, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[45,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[45,'bpibroadclasses'] = NA
    habcv.part[45,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[45,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[45,'HB'] = NA;
    habcv.part[45,'Patch'] =NA
    habcv.part[45,'PatchCompact'] = NA
    habcv.part[45,'AIC'] = gammod$aic
    habcv.part[45,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[45,'R2'] = round(gamanova$r.sq,3);
    
    gf46 = formula(presence ~ factor(bathyclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses))
    gammod = gam(gf46, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[46,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[46,'bpibroadclasses'] = NA
    habcv.part[46,'sedclasses'] = NA
    habcv.part[46,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[46,'HB'] = NA;
    habcv.part[46,'Patch'] = gamanova$s.table[1,4]
    habcv.part[46,'PatchCompact'] = NA
    habcv.part[46,'AIC'] = gammod$aic
    habcv.part[46,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[46,'R2'] = round(gamanova$r.sq,3);
    
    gf47 = formula(presence ~ factor(bathyclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
    gammod = gam(gf47, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[47,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[47,'bpibroadclasses'] = NA
    habcv.part[47,'sedclasses'] = NA
    habcv.part[47,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[47,'HB'] = NA;
    habcv.part[47,'Patch'] =NA
    habcv.part[47,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[47,'AIC'] = gammod$aic
    habcv.part[47,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[47,'R2'] = round(gamanova$r.sq,3);
    
    gf48 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + factor(bpifineclasses)) 
    gammod = gam(gf48, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[48,'bathyclasses'] = NA
    habcv.part[48,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[48,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[48,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[48,'HB'] = NA;
    habcv.part[48,'Patch'] = NA
    habcv.part[48,'PatchCompact'] = NA
    habcv.part[48,'AIC'] = gammod$aic
    habcv.part[48,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[48,'R2'] = round(gamanova$r.sq,3);
    
    gf49 = formula(presence ~ factor(bpibroadclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses)) 
    gammod = gam(gf49, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[49,'bathyclasses'] = NA
    habcv.part[49,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[49,'sedclasses'] = NA
    habcv.part[49,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[49,'HB'] = NA;
    habcv.part[49,'Patch'] = gamanova$s.table[1,4]
    habcv.part[49,'PatchCompact'] = NA
    habcv.part[49,'AIC'] = gammod$aic
    habcv.part[49,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[49,'R2'] = round(gamanova$r.sq,3);
    
    gf50 = formula(presence ~ factor(bpibroadclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
    gammod = gam(gf50, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[50,'bathyclasses'] = NA
    habcv.part[50,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[50,'sedclasses'] = NA
    habcv.part[50,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[50,'HB'] = NA;
    habcv.part[50,'Patch'] = NA
    habcv.part[50,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[50,'AIC'] = gammod$aic
    habcv.part[50,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[50,'R2'] = round(gamanova$r.sq,3);
    
    gf51 = formula(presence ~ factor(sedclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses))
    gammod = gam(gf51, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[51,'bathyclasses'] = NA
    habcv.part[51,'bpibroadclasses'] = NA
    habcv.part[51,'sedclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[51,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[51,'HB'] = NA;
    habcv.part[51,'Patch'] = gamanova$s.table[1,4]
    habcv.part[51,'PatchCompact'] = NA
    habcv.part[51,'AIC'] = gammod$aic
    habcv.part[51,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[51,'R2'] = round(gamanova$r.sq,3);
    
    gf52 = formula(presence ~ factor(sedclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
    gammod = gam(gf52, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[52,'bathyclasses'] = NA
    habcv.part[52,'bpibroadclasses'] = NA
    habcv.part[52,'sedclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[52,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[52,'HB'] = NA;
    habcv.part[52,'Patch'] = NA
    habcv.part[52,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[52,'AIC'] = gammod$aic
    habcv.part[52,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[52,'R2'] = round(gamanova$r.sq,3);
    
    gf53 = formula(presence ~ s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses))
    gammod = gam(gf53, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[53,'bathyclasses'] = NA
    habcv.part[53,'bpibroadclasses'] = NA
    habcv.part[53,'sedclasses'] = NA
    habcv.part[53,'bpifineclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[53,'HB'] = NA;
    habcv.part[53,'Patch'] = gamanova$s.table[1,4]
    habcv.part[53,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[53,'AIC'] = gammod$aic
    habcv.part[53,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[53,'R2'] = round(gamanova$r.sq,3);
    
    gf54 = formula(presence ~ factor(bathyclasses) + factor(bpifineclasses))
    gammod = gam(gf54, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[54,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[54,'bpibroadclasses'] = NA
    habcv.part[54,'sedclasses'] = NA
    habcv.part[54,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[54,'HB'] = NA;
    habcv.part[54,'Patch'] = NA
    habcv.part[54,'PatchCompact'] = NA
    habcv.part[54,'AIC'] = gammod$aic
    habcv.part[54,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[54,'R2'] = round(gamanova$r.sq,3);
    
    gf55 = formula(presence ~ factor(bpibroadclasses) + factor(bpifineclasses))
    gammod = gam(gf55, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[55,'bathyclasses'] = NA
    habcv.part[55,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[55,'sedclasses'] = NA
    habcv.part[55,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[55,'HB'] = NA;
    habcv.part[55,'Patch'] = NA
    habcv.part[55,'PatchCompact'] = NA
    habcv.part[55,'AIC'] = gammod$aic
    habcv.part[55,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[55,'R2'] = round(gamanova$r.sq,3);
    
    gf56 = formula(presence ~ factor(sedclasses) + factor(bpifineclasses))
    gammod = gam(gf56, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[56,'bathyclasses'] = NA
    habcv.part[56,'bpibroadclasses'] = NA
    habcv.part[56,'sedclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[56,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[56,'HB'] = NA;
    habcv.part[56,'Patch'] = NA
    habcv.part[56,'PatchCompact'] = NA
    habcv.part[56,'AIC'] = gammod$aic
    habcv.part[56,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[56,'R2'] = round(gamanova$r.sq,3);
    
    gf57 = formula(presence ~ s(log(Patch),k=k.num) + factor(bpifineclasses))
    gammod = gam(gf57, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[57,'bathyclasses'] = NA
    habcv.part[57,'bpibroadclasses'] = NA
    habcv.part[57,'sedclasses'] = NA
    habcv.part[57,'bpifineclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[57,'HB'] = NA;
    habcv.part[57,'Patch'] = gamanova$s.table[1,4]
    habcv.part[57,'PatchCompact'] = NA
    habcv.part[57,'AIC'] = gammod$aic
    habcv.part[57,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[57,'R2'] = round(gamanova$r.sq,3);
    
    gf58 = formula(presence ~ s(PatchCompact,k=k.num) + factor(bpifineclasses))
    gammod = gam(gf58, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[58,'bathyclasses'] = NA
    habcv.part[58,'bpibroadclasses'] = NA
    habcv.part[58,'sedclasses'] = NA
    habcv.part[58,'bpifineclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[58,'HB'] = NA;
    habcv.part[58,'Patch'] = NA
    habcv.part[58,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[58,'AIC'] = gammod$aic
    habcv.part[58,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[58,'R2'] = round(gamanova$r.sq,3);
    
    gf59 = formula(presence ~ factor(bpifineclasses))
    gammod = gam(gf59, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[59,'bathyclasses'] = NA
    habcv.part[59,'bpibroadclasses'] = NA
    habcv.part[59,'sedclasses'] = NA
    habcv.part[59,'bpifineclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[59,'HB'] = NA;
    habcv.part[59,'Patch'] = NA
    habcv.part[59,'PatchCompact'] = NA
    habcv.part[59,'AIC'] = gammod$aic
    habcv.part[59,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[59,'R2'] = round(gamanova$r.sq,3);

    # HB
    gf60 = formula(presence ~ factor(HB))
    gammod = gam(gf60, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[60,'bathyclasses'] = NA
    habcv.part[60,'bpibroadclasses'] = NA
    habcv.part[60,'sedclasses'] = NA
    habcv.part[60,'bpifineclasses'] = NA
    habcv.part[60,'HB'] = gamanova$pTerms.table[1,3];
    habcv.part[60,'Patch'] = NA
    habcv.part[60,'PatchCompact'] = NA
    habcv.part[60,'AIC'] = gammod$aic
    habcv.part[60,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[60,'R2'] = round(gamanova$r.sq,3);
    
    gf61 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(HB))
    gammod = gam(gf61, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[61,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[61,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[61,'sedclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[61,'bpifineclasses'] = NA
    habcv.part[61,'HB'] = gamanova$pTerms.table[4,3];
    habcv.part[61,'Patch'] = gamanova$s.table[1,4]
    habcv.part[61,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[61,'AIC'] = gammod$aic
    habcv.part[61,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[61,'R2'] = round(gamanova$r.sq,3);
    
    gf62 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + factor(HB))
    gammod = gam(gf62, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[62,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[62,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[62,'sedclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[62,'bpifineclasses'] = NA
    habcv.part[62,'HB'] = gamanova$pTerms.table[4,3];
    habcv.part[62,'Patch'] = gamanova$s.table[1,4]
    habcv.part[62,'PatchCompact'] = NA
    habcv.part[62,'AIC'] = gammod$aic
    habcv.part[62,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[62,'R2'] = round(gamanova$r.sq,3);
    
    gf63 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(PatchCompact,k=k.num) + factor(HB))
    gammod = gam(gf63, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[63,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[63,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[63,'sedclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[63,'bpifineclasses'] = NA
    habcv.part[63,'HB'] = gamanova$pTerms.table[4,3];
    habcv.part[63,'Patch'] = NA
    habcv.part[63,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[63,'AIC'] = gammod$aic
    habcv.part[63,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[63,'R2'] = round(gamanova$r.sq,3);
    
    gf64 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(HB)) 
    gammod = gam(gf64, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[64,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[64,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[64,'sedclasses'] = NA
    habcv.part[64,'bpifineclasses'] = NA
    habcv.part[64,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[64,'Patch'] = gamanova$s.table[1,4]
    habcv.part[64,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[64,'AIC'] = gammod$aic
    habcv.part[64,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[64,'R2'] = round(gamanova$r.sq,3);
    
    gf65 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(HB))
    gammod = gam(gf65, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[65,'bathyclasses'] = NA
    habcv.part[65,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[65,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[65,'bpifineclasses'] = NA
    habcv.part[65,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[65,'Patch'] = gamanova$s.table[1,4]
    habcv.part[65,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[65,'AIC'] = gammod$aic
    habcv.part[65,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[65,'R2'] = round(gamanova$r.sq,3);
    
    gf66 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + factor(HB)) 
    gammod = gam(gf66, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[66,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[66,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[66,'sedclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[66,'bpifineclasses'] = NA
    habcv.part[66,'HB'] = gamanova$pTerms.table[4,3];
    habcv.part[66,'Patch'] = NA
    habcv.part[66,'PatchCompact'] = NA
    habcv.part[66,'AIC'] = gammod$aic
    habcv.part[66,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[66,'R2'] = round(gamanova$r.sq,3);
    
    gf67 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(log(Patch),k=k.num) + factor(HB)) 
    gammod = gam(gf67, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[67,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[67,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[67,'sedclasses'] = NA
    habcv.part[67,'bpifineclasses'] = NA
    habcv.part[67,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[67,'Patch'] = gamanova$s.table[1,4]
    habcv.part[67,'PatchCompact'] = NA
    habcv.part[67,'AIC'] = gammod$aic
    habcv.part[67,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[67,'R2'] = round(gamanova$r.sq,3);
    
    gf68 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(PatchCompact,k=k.num) + factor(HB))
    gammod = gam(gf68, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[68,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[68,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[68,'sedclasses'] = NA
    habcv.part[68,'bpifineclasses'] = NA
    habcv.part[68,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[68,'Patch'] = NA
    habcv.part[68,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[68,'AIC'] = gammod$aic
    habcv.part[68,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[68,'R2'] = round(gamanova$r.sq,3);
    
    gf69 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + factor(HB))
    gammod = gam(gf69, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[69,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[69,'bpibroadclasses'] = NA
    habcv.part[69,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[69,'bpifineclasses'] = NA
    habcv.part[69,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[69,'Patch'] = gamanova$s.table[1,4]
    habcv.part[69,'PatchCompact'] = NA
    habcv.part[69,'AIC'] = gammod$aic
    habcv.part[69,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[69,'R2'] = round(gamanova$r.sq,3);
    
    gf70 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + s(PatchCompact,k=k.num) + factor(HB))
    gammod = gam(gf70, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[70,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[70,'bpibroadclasses'] = NA
    habcv.part[70,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[70,'bpifineclasses'] = NA
    habcv.part[70,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[70,'Patch'] = NA
    habcv.part[70,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[70,'AIC'] = gammod$aic
    habcv.part[70,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[70,'R2'] = round(gamanova$r.sq,3);
    
    gf71 = formula(presence ~ factor(bathyclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(HB))
    gammod = gam(gf71, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[71,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[71,'bpibroadclasses'] = NA
    habcv.part[71,'sedclasses'] = NA
    habcv.part[71,'bpifineclasses'] = NA
    habcv.part[71,'HB'] = gamanova$pTerms.table[2,3];
    habcv.part[71,'Patch'] = gamanova$s.table[1,4]
    habcv.part[71,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[71,'AIC'] = gammod$aic
    habcv.part[71,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[71,'R2'] = round(gamanova$r.sq,3);
    
    gf72 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + factor(HB))
    gammod = gam(gf72, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[72,'bathyclasses'] = NA
    habcv.part[72,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[72,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[72,'bpifineclasses'] = NA
    habcv.part[72,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[72,'Patch'] = gamanova$s.table[1,4]
    habcv.part[72,'PatchCompact'] = NA
    habcv.part[72,'AIC'] = gammod$aic
    habcv.part[72,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[72,'R2'] = round(gamanova$r.sq,3);
    
    gf73 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(PatchCompact,k=k.num) + factor(HB))
    gammod = gam(gf73, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[73,'bathyclasses'] = NA
    habcv.part[73,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[73,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[73,'bpifineclasses'] = NA
    habcv.part[73,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[73,'Patch'] = NA
    habcv.part[73,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[73,'AIC'] = gammod$aic
    habcv.part[73,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[73,'R2'] = round(gamanova$r.sq,3);
    
    gf74 = formula(presence ~ factor(bpibroadclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(HB)) 
    gammod = gam(gf74, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[74,'bathyclasses'] = NA
    habcv.part[74,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[74,'sedclasses'] = NA
    habcv.part[74,'bpifineclasses'] = NA
    habcv.part[74,'HB'] = gamanova$pTerms.table[2,3];
    habcv.part[74,'Patch'] = gamanova$s.table[1,4]
    habcv.part[74,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[74,'AIC'] = gammod$aic
    habcv.part[74,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[74,'R2'] = round(gamanova$r.sq,3);
    
    gf75 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(HB))
    gammod = gam(gf75, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[75,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[75,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[75,'sedclasses'] = NA
    habcv.part[75,'bpifineclasses'] = NA
    habcv.part[75,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[75,'Patch'] = NA
    habcv.part[75,'PatchCompact'] = NA
    habcv.part[75,'AIC'] = gammod$aic
    habcv.part[75,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[75,'R2'] = round(gamanova$r.sq,3);
    
    gf76 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + factor(HB))
    gammod = gam(gf76, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[76,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[76,'bpibroadclasses'] = NA
    habcv.part[76,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[76,'bpifineclasses'] = NA
    habcv.part[76,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[76,'Patch'] =NA
    habcv.part[76,'PatchCompact'] = NA
    habcv.part[76,'AIC'] = gammod$aic
    habcv.part[76,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[76,'R2'] = round(gamanova$r.sq,3);
    
    gf77 = formula(presence ~ factor(bathyclasses) + s(log(Patch),k=k.num) + factor(HB))
    gammod = gam(gf77, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[77,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[77,'bpibroadclasses'] = NA
    habcv.part[77,'sedclasses'] = NA
    habcv.part[77,'bpifineclasses'] = NA
    habcv.part[77,'HB'] = gamanova$pTerms.table[2,3];
    habcv.part[77,'Patch'] = gamanova$s.table[1,4]
    habcv.part[77,'PatchCompact'] = NA
    habcv.part[77,'AIC'] = gammod$aic
    habcv.part[77,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[77,'R2'] = round(gamanova$r.sq,3);
    
    gf78 = formula(presence ~ factor(bathyclasses) + s(PatchCompact,k=k.num) + factor(HB))
    gammod = gam(gf78, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[78,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[78,'bpibroadclasses'] = NA
    habcv.part[78,'sedclasses'] = NA
    habcv.part[78,'bpifineclasses'] = NA
    habcv.part[78,'HB'] = gamanova$pTerms.table[2,3];
    habcv.part[78,'Patch'] =NA
    habcv.part[78,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[78,'AIC'] = gammod$aic
    habcv.part[78,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[78,'R2'] = round(gamanova$r.sq,3);
    
    gf79 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + factor(HB)) 
    gammod = gam(gf79, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[79,'bathyclasses'] = NA
    habcv.part[79,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[79,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[79,'bpifineclasses'] = NA
    habcv.part[79,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[79,'Patch'] = NA
    habcv.part[79,'PatchCompact'] = NA
    habcv.part[79,'AIC'] = gammod$aic
    habcv.part[79,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[79,'R2'] = round(gamanova$r.sq,3);
    
    gf80 = formula(presence ~ factor(bpibroadclasses) + s(log(Patch),k=k.num) + factor(HB)) 
    gammod = gam(gf80, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[80,'bathyclasses'] = NA
    habcv.part[80,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[80,'sedclasses'] = NA
    habcv.part[80,'bpifineclasses'] = NA
    habcv.part[80,'HB'] = gamanova$pTerms.table[2,3];
    habcv.part[80,'Patch'] = gamanova$s.table[1,4]
    habcv.part[80,'PatchCompact'] = NA
    habcv.part[80,'AIC'] = gammod$aic
    habcv.part[80,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[80,'R2'] = round(gamanova$r.sq,3);
    
    gf81 = formula(presence ~ factor(bpibroadclasses) + s(PatchCompact,k=k.num) + factor(HB))
    gammod = gam(gf81, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[81,'bathyclasses'] = NA
    habcv.part[81,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[81,'sedclasses'] = NA
    habcv.part[81,'bpifineclasses'] = NA
    habcv.part[81,'HB'] = gamanova$pTerms.table[2,3];
    habcv.part[81,'Patch'] = NA
    habcv.part[81,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[81,'AIC'] = gammod$aic
    habcv.part[81,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[81,'R2'] = round(gamanova$r.sq,3);
    
    gf82 = formula(presence ~ factor(sedclasses) + s(log(Patch),k=k.num) + factor(HB))
    gammod = gam(gf82, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[82,'bathyclasses'] = NA
    habcv.part[82,'bpibroadclasses'] = NA
    habcv.part[82,'sedclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[82,'bpifineclasses'] = NA
    habcv.part[82,'HB'] = gamanova$pTerms.table[2,3];
    habcv.part[82,'Patch'] = gamanova$s.table[1,4]
    habcv.part[82,'PatchCompact'] = NA
    habcv.part[82,'AIC'] = gammod$aic
    habcv.part[82,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[82,'R2'] = round(gamanova$r.sq,3);
  
    gf83 = formula(presence ~ factor(sedclasses) + s(PatchCompact,k=k.num) + factor(HB))
    gammod = gam(gf83, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[83,'bathyclasses'] = NA
    habcv.part[83,'bpibroadclasses'] = NA
    habcv.part[83,'sedclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[83,'bpifineclasses'] = NA
    habcv.part[83,'HB'] = gamanova$pTerms.table[2,3];
    habcv.part[83,'Patch'] = NA
    habcv.part[83,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[83,'AIC'] = gammod$aic
    habcv.part[83,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[83,'R2'] = round(gamanova$r.sq,3);
  
    gf84 = formula(presence ~ s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(HB))
    gammod = gam(gf84, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[84,'bathyclasses'] = NA
    habcv.part[84,'bpibroadclasses'] = NA
    habcv.part[84,'sedclasses'] = NA
    habcv.part[84,'bpifineclasses'] = NA
    habcv.part[84,'HB'] = gamanova$pTerms.table[1,3];
    habcv.part[84,'Patch'] = gamanova$s.table[1,4]
    habcv.part[84,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[84,'AIC'] = gammod$aic
    habcv.part[84,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[84,'R2'] = round(gamanova$r.sq,3);
    
    gf85 = formula(presence ~ factor(bathyclasses) + factor(HB))
    gammod = gam(gf85, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[85,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[85,'bpibroadclasses'] = NA
    habcv.part[85,'sedclasses'] = NA
    habcv.part[85,'bpifineclasses'] = NA
    habcv.part[85,'HB'] = gamanova$pTerms.table[2,3];
    habcv.part[85,'Patch'] = NA
    habcv.part[85,'PatchCompact'] = NA
    habcv.part[85,'AIC'] = gammod$aic
    habcv.part[85,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[85,'R2'] = round(gamanova$r.sq,3);
    
    gf86 = formula(presence ~ factor(bpibroadclasses) + factor(HB))
    gammod = gam(gf86, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[86,'bathyclasses'] = NA
    habcv.part[86,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[86,'sedclasses'] = NA
    habcv.part[86,'bpifineclasses'] = NA
    habcv.part[86,'HB'] = gamanova$pTerms.table[2,3];
    habcv.part[86,'Patch'] = NA
    habcv.part[86,'PatchCompact'] = NA
    habcv.part[86,'AIC'] = gammod$aic
    habcv.part[86,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[86,'R2'] = round(gamanova$r.sq,3);
  
    gf87 = formula(presence ~ factor(sedclasses) + factor(HB))
    gammod = gam(gf87, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[87,'bathyclasses'] = NA
    habcv.part[87,'bpibroadclasses'] = NA
    habcv.part[87,'sedclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[87,'bpifineclasses'] = NA
    habcv.part[87,'HB'] = gamanova$pTerms.table[2,3];
    habcv.part[87,'Patch'] = NA
    habcv.part[87,'PatchCompact'] = NA
    habcv.part[87,'AIC'] = gammod$aic
    habcv.part[87,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[87,'R2'] = round(gamanova$r.sq,3);
    
    gf88 = formula(presence ~ s(log(Patch),k=k.num) + factor(HB))
    gammod = gam(gf88, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[88,'bathyclasses'] = NA
    habcv.part[88,'bpibroadclasses'] = NA
    habcv.part[88,'sedclasses'] = NA
    habcv.part[88,'bpifineclasses'] = NA
    habcv.part[88,'HB'] = gamanova$pTerms.table[1,3];
    habcv.part[88,'Patch'] = gamanova$s.table[1,4]
    habcv.part[88,'PatchCompact'] = NA
    habcv.part[88,'AIC'] = gammod$aic
    habcv.part[88,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[88,'R2'] = round(gamanova$r.sq,3);
    
    gf89 = formula(presence ~ s(PatchCompact,k=k.num) + factor(HB))
    gammod = gam(gf89, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[89,'bathyclasses'] = NA
    habcv.part[89,'bpibroadclasses'] = NA
    habcv.part[89,'sedclasses'] = NA
    habcv.part[89,'bpifineclasses'] = NA
    habcv.part[89,'HB'] = gamanova$pTerms.table[1,3];
    habcv.part[89,'Patch'] = NA
    habcv.part[89,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[89,'AIC'] = gammod$aic
    habcv.part[89,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[89,'R2'] = round(gamanova$r.sq,3);
    
    ### Fine included
    
    gf90 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf90, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[90,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[90,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[90,'sedclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[90,'bpifineclasses'] = gamanova$pTerms.table[4,3]
    habcv.part[90,'HB'] = gamanova$pTerms.table[5,3];
    habcv.part[90,'Patch'] = gamanova$s.table[1,4]
    habcv.part[90,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[90,'AIC'] = gammod$aic
    habcv.part[90,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[90,'R2'] = round(gamanova$r.sq,3);
    
    gf91 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf91, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[91,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[91,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[91,'sedclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[91,'bpifineclasses'] = gamanova$pTerms.table[4,3]
    habcv.part[91,'HB'] = gamanova$pTerms.table[5,3];
    habcv.part[91,'Patch'] = gamanova$s.table[1,4]
    habcv.part[91,'PatchCompact'] = NA
    habcv.part[91,'AIC'] = gammod$aic
    habcv.part[91,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[91,'R2'] = round(gamanova$r.sq,3);
    
    gf92 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf92, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[92,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[92,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[92,'sedclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[92,'bpifineclasses'] = gamanova$pTerms.table[4,3]
    habcv.part[92,'HB'] = gamanova$pTerms.table[5,3];
    habcv.part[92,'Patch'] = NA
    habcv.part[92,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[92,'AIC'] = gammod$aic
    habcv.part[92,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[92,'R2'] = round(gamanova$r.sq,3);
    
    gf93 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB)) 
    gammod = gam(gf93, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[93,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[93,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[93,'sedclasses'] = NA
    habcv.part[93,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[93,'HB'] = gamanova$pTerms.table[4,3];
    habcv.part[93,'Patch'] = gamanova$s.table[1,4]
    habcv.part[93,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[93,'AIC'] = gammod$aic
    habcv.part[93,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[93,'R2'] = round(gamanova$r.sq,3);
    
    gf94 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf94, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[94,'bathyclasses'] = NA
    habcv.part[94,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[94,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[94,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[94,'HB'] = gamanova$pTerms.table[4,3];
    habcv.part[94,'Patch'] = gamanova$s.table[1,4]
    habcv.part[94,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[94,'AIC'] = gammod$aic
    habcv.part[94,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[94,'R2'] = round(gamanova$r.sq,3);
    
    gf95 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(sedclasses) + factor(bpifineclasses) + factor(HB)) 
    gammod = gam(gf95, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[95,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[95,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[95,'sedclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[95,'bpifineclasses'] = gamanova$pTerms.table[4,3]
    habcv.part[95,'HB'] = gamanova$pTerms.table[5,3];
    habcv.part[95,'Patch'] = NA
    habcv.part[95,'PatchCompact'] = NA
    habcv.part[95,'AIC'] = gammod$aic
    habcv.part[95,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[95,'R2'] = round(gamanova$r.sq,3);
    
    gf96 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses) + factor(HB)) 
    gammod = gam(gf96, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[96,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[96,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[96,'sedclasses'] = NA
    habcv.part[96,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[96,'HB'] = gamanova$pTerms.table[4,3];
    habcv.part[96,'Patch'] = gamanova$s.table[1,4]
    habcv.part[96,'PatchCompact'] = NA
    habcv.part[96,'AIC'] = gammod$aic
    habcv.part[96,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[96,'R2'] = round(gamanova$r.sq,3);
    
    gf97 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf97, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[97,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[97,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[97,'sedclasses'] = NA
    habcv.part[97,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[97,'HB'] = gamanova$pTerms.table[4,3];
    habcv.part[97,'Patch'] = NA
    habcv.part[97,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[97,'AIC'] = gammod$aic
    habcv.part[97,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[97,'R2'] = round(gamanova$r.sq,3);
    
    gf98 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf98, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[98,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[98,'bpibroadclasses'] = NA
    habcv.part[98,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[98,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[98,'HB'] = gamanova$pTerms.table[4,3];
    habcv.part[98,'Patch'] = gamanova$s.table[1,4]
    habcv.part[98,'PatchCompact'] = NA
    habcv.part[98,'AIC'] = gammod$aic
    habcv.part[98,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[98,'R2'] = round(gamanova$r.sq,3);
    
    gf99 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf99, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[99,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[99,'bpibroadclasses'] = NA
    habcv.part[99,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[99,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[99,'HB'] = gamanova$pTerms.table[4,3];
    habcv.part[99,'Patch'] = NA
    habcv.part[99,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[99,'AIC'] = gammod$aic
    habcv.part[99,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[99,'R2'] = round(gamanova$r.sq,3);
    
    gf100 = formula(presence ~ factor(bathyclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf100, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[100,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[100,'bpibroadclasses'] = NA
    habcv.part[100,'sedclasses'] = NA
    habcv.part[100,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[100,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[100,'Patch'] = gamanova$s.table[1,4]
    habcv.part[100,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[100,'AIC'] = gammod$aic
    habcv.part[100,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[100,'R2'] = round(gamanova$r.sq,3);
    
    gf101 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf101, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[101,'bathyclasses'] = NA
    habcv.part[101,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[101,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[101,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[101,'HB'] = gamanova$pTerms.table[4,3];
    habcv.part[101,'Patch'] = gamanova$s.table[1,4]
    habcv.part[101,'PatchCompact'] = NA
    habcv.part[101,'AIC'] = gammod$aic
    habcv.part[101,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[101,'R2'] = round(gamanova$r.sq,3);
    
    gf102 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf102, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[102,'bathyclasses'] = NA
    habcv.part[102,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[102,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[102,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[102,'HB'] = gamanova$pTerms.table[4,3];
    habcv.part[102,'Patch'] = NA
    habcv.part[102,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[102,'AIC'] = gammod$aic
    habcv.part[102,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[102,'R2'] = round(gamanova$r.sq,3);
    
    gf103 = formula(presence ~ factor(bpibroadclasses) + s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB)) 
    gammod = gam(gf103, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[103,'bathyclasses'] = NA
    habcv.part[103,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[103,'sedclasses'] = NA
    habcv.part[103,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[103,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[103,'Patch'] = gamanova$s.table[1,4]
    habcv.part[103,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[103,'AIC'] = gammod$aic
    habcv.part[103,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[103,'R2'] = round(gamanova$r.sq,3);
    
    gf104 = formula(presence ~ factor(bathyclasses) + factor(bpibroadclasses) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf104, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[104,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[104,'bpibroadclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[104,'sedclasses'] = NA
    habcv.part[104,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[104,'HB'] = gamanova$pTerms.table[4,3];
    habcv.part[104,'Patch'] = NA
    habcv.part[104,'PatchCompact'] = NA
    habcv.part[104,'AIC'] = gammod$aic
    habcv.part[104,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[104,'R2'] = round(gamanova$r.sq,3);
    
    gf105 = formula(presence ~ factor(bathyclasses) + factor(sedclasses) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf105, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[105,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[105,'bpibroadclasses'] = NA
    habcv.part[105,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[105,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[105,'HB'] = gamanova$pTerms.table[4,3];
    habcv.part[105,'Patch'] =NA
    habcv.part[105,'PatchCompact'] = NA
    habcv.part[105,'AIC'] = gammod$aic
    habcv.part[105,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[105,'R2'] = round(gamanova$r.sq,3);
    
    gf106 = formula(presence ~ factor(bathyclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf106, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[106,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[106,'bpibroadclasses'] = NA
    habcv.part[106,'sedclasses'] = NA
    habcv.part[106,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[106,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[106,'Patch'] = gamanova$s.table[1,4]
    habcv.part[106,'PatchCompact'] = NA
    habcv.part[106,'AIC'] = gammod$aic
    habcv.part[106,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[106,'R2'] = round(gamanova$r.sq,3);
    
    gf107 = formula(presence ~ factor(bathyclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf107, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[107,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[107,'bpibroadclasses'] = NA
    habcv.part[107,'sedclasses'] = NA
    habcv.part[107,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[107,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[107,'Patch'] =NA
    habcv.part[107,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[107,'AIC'] = gammod$aic
    habcv.part[107,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[107,'R2'] = round(gamanova$r.sq,3);
    
    gf108 = formula(presence ~ factor(bpibroadclasses) + factor(sedclasses) + factor(bpifineclasses) + factor(HB)) 
    gammod = gam(gf108, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[108,'bathyclasses'] = NA
    habcv.part[108,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[108,'sedclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[108,'bpifineclasses'] = gamanova$pTerms.table[3,3]
    habcv.part[108,'HB'] = gamanova$pTerms.table[4,3];
    habcv.part[108,'Patch'] = NA
    habcv.part[108,'PatchCompact'] = NA
    habcv.part[108,'AIC'] = gammod$aic
    habcv.part[108,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[108,'R2'] = round(gamanova$r.sq,3);
    
    gf109 = formula(presence ~ factor(bpibroadclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses) + factor(HB)) 
    gammod = gam(gf109, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[109,'bathyclasses'] = NA
    habcv.part[109,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[109,'sedclasses'] = NA
    habcv.part[109,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[109,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[109,'Patch'] = gamanova$s.table[1,4]
    habcv.part[109,'PatchCompact'] = NA
    habcv.part[109,'AIC'] = gammod$aic
    habcv.part[109,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[109,'R2'] = round(gamanova$r.sq,3);
    
    gf110 = formula(presence ~ factor(bpibroadclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf110, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[110,'bathyclasses'] = NA
    habcv.part[110,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[110,'sedclasses'] = NA
    habcv.part[110,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[110,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[110,'Patch'] = NA
    habcv.part[110,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[110,'AIC'] = gammod$aic
    habcv.part[110,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[110,'R2'] = round(gamanova$r.sq,3);
    
    gf111 = formula(presence ~ factor(sedclasses) + s(log(Patch),k=k.num) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf111, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[111,'bathyclasses'] = NA
    habcv.part[111,'bpibroadclasses'] = NA
    habcv.part[111,'sedclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[111,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[111,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[111,'Patch'] = gamanova$s.table[1,4]
    habcv.part[111,'PatchCompact'] = NA
    habcv.part[111,'AIC'] = gammod$aic
    habcv.part[111,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[111,'R2'] = round(gamanova$r.sq,3);
    
    gf112 = formula(presence ~ factor(sedclasses) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf112, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[112,'bathyclasses'] = NA
    habcv.part[112,'bpibroadclasses'] = NA
    habcv.part[112,'sedclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[112,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[112,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[112,'Patch'] = NA
    habcv.part[112,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[112,'AIC'] = gammod$aic
    habcv.part[112,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[112,'R2'] = round(gamanova$r.sq,3);
    
    gf113 = formula(presence ~ s(log(Patch),k=k.num) + s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf113, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[113,'bathyclasses'] = NA
    habcv.part[113,'bpibroadclasses'] = NA
    habcv.part[113,'sedclasses'] = NA
    habcv.part[113,'bpifineclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[113,'HB'] = gamanova$pTerms.table[2,3];
    habcv.part[113,'Patch'] = gamanova$s.table[1,4]
    habcv.part[113,'PatchCompact'] = gamanova$s.table[2,4]
    habcv.part[113,'AIC'] = gammod$aic
    habcv.part[113,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[113,'R2'] = round(gamanova$r.sq,3);
    
    gf114 = formula(presence ~ factor(bathyclasses) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf114, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[114,'bathyclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[114,'bpibroadclasses'] = NA
    habcv.part[114,'sedclasses'] = NA
    habcv.part[114,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[114,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[114,'Patch'] = NA
    habcv.part[114,'PatchCompact'] = NA
    habcv.part[114,'AIC'] = gammod$aic
    habcv.part[114,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[114,'R2'] = round(gamanova$r.sq,3);
    
    gf115 = formula(presence ~ factor(bpibroadclasses) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf115, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[115,'bathyclasses'] = NA
    habcv.part[115,'bpibroadclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[115,'sedclasses'] = NA
    habcv.part[115,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[115,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[115,'Patch'] = NA
    habcv.part[115,'PatchCompact'] = NA
    habcv.part[115,'AIC'] = gammod$aic
    habcv.part[115,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[115,'R2'] = round(gamanova$r.sq,3);
    
    gf116 = formula(presence ~ factor(sedclasses) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf116, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[116,'bathyclasses'] = NA
    habcv.part[116,'bpibroadclasses'] = NA
    habcv.part[116,'sedclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[116,'bpifineclasses'] = gamanova$pTerms.table[2,3]
    habcv.part[116,'HB'] = gamanova$pTerms.table[3,3];
    habcv.part[116,'Patch'] = NA
    habcv.part[116,'PatchCompact'] = NA
    habcv.part[116,'AIC'] = gammod$aic
    habcv.part[116,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[116,'R2'] = round(gamanova$r.sq,3);
    
    gf117 = formula(presence ~ s(log(Patch),k=k.num) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf117, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[117,'bathyclasses'] = NA
    habcv.part[117,'bpibroadclasses'] = NA
    habcv.part[117,'sedclasses'] = NA
    habcv.part[117,'bpifineclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[117,'HB'] = gamanova$pTerms.table[2,3];
    habcv.part[117,'Patch'] = gamanova$s.table[1,4]
    habcv.part[117,'PatchCompact'] = NA
    habcv.part[117,'AIC'] = gammod$aic
    habcv.part[117,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[117,'R2'] = round(gamanova$r.sq,3);
    
    gf118 = formula(presence ~ s(PatchCompact,k=k.num) + factor(bpifineclasses) + factor(HB))
    gammod = gam(gf118, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[118,'bathyclasses'] = NA
    habcv.part[118,'bpibroadclasses'] = NA
    habcv.part[118,'sedclasses'] = NA
    habcv.part[118,'bpifineclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[118,'HB'] = gamanova$pTerms.table[2,3];
    habcv.part[118,'Patch'] = NA
    habcv.part[118,'PatchCompact'] = gamanova$s.table[1,4]
    habcv.part[118,'AIC'] = gammod$aic
    habcv.part[118,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[118,'R2'] = round(gamanova$r.sq,3);
    
    gf119 = formula(presence ~ factor(bpifineclasses) + factor(HB))
    gammod = gam(gf119, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    gamanova = anova(gammod);
    habcv.part[119,'bathyclasses'] = NA
    habcv.part[119,'bpibroadclasses'] = NA
    habcv.part[119,'sedclasses'] = NA
    habcv.part[119,'bpifineclasses'] = gamanova$pTerms.table[1,3]
    habcv.part[119,'HB'] = gamanova$pTerms.table[2,3];
    habcv.part[119,'Patch'] = NA
    habcv.part[119,'PatchCompact'] = NA
    habcv.part[119,'AIC'] = gammod$aic
    habcv.part[119,'devexp'] = round(gamanova$dev.expl,3)
    habcv.part[119,'R2'] = round(gamanova$r.sq,3);
    for (i in 2:7){
       habcv.part[which(habcv.part[,i]<.05),i] = '*'
       habcv.part[which(habcv.part[,i]<.01),i] = '**'
       habcv.part[which(habcv.part[,i]<.001),i] = '***' 
    }
    
    #write.csv(habcv.part,file=paste('C:/Users/brian.grieve/Documents/coca/Data/cval',Season,'/spp',spp,'/habCVpartHB_',part,'.csv',sep=''),row.names = FALSE)
    
   
    bestmod = which(habcv.part$AIC==min(habcv.part$AIC))
    # check for equal AIC values. This never occurred. 
    if (length(bestmod)==1){gf = gflist[[bestmod]]; habcval.df[part,'modelnumber'] = bestmod;}
    if (length(bestmod)>1){gf = gflist[[bestmod[1]]]; habcval.df[part,'modelnumber'] = bestmod[1]; warning('multiple models')}
    
    # refit best model
    gammod = gam(gf, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
    # predict witheld data
    testpred = predict.gam(gammod, newdata=fishtest, type = 'response');
    if (length(which(is.na(testpred)))>0){
      fishtest = fishtest[-which(is.na(testpred)),]    
      testpred = testpred[-which(is.na(testpred))]
    }

    habcval.mod[[part]] = gammod;
    habcval.df[part,'AIC'] = gammod$aic;
    habcval.df[part,'dev.exp'] = summary(gammod)$dev.exp;
    habcval.df[part,'r.sq'] = summary(gammod)$r.sq;
   
    
    # Pull out temperature coefficients
    if (Season==2){jlcoefs = read.csv(paste('C:/Users/brian.grieve/Documents/coca/Rich_GAMS/tmb/springcoefs/tmbcoefs',spp,'.csv',sep=''),header=TRUE);}
    if (Season==4){jlcoefs = read.csv(paste('C:/Users/brian.grieve/Documents/coca/Rich_GAMS/tmb/fallcoefs/tmbcoefs',spp,'.csv',sep=''),header=TRUE);}
    
    Er = jlcoefs[part,'Er']; 
    Ed = jlcoefs[part,'Ed']; 
    cc = jlcoefs[part,'cc']; 
    Topt.c = jlcoefs[part,'Topt.c']; 
    ccexp = jlcoefs[part,'ccexp']; 
    k = jlcoefs[part,'k'];
        
    # predict abundance with Johnson Lewin equation and a preliminary combination of the models
    fishtest$jlpred = B_Afunction(T.c=fishtest$BTEMP, Topt.c=Topt.c, Er=Er, Ed=Ed, cc=cc)
    fishtest$jlpred = sqrt(fishtest$jlpred/max(fishtest$jlpred,na.rm=TRUE))
    fishtest$multipred = fishtest$jlpred*sqrt(testpred);
    
    PAdata = cbind(1:nrow(fishtest),fishtest$presence,as.numeric(fishtest$multipred))
    thresh = optimal.thresholds(PAdata,na.rm = TRUE)
    PAdata = cbind(PAdata,as.numeric(fishtest$multipred>=thresh[3,2])) #optimal threshold to maximize TSS
    CMX4 = cmx(PAdata,which.model = 2);
    
    
    habcval.df[part,'auc'] = PresenceAbsence::auc(PAdata,which.model = 1)[[1]]
    habcval.df[part,'maxtss_thresh'] = thresh[3,2]
    habcval.df[part,'maxtss_sens'] = sensitivity(CMX4)[1]
    habcval.df[part,'maxtss_spec'] = specificity(CMX4)[1]
    habcval.df[part,'mult_cor'] = cor(fishtest$multipred,fishtest$num_sa_night)
    habcval.df[part,'maxtss_cor'] = cor(PAdata[,4]*fishtest$jlpred,fishtest$num_sa_night)
    habcval.df[part,'mult_rmse'] = as.numeric(postResample(fishtest$multipred,fishtest$num_sa_night)[1])
    habcval.df[part,'maxtss_rmse'] = as.numeric(postResample(PAdata[,4]*fishtest$jlpred,fishtest$num_sa_night)[1])

    
    }
  write.csv(habcval.df,file=paste('C:/Users/brian.grieve/Documents/coca/Data/cval',Season,'/spp',spp,'/habCValDFHB',spp,'.csv',sep=''),row.names = FALSE)
}
}
# Summary data frame of all 25 species
habmod.all = data.frame(matrix(data=NA,nrow = 25,ncol=ncol(habcval.df)+1)); names(habmod.all) = c('Commonname',names(habcval.df))
for (Season in c(2,4)){
  for (spp in u_svspp[c(1:25)]){
    ind.vec = which(u_svspp==spp);
    habcval.df = read.csv(file=paste('C:/Users/brian.grieve/Documents/coca/Data/cval',Season,'/spp',spp,'/habCValDFHB',spp,'.csv',sep=''))
    habmod.all[ind.vec,1] = species_conversions[ind.vec,'Common name']
    habmod.all[ind.vec,2:21] = habcval.df[which(habcval.df$AIC==min(habcval.df$AIC,na.rm=TRUE)),]
    habmod.all[,7:21] = round(habmod.all[,7:21],3)
    habmod.all$tss_max = habmod.all$maxtss_sens+habmod.all$maxtss_spec-1;
    habmod.all$tss_prev = habmod.all$prev_sens+habmod.all$prev_spec-1;
  }
    write.csv(habmod.all,file=paste('C:/Users/brian.grieve/Documents/coca/Data/cval',Season,'/habmodallHB',Season,'.csv',sep=''),row.names=FALSE)
} 
    
# Threshold determination and model combination evaluation
thresh.df = data.frame(matrix(data=0,nrow=25,ncol=10)); 
names(thresh.df) = c('max_thresh_nn','max_tss_nn','max_thresh_22','max_tss_22','max_thresh_nsd','max_tss_nsd','max_thresh_2sd','max_tss_2sd','max_thresh_22sd','max_tss_22sd')
for (spp in u_svspp[c(1:25)]){
    # load CVal
  habcval.df = read.csv(paste('C:/Users/brian.grieve/Documents/coca/Data/cval',Season,'/spp',spp,'/habCValDFHB',spp,'.csv',sep=''),header = TRUE)
  spname = (species_conversions[which(species_conversions$SVSPP==spp),'Common name']); print(spname)
  ind.vec = which(u_svspp==spp)
  fishdata = alltrawl[which(alltrawl$SVSPP==spp & alltrawl$SEASON==Season),];
  fishdata[which(fishdata$sedclasses==128),'sedclasses'] = NA;
  if(patchcor[ind.vec,'usepatch']=='Perim'){fishdata$Patch = fishdata$PatchPerim; habvars2$Patch = habvars2$PatchPerim}
  if(patchcor[ind.vec,'usepatch']=='Area'){fishdata$Patch = fishdata$PatchArea; habvars2$Patch = habvars2$PatchArea}
  i.na = which(is.na(fishdata[,'presence'])|is.na(fishdata$BTEMP)|is.na(fishdata$PatchArea)|is.na(fishdata$PatchPerim)|is.na(fishdata$PatchCompact)|is.na(fishdata$bathyclasses)|is.na(fishdata$bpibroadclasses)|is.na(fishdata$sedclasses)|is.na(fishdata$HB))
  if(length(i.na)>0){
    fishdata = fishdata[-i.na,]}
  
  part = which(habcval.df$AIC==min(habcval.df$AIC,na.rm=TRUE));
  fishtrain = fishdata[which(fishdata$partition!=part),]
  fishtest = fishdata[which(fishdata$partition==part),]
  # final runs use all data
  fishtest = fishdata; fishtrain = fishdata;
  
  if (Season==2){jlcoefs = read.csv(paste('C:/Users/brian.grieve/Documents/coca/Rich_GAMS/tmb/springcoefs/tmbcoefs',spp,'.csv',sep=''),header=TRUE);}
  if (Season==4){jlcoefs = read.csv(paste('C:/Users/brian.grieve/Documents/coca/Rich_GAMS/tmb/fallcoefs/tmbcoefs',spp,'.csv',sep=''),header=TRUE);}
  Er = jlcoefs[part,'Er']; 
  Ed = jlcoefs[part,'Ed']; 
  cc = jlcoefs[part,'cc']; 
  Topt.c = jlcoefs[part,'Topt.c']; 
  ccexp = jlcoefs[part,'ccexp']; 
  k = jlcoefs[part,'k']
  
  bestmod = habcval.df[part,'modelnumber']
  gf = gflist[[bestmod]]
  # Using all data for final model
  gammod = gam(gf, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
  
  # Determine what combination of models works best. Candidates were normalized, square rooted, or normalized to 2SD
  testpred = predict.gam(gammod, newdata=fishtest, type = 'response');
  maxsd = mean(testpred,na.rm = TRUE) + 2*sd(testpred,na.rm = TRUE)
  testpred_norm = testpred/max(testpred,na.rm = TRUE) 
  testpred_sd = testpred
  testpred_sd[which(testpred_sd>maxsd)] = maxsd;
  testpred_sd = testpred_sd/maxsd
  
  
  maxabun = B_Afunction(T.c = Topt.c, Topt.c = Topt.c, Er=Er, Ed=Ed, cc=cc)
  fishtest$jlpred = B_Afunction(T.c=fishtest$BTEMP, Topt.c=Topt.c, Er=Er, Ed=Ed, cc=cc)
  fishtest$jlpred_norm = fishtest$jlpred/maxabun
  fishtest$jlpred_sqrt = sqrt(fishtest$jlpred/maxabun)
  
  fishtest$jlnorm_habnorm = fishtest$jlpred_norm*testpred_norm;
  fishtest$jlsqrt_habsqrt = fishtest$jlpred_sqrt*sqrt(testpred_norm)
  fishtest$jlnorm_habsd = fishtest$jlpred_norm*testpred_sd
  fishtest$jlsqrt_habsd = fishtest$jlpred_sqrt*testpred_sd
  fishtest$jlsqrt_habsdsqrt = fishtest$jlpred_sqrt*sqrt(testpred_sd)
  
  
  PAdata = cbind(1:nrow(fishtest),fishtest$presence,fishtest$jlnorm_habnorm,fishtest$jlsqrt_habsqrt,fishtest$jlnorm_habsd,fishtest$jlsqrt_habsd,fishtest$jlsqrt_habsdsqrt)
  thresh.df[ind.vec,'max_thresh_nn'] = optimal.thresholds(PAdata,which.model = 1)[3,2];
  thresh.df[ind.vec,'max_thresh_22'] = optimal.thresholds(PAdata,which.model = 2)[3,2];
  thresh.df[ind.vec,'max_thresh_nsd'] = optimal.thresholds(PAdata,which.model = 3)[3,2];
  thresh.df[ind.vec,'max_thresh_2sd'] = optimal.thresholds(PAdata,which.model = 4)[3,2];
  thresh.df[ind.vec,'max_thresh_22sd'] = optimal.thresholds(PAdata,which.model = 5)[3,2];

  
  for (mod in 1:5){
    cmx = cmx(PAdata, threshold = thresh.df[ind.vec,mod*2-1],which.model = mod)
    thresh.df[ind.vec,mod*2] = sensitivity(cmx,st.dev=FALSE) + specificity(cmx,st.dev=FALSE) -1}
  
}
# Determine which threshold is best
for (i in 1:25){
  thresh.df$which[i] = which(thresh.df[i,c(2,4,6,8,10)]==max(thresh.df[i,c(2,4,6,8,10)]))}

write.csv(thresh.df,file=paste('C:/Users/brian.grieve/Documents/coca/Data/cval',Season,'/threshdf',Season,'.csv',sep=''),row.names=FALSE,col.names = TRUE)

 


# Final evaluation tables of ideal combined model
eval.df = data.frame(matrix(data=0,nrow=25,ncol=13)); 
names(eval.df) = c('SVSPP','Season','modelnumber','dev.exp','r.sq','auc','tss','sens','spec','TP','TN','FP','FN')
for (spp in u_svspp[c(1:25)]){
  habcval.df = read.csv(paste('C:/Users/brian.grieve/Documents/coca/Data/cval',Season,'/spp',spp,'/habCValDFHB',spp,'.csv',sep=''),header = TRUE)
  spname = (species_conversions[which(species_conversions$SVSPP==spp),'Common name']); print(spname)
  ind.vec = which(u_svspp==spp)
  fishdata = alltrawl[which(alltrawl$SVSPP==spp & alltrawl$SEASON==Season),];
  fishdata[which(fishdata$sedclasses==128),'sedclasses'] = NA;
  if(patchcor[ind.vec,'usepatch']=='Perim'){fishdata$Patch = fishdata$PatchPerim; habvars2$Patch = habvars2$PatchPerim}
  if(patchcor[ind.vec,'usepatch']=='Area'){fishdata$Patch = fishdata$PatchArea; habvars2$Patch = habvars2$PatchArea}
   i.na = which(is.na(fishdata[,'presence'])|is.na(fishdata$BTEMP)|is.na(fishdata$PatchArea)|is.na(fishdata$PatchPerim)|is.na(fishdata$PatchCompact)|is.na(fishdata$bathyclasses)|is.na(fishdata$bpibroadclasses)|is.na(fishdata$sedclasses)|is.na(fishdata$HB))
  if(length(i.na)>0){
    fishdata = fishdata[-i.na,]}
  
  part = which(habcval.df$AIC==min(habcval.df$AIC,na.rm=TRUE));
  fishtrain = fishdata[which(fishdata$partition!=part),]
  fishtest = fishdata[which(fishdata$partition==part),]
  # final runs use all data
  fishtest = fishdata; fishtrain = fishdata;
  
  if (Season==2){jlcoefs = read.csv(paste('C:/Users/brian.grieve/Documents/coca/Rich_GAMS/tmb/springcoefs/tmbcoefs',spp,'.csv',sep=''),header=TRUE);}
  if (Season==4){jlcoefs = read.csv(paste('C:/Users/brian.grieve/Documents/coca/Rich_GAMS/tmb/fallcoefs/tmbcoefs',spp,'.csv',sep=''),header=TRUE);}
  Er = jlcoefs[part,'Er']; 
  Ed = jlcoefs[part,'Ed']; 
  cc = jlcoefs[part,'cc']; 
  Topt.c = jlcoefs[part,'Topt.c']; 
  ccexp = jlcoefs[part,'ccexp']; 
  k = jlcoefs[part,'k']
  
  bestmod = habcval.df[part,'modelnumber']
  gf = gflist[[bestmod]]
  # Using all data for final model
  gammod = gam(gf, data = fishtrain, family = 'binomial', select = TRUE, method = 'REML', gamma = 1.4)
  testpred = predict.gam(gammod, newdata=fishtest, type = 'response');
  maxsd = mean(testpred,na.rm = TRUE) + 2*sd(testpred,na.rm = TRUE)
  testpred_norm = testpred/max(testpred,na.rm = TRUE) 
  testpred_sd = testpred
  testpred_sd[which(testpred_sd>maxsd)] = maxsd;
  testpred_sd = testpred_sd/maxsd
  
  
  maxabun = B_Afunction(T.c = Topt.c, Topt.c = Topt.c, Er=Er, Ed=Ed, cc=cc)
  fishtest$jlpred = B_Afunction(T.c=fishtest$BTEMP, Topt.c=Topt.c, Er=Er, Ed=Ed, cc=cc)
  fishtest$jlpred_sqrt = sqrt(fishtest$jlpred/maxabun)
  fishtest$jlsqrt_habsd = fishtest$jlpred_sqrt*testpred_sd

  
  
  PAdata = cbind(1:nrow(fishtest),fishtest$presence,fishtest$jlsqrt_habsd)
  thresh = optimal.thresholds(PAdata,which.model = 1)[3,2];
  CMX = cmx(PAdata, threshold = thresh)
  
  eval.df[ind.vec,'SVSPP'] = spp
  eval.df[ind.vec,'Season'] = habcval.df[part,'Season']
  eval.df[ind.vec,'modelnumber'] = habcval.df[part,'modelnumber']
  eval.df[ind.vec,'dev.exp'] = habcval.df[part,'dev.exp']
  eval.df[ind.vec,'r.sq'] = habcval.df[part,'r.sq']
  eval.df[ind.vec,'auc'] = habcval.df[part,'auc']
  eval.df[ind.vec,'tss'] = sensitivity(CMX, st.dev = FALSE) + specificity(CMX, st.dev = FALSE) - 1
  eval.df[ind.vec,'sens'] = sensitivity(CMX, st.dev = FALSE)
  eval.df[ind.vec,'spec'] = specificity(CMX, st.dev = FALSE)
  eval.df[ind.vec,'TP'] = CMX[1,1]
  eval.df[ind.vec,'TN'] = CMX[2,2]
  eval.df[ind.vec,'FP'] = CMX[1,2]
  eval.df[ind.vec,'FN'] = CMX[2,1]
    
}
eval.df[,c(4:9)] = round(eval.df[,c(4:9)],3) 

write.csv(eval.df,file=paste('C:/Users/brian.grieve/Documents/coca/Data/cval',Season,'/evaldf',Season,'.csv',sep=''),row.names=FALSE,col.names = TRUE)


# Project and map final habitat models, save GAMs
PAgams = list(); 

for (spp in u_svspp){
  spname = (species_conversions[[which(species_conversions$SVSPP==spp),'Common name']]); print(spname)
  ind.vec = which(u_svspp==spp)
  ind.fish = which(alltrawl$SVSPP==spp & alltrawl$SEASON==Season)
  fishdata = alltrawl[ind.fish,];
  fishdata[which(fishdata$sedclasses==128),'sedclasses'] = NA;
  if(patchcor[ind.vec,'usepatch']=='Perim'){fishdata$Patch = fishdata$PatchPerim; habvars2$Patch = habvars2$PatchPerim}
  if(patchcor[ind.vec,'usepatch']=='Area'){fishdata$Patch = fishdata$PatchArea; habvars2$Patch = habvars2$PatchArea}
  
  i.na = which(is.na(fishdata[,resp])|is.na(fishdata$BTEMP)|is.na(fishdata$PatchArea)|is.na(fishdata$PatchPerim)|is.na(fishdata$PatchCompact)|is.na(fishdata$bathyclasses)|is.na(fishdata$bpibroadclasses)|is.na(fishdata$bpifineclasses)|is.na(fishdata$sedclasses)|is.na(fishdata$HB))
  if(length(i.na)>0){
    fishdata = fishdata[-i.na,]}

  habmod.part = habmod.all[ind.vec,'partition']
  habmod.gf = gflist[[habmod.all[ind.vec,'modelnumber']]]
  fishtrain = fishdata[which(fishdata$partition!=habmod.part),]
  fishtest = fishdata[which(fishdata$partition==habmod.part),]
  
  # final run uses all data
  gampa = gam(habmod.gf, data = fishdata, family = binomial, select = TRUE, method = 'REML', gamma = 1.4)
  fishdata$habpred_raw = predict.gam(gampa,fishdata,type='response')

  # Save GAM
  PAgams[[ind.vec]] = gampa;
  
  # Save Habitat Projection
  papred = predict.gam(gampa, newdata=habvars2, type = 'response');
  papred = matrix(data=papred,nrow=9869,ncol=10704);
  papred[which(is.na(papred))] = -999;
  write.table(papred, file= paste('C:/Users/brian.grieve/Documents/coca/Benthic_maps/s',Season,'/binomial/pred',spp,'.txt',sep='') , row.names = FALSE, col.names = FALSE)
  
  rm(gampa,papred); gc()
}
save(PAgams,file = paste('C:/Users/brian.grieve/Documents/coca/Data/FinalGAMs',Season,'.Rdata',sep=''))

