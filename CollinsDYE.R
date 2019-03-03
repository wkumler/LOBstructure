#CollinsDYE.R
##Attempting to run DYEatom 6d_TLE_ESI cruise data through LOBSTAHS
##using settings from PrepOrbidata in the VM toolbox

library(LOBSTAHS)
library(BiocParallel)



#Setup and variable declaration

pol <- "negative"
sim_pol <- unlist(strsplit(pol, ""))[1:3]
sim_pol <- paste0(toupper(sim_pol[1]), sim_pol[2], sim_pol[3], collapse = "")

filesource <- 
paste0("/global/scratch/wkumler/DYEatom_data/6a_TLE_ESI/mzXML_", 
tolower(sim_pol))
filesource

##Parallel stuff
ncores <- detectCores()-1
param <- MulticoreParam(workers = ncores)
register(param)


##List of filenames
mzXML_files <- list.files(filesource, full.names = T, recursive = T)
mzXML_files <- mzXML_files[1:8]

#Start timing
tmx <- Sys.time()





#XCMS

##Create xset object
xset <- xcmsSet(mzXML_files,
                method="centWave",
                BPPARAM = bpparam(),
                prefilter=c(3,7500), 
                profparam=list(step=0.001),
                ppm=2.5,
                peakwidth=c(10,45),
                fitgauss=TRUE,
                noise=500,
                mzdiff=0.005,
                snthresh=10,
                integrate=1,
                mzCenterFun=c("wMean"))
xset
Sys.time()-tmx

##Grouping
xset <- group(xset,
              method="density",
              bw = 5, #15?
              minfrac = 0.25,
              minsamp = 2,
              mzwid = 0.015, #0.001?
              max = 50)
Sys.time()-tmx
xset <- retcor(xset,
               method = "obiwarp",
               plottype = "none",
               profStep = 1,
               response = 1,
               distFunc = "cor_opt",
               gapInit = NULL,
               gapExtend = NULL,
               factorDiag = 2,
               factorGap = 1,
               localAlignment = 0,
               initPenalty = 0)
Sys.time()-tmx
xset <- group(xset,
              method="density",
              bw=5,
              minfrac=0.25,
              minsamp = 2,
              mzwid=0.015,
              max=50)
xset
Sys.time()-tmx

##Fill peaks
register(MulticoreParam(workers=2))
xset.fin <- fillPeaks.chrom(xset)

Sys.time()-tmx



#CAMERA

xset_a = annotate(xset.fin, 
                  polarity=pol,
                  quick=F,
                  sample=NA,
                  sigma=6,
                  perfwhm=0.6,
                  cor_eic_th=0.75,
                  graphMethod="hcs",
                  pval=0.05,
                  calcCiS=T,
                  calcIso=T,
                  calcCaS=F,
                  maxcharge=4,
                  maxiso=4,
                  minfrac=0.5,
                  psg_list=NULL,
                  rules=NULL,
                  multiplier=3,
                  max_peaks=100,
                  intval="into",
                  ppm=2.5,
                  mzabs=0.0015)
Sys.time()-tmx

#LOBSTAHS

##Load database
data(default.LOBdbase)

##Screen peaks
LOB <- doLOBscreen(xsA=xset_a, match.ppm = 2.5, 
rt.restrict=T, exclude.oddFA=F)
Sys.time()-tmx

##Write file out
write.csv(LOB@peakdata, file = paste0("CollinsDYE_", sim_pol, 
".csv"))
