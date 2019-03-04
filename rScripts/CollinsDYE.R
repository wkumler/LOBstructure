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

#Start timing
tmx <- Sys.time()





#XCMS

##Create xset object
xset_method = "centWave"
xset_BPPARAM = bpparam()
xset_prefilter=c(3,7500) 
xset_profparam=list(step=0.001)
xset_ppm=2.5
xset_peakwidth=c(10,45)
xset_fitgauss=TRUE
xset_noise=500
xset_mzdiff=0.005
xset_snthresh=10
xset_integrate=1
xset_mzCenterFun=c("wMean")
xset <- xcmsSet(mzXML_files,
                method=xset_method,
                BPPARAM = xset_BPPARAM,
                prefilter=xset_prefilter, 
                profparam=xset_profparam,
                ppm=xset_ppm,
                peakwidth=xset_peakwidth,
                fitgauss=xset_fitgauss,
                noise=xset_noise,
                mzdiff=xset_mzdiff,
                snthresh=xset_snthresh,
                integrate=xset_integrate,
                mzCenterFun=xset_mzCenterFun)
xset
xset_time <- Sys.time()-tmx

##Grouping
group_method = "density"
group_bw <- 5 #15? Decreased from default = 30
group_minfrac = 0.25 #Decreased from default = 0.5
group_minsamp = 2 #Increased from default = 1
group_mzwid = 0.015 #0.001? Decreased from default = 0.25
group_max = 50 #Default

rc_method = "obiwarp"
rc_plottype = "none"
rc_profStep = 1 #Default
rc_response = 1 #Default
rc_distFunc = "cor_opt" #Default
rc_gapInit = NULL #Default
rc_gapExtend = NULL #Default
rc_factorDiag = 2 #Default
rc_factorGap = 1 #Default
rc_localAlignment = 0 #Default
rc_initPenalty = 0 #Default

xset <- group(xset,
              method=group_method,
              bw = group_bw,
              minfrac = group_minfrac,
              minsamp = group_minsamp,
              mzwid = group_mzwid,
              max = group_max)
g1_time <- Sys.time()-tmx
xset <- retcor(xset,
               method = rc_method,
               plottype = rc_plottype,
               profStep = rc_profStep,
               response = rc_response,
               distFunc = rc_distFunc,
               gapInit = rc_gapInit,
               gapExtend = rc_gapExtend,
               factorDiag = rc_factorDiag,
               factorGap = rc_factorGap,
               localAlignment = rc_localAlignment,
               initPenalty = rc_initPenalty)
rc_time <- Sys.time()-tmx
xset <- group(xset,
              method=group_method,
              bw = group_bw,
              minfrac = group_minfrac,
              minsamp = group_minsamp,
              mzwid = group_mzwid,
              max = group_max)
xset
g2_time <- Sys.time()-tmx

##Fill peaks
register(MulticoreParam(workers=2))
xset.fin <- fillPeaks.chrom(xset)

Sys.time()-tmx



#CAMERA

cm_polarity = pol
cm_quick = F #Default
cm_sample = NA #Default
cm_sigma = 6 #Default
cm_perfwhm = 0.6 #Default
cm_cor_eic_th = 0.75 #Default
cm_graphMethod = "hcs" #Default
cm_pval=0.05 #Default
cm_calcCiS=T #Default
cm_calcIso=T #Changed from default FALSE to TRUE
cm_calcCaS=F #Default
cm_maxcharge=4 #Increased from default 3 to 4
cm_maxiso=4 #Default
cm_minfrac=0.5 #Default
cm_psg_list=NULL #Default
cm_rules=NULL #Default
cm_multiplier=3 #Default
cm_max_peaks=100 #Default
cm_intval="into" #Default
cm_ppm=2.5 #Decreased from default 5
cm_mzabs=0.0015 #Decreased from default 0.015

xset_a = annotate(xset.fin, 
                  polarity=cm_polarity,
                  quick=cm_quick,
                  sample=cm_sample,
                  sigma=cm_sigma,
                  perfwhm=cm_perfwhm,
                  cor_eic_th=cm_cor_eic_th,
                  graphMethod=cm_graphMethod,
                  pval=cm_pval,
                  calcCiS=cm_calcCiS,
                  calcIso=cm_calcIso,
                  calcCaS=cm_calcCaS,
                  maxcharge=cm_maxcharge,
                  maxiso=cm_maxiso,
                  minfrac=cm_minfrac,
                  psg_list=cm_psg_list,
                  rules=cm_rules,
                  multiplier=cm_multiplier,
                  max_peaks=cm_max_peaks,
                  intval=cm_intval,
                  ppm=cm_ppm,
                  mzabs=cm_mzabs)
cm_time <- Sys.time()-tmx

#LOBSTAHS

##Load database
data(default.LOBdbase)

##Screen peaks
LOB <- doLOBscreen(xsA=xset_a, match.ppm = 2.5, 
rt.restrict=F, exclude.oddFA=F)
Sys.time()-tmx

##Write report file
report_name <- paste0("../Outputs/Collins_", Sys.getenv("SLURM_JOBID"))

sink(report_name)
cat("===========================")
cat("LOBSTAHS OUTPUT REPORT\n")
paste("FOR JOB", Sys.getenv("SLURM_JOBID"))
cat("===========================")

cat("XCMS SETTINGS:\n")
cat("Files:", mzXML_files)
cat("Method:", xset_method)
cat("BPPARAM:", xset_BPPARAM)
cat("Prefilter:", xset_prefilter)
cat("Profparam:", xset_profparam)
cat("ppm:", xset_ppm)
cat("Peakwidth:", xset_peakwidth)
cat("Fitgauss:", xset_fitgauss)
cat("Noise:", xset_noise)
cat("mzdiff:", xset_mzdiff)
cat("snthresh:", xset_snthresh)
cat("Integrate:", xset_integrate)
cat("mzCenterFun", xset_mzCenterFun, "\n\n")
cat("Time taken:", xset_time)

cat("\n\n\n\n\n")

cat("GROUPING SETTINGS:\n")
cat("Method:", group_method)
cat("BW:", group_bw)
cat("MinFrac:", group_minfrac)
cat("MinSamp:", group_minsamp)
cat("mzwid:", group_mzwid)
cat("Max:", group_max, "\n\n")
cat("Time taken \#1:", time_g1-xset_time)
cat("Time taken \#2:", time_g2-rc_time)

cat("\n\n\n\n\n")

cat("RETCOR SETTINGS:\n")
cat("Method:", rc_method)
cat("plottype:", rc_plottype)
cat("profStep:", rc_profStep)
cat("response:", rc_response)
cat("distFunc:", rc_distFunc)
cat("gapInit:", rc_gapInit)
cat("gapExtend:", rc_gapExtend)
cat("factorDiag:", rc_factorDiag)
cat("factorGap:", rc_factorGap)
cat("localAlignment:", rc_localAlignment)
cat("initPenalty:", rc_initPenalty, "\n\n")
cat("Time taken: ", rc_time-time_g1)

cat("\n\n\n\n\n")

cat("CAMERA SETTINGS:\n")
cat("polarity",cm_polarity)
cat("quick", cm_quick)
cat("sample",cm_sample)
cat("sigma",cm_sigma)
cat("perfwhm",cm_perfwhm)
cat("cor_eic_th", cm_cor_eic_th)
cat("graphMethod",cm_graphMethod)
cat("pval",cm_pval)
cat("calcCiS", cm_calcCiS)
cat("calcIso", cm_calcIso)
cat("calcCaS", cm_calcCaS)
cat("maxcharge", cm_maxcharge)
cat("maxiso", cm_maxiso)
cat("minfrac", cm_minfrac)
cat("psg_list", cm_psg_list)
cat("rules", cm_rules)
cat("multiplier", cm_multiplier)
cat("max_peaks", cm_max_peaks)
cat("intval", cm_intval)
cat("ppm", cm_ppm)
cat("mzabs", cm_mzabs, "\n")
cat("Time taken: ", cm_time-time_g2)

cat("\n\n\n\n\n")

cat("LOBSTAHS SETTINGS:\n")
cat("match.ppm", 2.5)
cat("rt.restrict", F)
cat("exclude.oddFA",F)

sink()
