#' JSQI
#' 
#' Signal quality index calculator.
#' Original MATLAB code by James Sun
#' Ported to R by Ran Liu on June 18, 2019
#' 
#' @param features Features computed using abpfeature()
#' @param onset Onsets computed using wabp()
#' @param abp ABP waveform (125Hz sampled, mmHg units of measure)
#' 
#' @return The following signal quality metrics:
#' 1. Logical OR of cols 2 through 10
#' 2. P not physiologic 
#' 3. MAP not physiologic
#' 4. HR not physiologic
#' 5. PP not physiologic
#' 6. abnormal Psys
#' 7. abnormal Pdias
#' 8. abnormal period
#' 9. abnormal P(onset)
#' 10. noisy beat
#' 
#' @export
jsqi = function(features, onset, abp) {

    rangeP = c(20, 300)
    rangeMAP = c(30, 200)
    rangeHR = c(20, 200)
    rangePP = c(20, Inf)

    dPsys = 20
    dPdias = 20
    dPeriod = 62.5
    dPOnset = 20

    noise = -3

    Psys = features[,2]
    Pdias = features[,4]
    PP = features[,5]
    MAP = features[,6]
    BeatPeriod = features[,7]
    mean_dyneg = features[,8]
    HR = 60*125/BeatPeriod

    badP = which(Pdias<rangeP[1]|Psys>rangeP[2])
    badMAP = which(MAP<rangeMAP[1]|MAP>rangeMAP[2])
    badHR = which(HR<rangeHR[1]|HR>rangeHR[2])
    badPP = which(PP<rangePP[1])

    jerkPSys = 1 + which(abs(Psys[-1]-Psys[-length(Psys)]) > dPsys)
    jerkPdias = which(abs(Pdias[-1]-Pdias[-length(Pdias)]) > dPdias)
    jerkPeriod = 1 + which(abs(BeatPeriod[-1]-BeatPeriod[-length(BeatPeriod)]) > dPeriod)
    POnset = abp[onset]
    jerkPOnset = which(abs(POnset[-1]-POnset[-length(POnset)]) > dPOnset)

    noisy = which(mean_dyneg < noise)

    bq = array(F,dim=c(length(onset),10))
    bq[badP,2] = T
    bq[badMAP,3] = T
    bq[badHR,4] = T
    bq[badPP,5] = T
    bq[jerkPSys,6] = T
    bq[jerkPdias,7] = T
    bq[jerkPeriod,8] = T
    bq[jerkPOnset,9] = T
    bq[noisy,10] = T

    bq[,1] = apply(bq[,-1],1,any)

    # make all "...101..." into "...111..."
    z = bq[,1]
    z = z[-1]-z[-length(z)]
    z = z[-1]-z[-length(z)]
    bq[which(z==2)+1,1] = T

    return(bq)

}