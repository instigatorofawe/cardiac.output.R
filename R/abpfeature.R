#' ABPFEATURE
#'
#' ABP waveform feature extractor.
#' Original MATLAB code by James Sun on Nov 19, 2005.
#' Ported to R by Ran Liu on June 24, 2018.
#'
#' @param abp Arterial blood pressure waveform
#' @param OnsetTimes Onset times (in indices)
#'
#' @return List of ABP features:
#' 1. Time of systole (samples)
#' 2. Systolic BP (mmHg)
#' 3. Time of diastole (samples)
#' 4. Diastolic BP (mmHg)
#' 5. Pulse Pressure (mmHg)
#' 6. Mean ABP (mmHg)
#' 7. Beat Period (samples)
#' 8. mean_dyneg TODO
#' 9. End of systole time, 0.3*sqrt(RR) method
#' 10. Area under systole, 0.3*sqrt(RR) method TODO
#' 11. End of systole time, 1st min-slope method TODO
#' 12. Area under systole, 1st min-slope method TODO
#' @export
abpfeature = function(abp, OnsetTimes) {
    Window = 40
    OT = OnsetTimes[-length(OnsetTimes)]
    BeatQty = length(OT)

    MinDomain = array(rep(0,BeatQty*Window),dim=c(BeatQty,Window))
    MaxDomain = array(rep(0,BeatQty*Window),dim=c(BeatQty,Window))

    for (i in 1:Window) {
        MinDomain[,i] = OT-i+1
        MaxDomain[,i] = OT+i-1
    }

    MinDomain[MinDomain<1] = 1
    MaxDomain[MaxDomain<1] = 1

    P_dias=apply(array(abp[MinDomain],dim=c(BeatQty,Window)),1,min)
    Dindex=apply(array(abp[MinDomain],dim=c(BeatQty,Window)),1,which.min)
    P_sys=apply(array(abp[MaxDomain],dim=c(BeatQty,Window)),1,max)
    Sindex=apply(array(abp[MaxDomain],dim=c(BeatQty,Window)),1,which.max)

    DiasTime = sapply(1:BeatQty, function(x) MinDomain[x,Dindex[x]])
    SysTime = sapply(1:BeatQty, function(x) MaxDomain[x,Sindex[x]])

    PP = P_sys - P_dias
    BeatPeriod = OnsetTimes[-1]-OnsetTimes[-length(OnsetTimes)]

    dyneg = abp[-1]-abp[-length(abp)]
    dyneg[dyneg>0] = 0

    MAP = sapply(1:BeatQty, function(x) mean(abp[OnsetTimes[x]:OnsetTimes[x+1]]))

    RR = BeatPeriod/125
    sys_duration = 0.3*sqrt(RR)
    EndOfSys1 = round(OT + sys_duration*125)

    return(list(time.systole=SysTime,systolic.bp=P_sys,time.diastole=DiasTime,diastolic.bp=P_dias,
        pulse.pressure=PP, mean.bp=MAP, beat.period=BeatPeriod, end.of.systole.rr=EndOfSys1))
}
