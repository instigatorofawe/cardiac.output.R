#'
#' ESTIMATE_CO_PARLIKAR
#'
#' @param abp Arterial blood pressure waveform (125 Hz, mmHg)
#' @param feat Features computed using abpfeature()
#' @param onsets Beat onsets computed using wabp()
#' @param window.radius Window under which we compute least-squares estimate of the time constant
#' @return Vector of uncalibrated beat-to-beat estimates of cardiac output
#' @export
estimate.co.parlikar = function(abp,feat,onsets,window.radius) {
    deltaP = abp[onsets[-1]] - abp[onsets[-length(onsets)]]
    tau.raw = (2 * (feat[,6]-feat[,4]) - deltaP) / feat[,7]
    tau.indices = lapply(1:dim(feat)[1], function(x) max(x-window.radius,1):min(x+window.radius,dim(feat)[1]))
    tau = sapply(tau.indices, function(x) sum(feat[x,6]*feat[x,6])/sum(feat[x,6]*tau.raw[x]))
    return(deltaP/feat[,7]+feat[,6]/tau)
}
