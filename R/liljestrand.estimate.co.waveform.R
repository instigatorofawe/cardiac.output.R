#' liljestrand.estimate.co.waveform
#' 
#' Returns a function mapping t to cardiac output.
#' 
#' @param abp Arterial blood pressure waveform
#' @param i0 Average value of cardiac output over interval.
#' @return Cardiac output waveform estimated using Liljestrand estimator.
#' @export
liljestrand.estimate.co.waveform = function(abp, i0) {
    onsets = wabp(abp)
    features = abpfeature(abp,onsets)

    durations = (features$end.of.systole.rr-onsets[-length(onsets)])/125
    hr = 125*60/features$beat.period
    relative.co = features$pulse.pressure/(features$systolic.bp+features$diastolic.bp)*hr
    fractions = relative.co/sum(relative.co)
    total = mean(length(onsets)/hr)*i0

    stroke.volumes = fractions*total
    amplitudes = stroke.volumes / durations

    cardiac.output = rep(0,length(abp))
    for (c in 1:(length(onsets)-1)) {
        cardiac.output[onsets[c]:features$end.of.systole.rr[c]] = amplitudes[c]
    }

    return(cardiac.output)
}