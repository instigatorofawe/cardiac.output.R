#' create.co.mapping
#' 
#' Given a 125 Hz CO waveform, returns a function mapping time to cardiac output.
#' @param co Cardiac output waveform.
#' @return Function mapping t to cardiac output
#' @export
create.co.mapping = function(co) {
    cardiac.output = function(t) {
        times = seq(from=0,by=1/125,to=(length(co)-1)/125)
        indices = sapply(t, function(x) which.min(abs(x-times)))
        return(co[indices])
    }

    return(cardiac.output)
}