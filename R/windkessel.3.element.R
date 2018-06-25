#' windkessel.three.element
#'
#' Generates data according to the 3-element windkessel model
#' @param T end time
#' @param N number of discrete time steps for the interval [0,T]
#' @param r1 Resistance 1
#' @param r2 Resistance 2
#' @param c Capacitance
#' @param p0 Initial value of pressure
#' @param i Function that maps time to cardiac output (mL/s)
#' @return Waveform of ABP over time
#' @export
windkessel.three.element <- function(T, N, r1, r2, c, p0, cardiac.output) {
  times = seq(0, T, by=T/(N-1))

  state_evolution = function(t, state, params) {
    with(as.list(c(state)), {
      dx = cardiac.output(t)/c - x/(r2*c)
      return(list(dx))
    })
  }

  result = deSolve::ode(y = c(x=p0), times = times, func = state_evolution, parms=NULL, rtol = 1e-10, atol = 1e-10)
  cardiac_output = cardiac.output(times)

  result[,2] = result[,2] + cardiac_output * r1

  return(result)

}
