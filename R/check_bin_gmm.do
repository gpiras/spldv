*** Check GMM probit model in Stata ***

use http://www.stata-press.com/data/r13/auto

global xb "{b1}*gear_ratio + {b2}*length + {b3}*headroom + {b0}"
global Phi "normal($xb)"

gmm (foreign - $Phi), instruments(gear_ratio length headroom) onestep winitial(unadjusted) vce(unadjusted)
estimates store gmm1
gmm (foreign - $Phi), instruments(gear_ratio length headroom) onestep winitial(unadjusted) 
estimates store gmm2
estimates table gmm1 gmm2, b se

