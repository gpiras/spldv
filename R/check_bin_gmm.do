*** Check GMM probit model in Stata ***

use http://www.stata-press.com/data/r13/auto

global xb "{b1}*gear_ratio + {b2}*length + {b3}*headroom + {b0}"
global Phi "normal($xb)"

gmm (foreign - $Phi), instruments(gear_ratio length headroom) onestep winitial(unadjusted) vce(unadjusted)
estimates store gmm1
gmm (foreign - $Phi), instruments(gear_ratio length headroom) onestep winitial(unadjusted) vce(robust)
estimates store gmm2
gmm (foreign - $Phi), instruments(gear_ratio length headroom) twostep winitial(unadjusted) vce(unadjusted) wmatrix(robust)
estimates store gmm3
gmm (foreign - $Phi), instruments(gear_ratio length headroom) twostep winitial(unadjusted) vce(robust)     wmatrix(robust)
estimates store gmm4
gmm (foreign - $Phi), instruments(gear_ratio length headroom) twostep winitial(unadjusted) vce(unadjusted) wmatrix(unadjusted)
estimates store gmm5
gmm (foreign - $Phi), instruments(gear_ratio length headroom) twostep winitial(unadjusted) vce(robust)     wmatrix(unadjusted)
estimates store gmm6
estimates table gmm1 gmm2 gmm3 gmm4 gmm5 gmm6, b se


