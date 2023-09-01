// Log
capture log close
log using "case-study-log", replace text

// Data
use dataOvarian.dta, clear

// Describe data
codebook
describe
list if patientID <= 5

// stset for progression-free survival
stset timeS, failure(statusS == 1)

// mestreg
mestreg trt || trialID:, dist(exp)
estimates store me_rint_exp
mestreg trt || trialID:, dist(wei)
estimates store me_rint_wei

// Compare
estimates stats *

// Weibull is better
// We now add a random treatment effect too
mestreg trt || trialID: trt, dist(wei) nohr cov(unstr)
estimates store me_rboth_wei
estimates stats me_rint_wei me_rboth_wei
lrtest me_rboth_wei me_rint_wei

// Adding a random treatment effect does not significantly improve fit
// Restore 'better' model and print results
estimates restore me_rint_wei
mestreg
mestreg, nohr
mestreg, notable
mestreg, noheader nohr
mestreg, noheader

//
preserve
replace _t = 2
predict S_cond, surv conditional
predict S_fixed, surv conditional(fixedonly)
predict S_marg, surv marginal
list patientID trialID trt S_* ///
	if patientID == 1 | patientID == 4 | patientID == 129 | patientID == 130
restore

// Hazard predictions
// 'fixedonly':
stcurve, hazard at1(trt=0) at2(trt=1) fixedonly name("h_fixed", replace) outfile(h_fixed, replace) ///
	title("") legend(rows(1) position(12)) xsize(16in) ysize(10in) scale(1.2)
graph export "h-fixed.pdf", replace
// 'marginal':
stcurve, hazard at1(trt=0) at2(trt=1) marginal name("h_marginal", replace) outfile(h_marginal, replace) ///
	title("") legend(rows(1) position(12)) xsize(16in) ysize(10in) scale(1.2)
graph export "h-marginal.pdf", replace

// Survival predictions
// 'fixedonly':
stcurve, survival at1(trt=0) at2(trt=1) fixedonly name("S_fixed", replace) outfile(S_fixed, replace) ///
	title("") legend(rows(1) position(12)) xsize(16in) ysize(10in) scale(1.2)
graph export "S-fixed.pdf", replace
// 'marginal':
stcurve, survival at1(trt=0) at2(trt=1) marginal name("S_marginal", replace) outfile(S_marginal, replace) ///
	title("") legend(rows(1) position(12)) xsize(16in) ysize(10in) scale(1.2)
graph export "S-marginal.pdf", replace

// Predict random effects (and SEs)
predict b_trialID, reffects reses(se_b_trialID)

// Tag distinct random effects
// This is useful for plotting...
bysort trialID: gen plot_tag = _n == 1

// Rank the predicted random effects
egen trialID_rank = rank(b_trialID) if plot_tag == 1

// Make confidence intervals for random effect values
gen b_trialID_lci = b_trialID - invnormal(1 - 0.05 / 2) * se_b_trialID
gen b_trialID_uci = b_trialID + invnormal(1 - 0.05 / 2) * se_b_trialID

// Generate colour for plotting
gen col = 0
replace col = -1 if b_trialID_uci < 0
replace col = 1 if b_trialID_lci > 0

// Plot ranks for BLUPs
twoway ///
	(rcap b_trialID_lci b_trialID_uci trialID_rank if col == 0, color(black)) ///
	(rcap b_trialID_lci b_trialID_uci trialID_rank if col == 1, color(stred)) ///
	(rcap b_trialID_lci b_trialID_uci trialID_rank if col == -1, color(stgreen)) ///
	(scatter b_trialID trialID_rank if col == 0, color(black)) ///
	(scatter b_trialID trialID_rank if col == 1, color(stred)) ///
	(scatter b_trialID trialID_rank if col == -1, color(stgreen)) ///
	, legend(off) yline(0, lpattern(dash) lcolor(gray)) name("b_trialID_ranks", replace) ///
	xtitle("Rank") ytitle("Predicted BLUPs") xsize(16in) ysize(10in) scale(1.2)
graph export "rank-BLUPs.pdf", replace

// Use a more interpretable version for BLUPs
// E.g., 2-years baseline survival
predictnl eta = log(-log(1 - exp(-exp(_b[_t:_cons] + b_trialID) * (2^(exp(_b[/:ln_p])))))), ci(eta_lci eta_uci)
gen S_2Y_trialID = 1 - exp(-exp(eta))
gen S_2Y_trialID_lci = 1 - exp(-exp(eta_lci))
gen S_2Y_trialID_uci = 1 - exp(-exp(eta_uci))

// Reference (b = 0)
local etaref = exp(-exp(_b[_t:_cons]) * (2^(exp(_b[/:ln_p]))))

// Plot it
twoway ///
	(rcap S_2Y_trialID_lci S_2Y_trialID_uci trialID_rank if col == 0, color(black)) ///
	(rcap S_2Y_trialID_lci S_2Y_trialID_uci trialID_rank if col == 1, color(stred)) ///
	(rcap S_2Y_trialID_lci S_2Y_trialID_uci trialID_rank if col == -1, color(stgreen)) ///
	(scatter S_2Y_trialID trialID_rank if col == 0, color(black)) ///
	(scatter S_2Y_trialID trialID_rank if col == 1, color(stred)) ///
	(scatter S_2Y_trialID trialID_rank if col == -1, color(stgreen)) ///
	, legend(off) yline(`etaref', lcolor(gray) lpattern(dash)) name("S_2Y_trialID_ranks", replace) ///
	xtitle("Rank") ytitle("Baseline Survival at 2 Years") xsize(16in) ysize(10in) scale(1.2)
graph export "rank-S-2Y.pdf", replace

// Contextual effects
matrix list e(b)

// Variance:
display _b[/:var(_cons[trialID])]

// Plot:
twoway ///
	(histogram b_trialID) ///
	(function y=normalden(x, sqrt(_b[/:var(_cons[trialID])])), range(-1 1)) ///
	, xtitle(Predicted BLUPs) ytitle(Density) legend(order(1 "Density" 2 "Fitted Normal Distribution") rows(1) position(12)) name("BLUPs_distribution", replace) xsize(16in) ysize(10in) scale(1.2)
graph export "BLUPs-distribution.pdf", replace

// Median hazard ratio:
nlcom exp(sqrt(2 * _b[/:var(_cons[trialID])]) * invnormal(0.75))

// For comparison:
mestreg ib1.trt || trialID:, dist(wei) cov(unstr)
estimates restore me_rint_wei

// IqHR
nlcom exp(sqrt(_b[/:var(_cons[trialID])]) * (invnormal(0.875) - invnormal(0.125)))

// Combine plots
// Hazard:
use h_fixed, clear
gen type = "fixedonly"
append using h_marginal
replace type = "marginal" if type == ""
twoway ///
	(line haz1 _t if type == "fixedonly", sort lcolor(stblue) lpattern(solid)) ///
	(line haz2 _t if type == "fixedonly", sort lcolor(stred) lpattern(solid)) ///
	(line haz1 _t if type == "marginal", sort lcolor(stblue) lpattern(dash)) ///
	(line haz2 _t if type == "marginal", sort lcolor(stred) lpattern(dash)) ///
	, legend(position(12) rows(1) order(1 "trt=0, cond." 2 "trt=1, cond." 3 "trt=0, marg." 4 "trt=1, marg.")) ///
	xtitle("Years (follow-up time)") ytitle("Hazard") xsize(16in) ysize(10in) scale(1.2)
graph export "h-comp.pdf", replace
// Hazard ratios:
gen hr = haz2 / haz1
twoway ///
	(line hr _t if type == "fixedonly", sort lcolor(stblue) lpattern(solid)) ///
	(line hr _t if type == "marginal", sort lcolor(stblue) lpattern(dash)) ///
	, legend(position(12) rows(1) order(1 "conditional" 2 "marginal")) ///
	xtitle("Years (follow-up time)") ytitle("Hazard Ratio") xsize(16in) ysize(10in) scale(1.2)
graph export "hr-comp.pdf", replace

// Survival:
use S_fixed, clear
gen type = "fixedonly"
append using S_marginal
replace type = "marginal" if type == ""
twoway ///
	(line surv1 _t if type == "fixedonly", sort lcolor(stblue) lpattern(solid)) ///
	(line surv2 _t if type == "fixedonly", sort lcolor(stred) lpattern(solid)) ///
	(line surv1 _t if type == "marginal", sort lcolor(stblue) lpattern(dash)) ///
	(line surv2 _t if type == "marginal", sort lcolor(stred) lpattern(dash)) ///
	, legend(position(12) rows(1) order(1 "trt=0, cond." 2 "trt=1, cond." 3 "trt=0, marg." 4 "trt=1, marg.")) ///
	xtitle("Years (follow-up time)") ytitle("Survival") xsize(16in) ysize(10in) scale(1.2)
graph export "S-comp.pdf", replace

// Close log
log close
