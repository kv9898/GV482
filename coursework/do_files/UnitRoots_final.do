/**********************/
/*Clean and set memory*/
/**********************/
clear all
set matsize 5000
set more off

/******************/
/*Install programs*/
/******************/
cap ssc install xtabond2 
cap ssc install xtivreg2 
cap ssc install spmat 
cap ssc install spmack

/*********************/
/*Sets base directory*/
/*********************/
global project "/bbkinghome/pascual/winprofile/mydocs/replication_files_ddcg/"  /* Set base directory                         */

use "$project/DDCGdata_final.dta"

/*******************************************************************************************************************************/
/*******************************************************************************************************************************/
/*********************DEFINE REQUIRED PROGRAMS THAT WILL BE USED DURING THE EXECUTION OF THIS DO FILE **************************/
/*******************************************************************************************************************************/
/*******************************************************************************************************************************/
capture program drop llctest
program define llctest, rclass
syntax anything[, lags(integer 3) varalt mumain(real -.53796) sigmamain(real .85408)]
	local 0 `anything' 
	gettoken yvar 0 : 0 /*Variable to be analyzes*/
	gettoken excov 0 : 0, match(par) /*exogenous covariates*/
	gettoken bpvar 0 : 0, match(par) /*variable for balanced panel to estimate long run variance/ short run*/
	
/*Filter data*/
quietly: if `lags'>=1{
reg d.`yvar' l(1/`lags').d.`yvar' yy* i.wbcode2 `excov'
}
quietly: if `lags'==0{
reg d.`yvar' yy* i.wbcode2 `excov'
}
quietly: predict e if e(sample), resid


/*Main estimate of long-run variance/short run innovation variance*/
if `lags'==0{
return scalar SA=1
local SA=1
}
if `lags'==1{
return scalar SA=1/abs(1-(_b[LD.`yvar']))
local SA=1/abs(1-(_b[LD.`yvar']))
}
if `lags'==2{
return scalar SA=1/abs(1-(_b[LD.`yvar']+_b[L2D.`yvar']))
local  SA=1/abs(1-(_b[LD.`yvar']+_b[L2D.`yvar']))
}
if `lags'==3{
return scalar SA=1/abs(1-(_b[LD.`yvar']+_b[L2D.`yvar']+_b[L3D.`yvar']))
local SA=1/abs(1-(_b[LD.`yvar']+_b[L2D.`yvar']+_b[L3D.`yvar']))
}
if `lags'==4{
return scalar SA=1/abs(1-(_b[LD.`yvar']+_b[L2D.`yvar']+_b[L3D.`yvar']+_b[L4D.`yvar']))
local SA=1/abs(1-(_b[LD.`yvar']+_b[L2D.`yvar']+_b[L3D.`yvar']+_b[L4D.`yvar']))
}
if `lags'==5{
return scalar SA=1/abs(1-(_b[LD.`yvar']+_b[L2D.`yvar']+_b[L3D.`yvar']+_b[L4D.`yvar']+_b[L5D.`yvar']))
local  SA=1/abs(1-(_b[LD.`yvar']+_b[L2D.`yvar']+_b[L3D.`yvar']+_b[L4D.`yvar']+_b[L5D.`yvar']))
}
if `lags'==6{
return scalar SA=1/abs(1-(_b[LD.`yvar']+_b[L2D.`yvar']+_b[L3D.`yvar']+_b[L4D.`yvar']+_b[L5D.`yvar']+_b[L6D.`yvar']))
local SA=1/abs(1-(_b[LD.`yvar']+_b[L2D.`yvar']+_b[L3D.`yvar']+_b[L4D.`yvar']+_b[L5D.`yvar']+_b[L6D.`yvar']))
}
if `lags'==7{
return scalar SA=1/abs(1-(_b[LD.`yvar']+_b[L2D.`yvar']+_b[L3D.`yvar']+_b[L4D.`yvar']+_b[L5D.`yvar']+_b[L6D.`yvar']+_b[L7D.`yvar']))
local SA=1/abs(1-(_b[LD.`yvar']+_b[L2D.`yvar']+_b[L3D.`yvar']+_b[L4D.`yvar']+_b[L5D.`yvar']+_b[L6D.`yvar']+_b[L7D.`yvar']))
}


quietly: if `lags'>=1{
reg l.`yvar' l(1/`lags').d.`yvar' yy* i.wbcode2 `excov'
}
quietly: if `lags'==0{
reg l.`yvar'  yy* i.wbcode2 `excov'
}
quietly: predict v if e(sample), resid

/*Estimate country variances to control for country-level heteroskedasticity*/
quietly: if `lags'>=1{
reg d.`yvar' l.`yvar' l(1/`lags').d.`yvar' yy* i.wbcode2
}
quietly: if `lags'==0{
reg d.`yvar' l.`yvar' yy* i.wbcode2
}
quietly: predict residual if e(sample), resid
quietly: gen res2=residual^2
quietly: bysort wbcode2: egen vari=mean(res2)
quietly: gen sigmai=sqrt(vari)

/*Alternative routine for long run variance/short run obtained from llc original test*/
quietly: if "`varalt'"!=""{
preserve
keep if `bpvar'==1 /*only operates in balanced panel*/
xtunitroot llc `yvar', demean lags(`lags')
return scalar SA=r(sbar)
local SA=r(sbar)
restore
}

/*Compute panel test statistic*/
quietly: gen e_est=e/sigmai
quietly: gen v_est=v/sigmai

quietly: reg e_est v_est, noconstant
return scalar tdelta=_b[v_est]/_se[v_est]
return scalar stddelta=_se[v_est]
local tdelta=_b[v_est]/_se[v_est]
local stddelta=_se[v_est]

quietly: predict z if e(sample), resid
quietly: gen z2=z^2
quietly: sum z2
return scalar sigma2=r(mean)
return scalar num=r(N)
local sigma2=r(mean)
local num=r(N)

/*Obtain adjustment mu and sigma*/

local tadjmain=(`tdelta'-`num'*`SA'*`stddelta'*`mumain'/`sigma2')/`sigmamain'
return scalar pvalmain=normal(`tadjmain')
return scalar  tadjmain=(`tdelta'-`num'*`SA'*`stddelta'*`mumain'/`sigma2')/`sigmamain'



drop e v residual res2 vari sigmai e_est v_est z z2

end



/****************************/
/****************************/
/*Table 2: Row for LLC tests*/
/****************************/
/****************************/
llctest y (dem) (), lags(0)
display r(tadjmain)
display r(pvalmain)

llctest y (dem) (), lags(1)
display r(tadjmain)
display r(pvalmain)

llctest y (dem) (), lags(3)
display r(tadjmain)
display r(pvalmain)

llctest y (dem) (), lags(7)
display r(tadjmain)
display r(pvalmain)

/****************************/
/****************************/
/*Table 4: Row for LLC tests*/
/****************************/
/****************************/
llctest y (dem interfull*) (), lags(3)
display r(tadjmain)
display r(pvalmain)

llctest y (dem sov1 sov2 sov3 sov4) (), lags(3)
display r(tadjmain)
display r(pvalmain)

llctest y (dem l(1/4).unrest) (), lags(3)
display r(tadjmain)
display r(pvalmain)

llctest y (dem l(1/4).tradewb) (), lags(3)
display r(tadjmain)
display r(pvalmain)

llctest y (dem regdum*) (), lags(3)
display r(tadjmain)
display r(pvalmain)

llctest y (dem l(1/4).nfagdp) (), lags(3)
display r(tadjmain)
display r(pvalmain)
