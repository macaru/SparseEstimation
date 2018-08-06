/*************************************************************************
The macro "ose_scad" returns one-step estimator that was suggested by Zou and Li (2008)
The macro "lasso_shot" returns LASSO estimator using shooting algorithm that was suggested by Fu (1998)
The macros "diff" and "cv"returns the value of K folds cross validation for estimated coefficients
The macro "optlam" searchs the optimal smoothing parameter based on K folds cross validation

Example for usage of macro.
The name of dataset: dat
The name of response variable: y
The name of predict variables : col1 - col12
The min of smoothing parameter for grid search: 0.1
The max of smoothing parameter for grid search: 1
The length of grid for smoothing parameter lambda: 30 
%optlam(dat, y, col1 - col12, 0.1, 1, 30);

You could get 3 datasets optcoef, optlambda and seqcoef
The dataset optcoef is estimated coefficient
The dataset seqcoef shows coefficients when lambda tends to large
The dataset optlambda shows the value of K folds CV when lambda tends to large

[References]
[1] Fu (1998). Penalized regression: The bridge versus the lasso. 
   Journal of Computational and Graphical Statistics 7: 397 -416.
[2] Zou and Li (2008). One-step sparse estimates in nonconacve penalized likelihood models. 
   The Annals of Statistics 36: 267-288.

Last changed: 06 AUG 2018 by Masaru Kanba
*************************************************************************/

%macro ose_scad(dat, res, xvar, lambda, out);

data _null_;
  set &dat. end = end;
  if (end = 1) then call symputx("n", _n_);
run;

data _null_;
  nl =  2 * &n. * &lambda.;
  call symputx("nl", nl);
run;

/***** MLE *****/
ods output ParameterEstimates = _parm OverallANOVA = var PredictedValues = _yhat;
proc glm data = &dat.;
  model &res. = &xvar. / noint predicted;
run;
quit;

data _null_;
  set var;
  if (source = "Error") then call symputx("sig", ms);
run;

/*****************
    Step 1a. 
******************/

/* design matrix */
data X;
  set &dat.;
  keep &xvar.;
run;

proc iml;
  use X;
    read all into _X[colname = varname];
  _xs = sqrt(1/ (&sig. ) ) * _X;
  create _x_s from _xs[colname = varname];
  append from _xs;
quit;

data y_s;
  set _yhat;
  yhat = sqrt(1/(&sig.) ) * Predicted;
  keep yhat;
run;

/*****************
    Step 1b. 
******************/

/* p'(|beta|^0) */
data _p_dev;
  length u v $1000.;
  set _parm end = end;
  retain u v cu cv;
  counter + 1;
  lam = &lambda.;
  a = 3.7;
  alam = a * lam;
  abeta = abs(estimate);
  if (_n_ = 1) then do;
    u = ' ';
	v = ' ';
	cu = 0;
	cv = 0;
  end;
  if (abeta <= lam) then p = lam;
  if (lam < abeta < alam) then p = (alam - abeta) / (a - 1);
  if (alam <= abeta) then p = 0;

  if (p = 0) then u = trim(u) || ' ' || compress(parameter);
  if (p > 0) then v = trim(v) || ' ' || compress(parameter);
  if (p > 0) then lam_p = lam / p;
  if (p = 0) then cu = cu + 1;
  if (p > 0) then cv = cv + 1;

  if (end = 1) then do;
    call symputx('xu', u);
	call symputx('xv', v);
    call symputx('cu', cu);
	call symputx('cv', cv);
  end;
run;

data label_bu;
  set _p_dev;
  if (p = 0) then output;
  keep parameter counter;
run;

data label_bv;
  set _p_dev;
  if (p > 0) then output;
  keep parameter counter;
run;

/* lambda / p'(|beta|^0) */
data _p0;
  set _p_dev;
  if (p > 0) then output;
  keep lam_p;
run;

%put U: &xu. # V: &xv. ;

%if &cu. ^= 0 %then %do;
data xu_s;
  set _x_s;
  keep &xu.;
run;
%end;
%if &cv. ^= 0 %then %do;
data xv_s;
  set _x_s;
  keep &xv.;
run; 

proc iml;
  use xv_s;
    read all into _x[colname = varname];
    read all into _x2[colname = varname];
  use _p0;
    read all into _p[colname = varname];
  do i = 1 to &n.;
    do j = 1 to &cv.;
      _x2[i,j] = _x[i,j] * _p[j,1];
	end;
  end;
  create xv_s2 from _x2;
  append from _x2;
quit;
%end;


/*****************
    Step 1c. 
******************/

proc iml;
  %if &cu. ^= 0 and &cv. ^= 0 %then %do;
  use xv_s2;
    read all into xv_s2[colname = varname];
  use xu_s;
    read all into xu_s[colname = varname];
  use y_s;
    read all into y_s[colname = varname];
  Hu = xu_s * inv(t(xu_s) * xu_s) * t(xu_s);
  _Xv_ss = xv_s2 - Hu * xv_s2;
  _y_ss = y_s - Hu * y_s;
  create xv_ss from _xv_ss;
  append from _xv_ss;
  create y_ss from _y_ss;
  append from _y_ss;
  %end;
  %if &cu. = 0 %then %do;
  use xv_s2;
    read all into xv_s2[colname = varname];
  use y_s;
    read all into y_s[colname = varname];
  _Xv_ss = xv_s2;
  _y_ss = y_s;
  create xv_ss from _xv_ss;
  append from _xv_ss;
  create y_ss from _y_ss;
  append from _y_ss;
  %end;
quit;

/*****************
    Step 2. 
******************/

%if &cv. ^= 0 %then %do;
data _dat_lar;
  merge y_ss (rename = (COL1 = y)) xv_ss;
run;

%lasso_shot(y, col1 - col&cv., _dat_lar, &nl., _parm_lars);

data bv_s;
  set _parm_lars;
  keep estimate;
run;

data bv_s02;
  merge bv_s _p0;
  bv = estimate * lam_p;
run;

data bv_s03;
  set bv_s02;
  keep bv;
run;

data bv;
  merge label_bv bv_s02;
  beta = bv;
run;

%end;

/*****************
    Step 3. 
******************/

proc iml;
  %if &cu. ^= 0 and &cv. ^= 0 %then %do;
  use xv_s;
    read all into xv_s[colname = varname];
  use xu_s;
    read all into xu_s[colname = varname];
  use bv_s03;
    read all into bv_s[colname = varname];
  use y_s;
    read all into y_s[colname = varname];
  _bu = inv(t(xu_s) * xu_s) * t(xu_s) * (y_s - xv_s * bv_s);
  create bu_s from _bu;
  append from _bu;
  %end;
  %if &cv. = 0 %then %do;
  use xu_s;
    read all into xu_s[colname = varname];
  use y_s;
    read all into y_s[colname = varname];
  _bu = inv(t(xu_s) * xu_s) * t(xu_s) * y_s;
  create bu_s from _bu;
  append from _bu;
  %end;
quit;

%if &cu. ^= 0 %then %do;
data bu;
  merge label_bu bu_s;
  beta = col1;
run;
%end;

%if &cu. ^= 0 and &cv. ^= 0 %then %do;
data _beta;
  set bu bv;
  keep parameter counter beta;
run;
%end;

%if &cu. = 0 %then %do;
data _beta;
  set bv;
  keep parameter counter beta;
run;
%end;

%if &cv. = 0 %then %do;
data _beta;
  set bu;
  keep parameter counter beta;
run;
%end;

proc sort data = _beta out = _coef;
  by counter;
run;

data &out.;
  set _coef;
run;

proc datasets lib = work;
  delete _: x xv_s xv_s2 bu bu_s bv bv_s bv_s02 bv_s03 y_s label: hu: var xu_s xv_ss y_ss;
run;
quit;

%mend;

%macro lasso_shot(y, xvar, dat, lambda_l, coef_la);

data _null_;
  set &dat. end = end;
  array x{*} &xvar.;
  if (end = 1) then call symputx("m", dim(x) );
run;

ods output ParameterEstimates = _parm_ols;
proc glm data = &dat.;
  model &y. = &xvar. / noint predicted;
run;
quit;

proc transpose data = _parm_ols out = _parm_ols2;
  id parameter;
  var estimate;
run;

data old_beta;
  set _parm_ols2;
  array b{*} oldb1 - oldb&m.;
  array bhat{*} &xvar.;
  do j3 = 1 to &m.;
    b{j3} = bhat{j3};
  end;
  dummy = 1;
  keep dummy oldb1 - oldb&m.;
run;

data _temp01;
  set &dat.;
  dummy = 1;
run;

%let flg = 1;
%let cflg = 10;
%let difflg = 1;

%do %while( (&cflg. > 0.1) and (&difflg. = 1) );

%if (&flg. ^= 1) %then %do;
data old_beta;
  set new_beta;
  array b{*} oldb1 - oldb&m.;
  array bhat{*} bhat1 - bhat&m.;
  do j = 1 to &m.;
    b{j} = bhat{j};
  end;
  dummy = 1;
  keep dummy oldb1 - oldb&m.;
run;
%end;

data _temp02;
  merge _temp01 old_beta;
  by dummy;
run;

data _temp03;
  set _temp02;
  array x{*} &xvar.;
  array s{*} s1 - s&m.;
  array t{*} t1 - t&m.;
  array u{*} u1 - u&m.;
  array oldb{*} oldb1 - oldb&m.;

  do j = 1 to dim(x);
    s{j} = 2 * (x{j} ** 2);
	u{j} = -2 * x{j} * y;
	t0 = 0;
	do i = 1 to dim(x);
	  if (j ne  i) then t0 = t0 +  2 * x{j} * x{i} * oldb{i};
	end;
	t{j} = t0;
  end;
  keep &xvar. s1 - s&m. t1 - t&m. u1 - u&m.;
run;

proc means data = _temp03;
  var s1 - s&m. t1 - t&m. u1 - u&m.;
  output out = _temp04 (drop = _:) sum = s1 - s&m. t1 - t&m. u1 - u&m.;
run;

data new_beta;
  merge _temp04 old_beta;
  array _shat{*} _shat1 - _shat&m.;
  array oldb{*} oldb1 - oldb&m.;
  array bhat{*} bhat1 - bhat&m.;
  array s{*} s1 - s&m.;
  array t{*} t1 - t&m.;
  array u{*} u1 - u&m.;
  do j = 1 to &m.;
    _shat{j} =  t{j} + u{j};
	if (_shat{j} >= &lambda_l.) then bhat{j} = (&lambda_l. - _shat{j}) / s{j};
	if (_shat{j} <= -1 * &lambda_l.) then bhat{j} = (-1 * &lambda_l. - _shat{j}) / s{j};
	if (abs(_shat{j}) < &lambda_l.) then bhat{j} = 0;
	if (j = 1) then cflg = 0;
	cflg = cflg + abs(bhat{j} - oldb{j});
  end;
  difflg = 1;
  if ( (&cflg. - cflg) < 0) then difflg = 2;
  call symputx("cflg", cflg);
  call symputx("difflg", difflg);
run;

%let flg = 2;

%end;

proc transpose data = new_beta out = new_beta02;
  var bhat1 - bhat&m.;
run;

data &coef_la.;
  merge _parm_ols new_beta02;
  estimate = col1;
  keep parameter estimate;
run;

proc datasets lib = work;
  delete _temp: new: old_beta;
run;
quit;

%mend;

%macro diff(dat, res, xvar, beta, outd);

data _diff_y;
  set &dat.;
  keep &res.;
run;

data _diff_x;
  set &dat.;
  keep &xvar.;
run;

data beta;
  set &beta.;
  keep beta;
run;

proc iml;
  use _diff_y;
    read all into _diff_y[colname = varname];
  use _diff_x;
    read all into _diff_x[colname = varname];
  use beta;
    read all into beta[colname = varname];
  _cv = _diff_y - _diff_x * beta;
  __cv = t(_cv) * _cv;
  create ___cv from __cv;
  append from __cv;
quit;

data &outd.;
  set ___cv;
  rename col1 = cv;
run;

proc datasets lib = work;
  delete _cv __cv ___cv _diff: beta;
run;

%mend;


%macro cv(dat, res, xvar, K, lam, outcv);

data cvdata;
  set &dat.;
  s = rand('UNIFORM');
run;

proc sort data = cvdata; by s; run;


data _null_;
  set &dat. end = end;
  if (end = 1) then call symputx("_n", _n_);
run;

%do i = 1 %to &K.;

%if (&i. < &K.) %then %do;
data train&i. test&i.;
  set cvdata;
  low = (&i. - 1) * ceil(&_n. / &K. ) + 1;
  high = (&i. ) * ceil(&_n. / &K. );
  if (low <= _n_ <= high) then output test&i.;
  if (low > _n_  or _n_ > high) then output train&i.;
run;

%ose_scad(train&i., &res., %str(&xvar.), &lam., coef&i.);
%diff(test&i., &res., %str(&xvar.), coef&i., _cv_temp);

%end;

%if (&i. = &K.) %then %do;
data train&i. test&i.;
  set cvdata;
  low = (&i. - 1) * ceil(&_n. / &K. ) + 1;
  if (low <= _n_ <= &_n.) then output test&i.;
  if (low > _n_  or _n_ > &_n.) then output train&i.;
run;

%ose_scad(train&i., &res., %str(&xvar.), &lam., coef&i.);
%diff(test&i., &res., %str(&xvar.), coef&i., _cv_temp);

%end;

%if &i. = 1 %then %do;
data fcv;
  set _cv_temp;
run;
%end;

%if &i. > 1 %then %do;
data fcv;
  set fcv _cv_temp;
run;
%end;

%end;

proc means data = fcv;
  var cv;
  output out = &outcv. (drop = _:) sum = cv;
run;

proc datasets lib = work;
  delete cvdata _cv: test: train: fcv coef:;
run;

%mend;

%macro optlam(dat, res, xvar, minl, maxl, len);

data lambda;
  band = (&maxl. - &minl. ) / (&len. - 1);
  do i = 1 to &len.;
    lam = &minl + band * (i - 1);
	output;
  end;
run;

%do l = 1 %to &len.;

data _null_;
  set lambda;
  if (_n_ = &l.) then call symputx('lam', lam);
run;

%put Lambda &l. : &lam.;


%cv(&dat., &res., &xvar., 5, &lam., cv_t);
%ose_scad(&dat., &res., %str(&xvar.), &lam., coef);

proc transpose data = coef out = coef_t;
  id counter;
  idlabel Parameter;
run;

data coef_t;
  set coef_t;
  num = &l.;
  lam = &lam.;
run;

data cv_t;
  set cv_t;
  num = &l.;
run;

%if &l. = 1 %then %do;
data cv;
  set cv_t;
run;

data seqcoef;
  set coef_t;
run;

%end;

%if &l. > 1 %then %do;
data cv;
  set cv cv_t;
run;

data seqcoef;
  set seqcoef coef_t;
run;

%end;

%end;

data lambda02;
  merge lambda cv;
run;

data optlambda;
  set lambda02 end = end;
  retain mincv minlam flg;
  if (_n_ = 1) then flg = 1;
  if (flg = 1 and cv ^= .) then do;
    mincv = cv;
	minlam = lam;
	flg = 2;
  end;
  if (flg = 2 and mincv > cv > .) then do;
    mincv = cv;
    minlam = lam;
  end;
  if (end = 1) then call symputx("optlam", minlam);
run;

%ose_scad(&dat., &res., %str(&xvar.), &optlam., optcoef);

proc datasets lib = work;
  delete lambda: cv cv_t coef_t coef;
run;
quit;

%mend;
