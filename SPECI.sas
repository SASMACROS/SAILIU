
/************************************************************************************************************/
/*       NAME: SPECI.SAS                                         				  	    */
/*      TITLE: SAS Macro for Covariate Specification in Linear, Logistic, and Survival Regression           */
/*     AUTHOR: Sai Liu, MPH, Stanford University                                          	            */
/*         OS: Windows 7 Ultimate 64-bit								    */
/*   Software: SAS 9.4											    */
/*       DATE: 29 DEC 2016                                        					    */
/*DESCRIPTION: This program produces a one-page report to help users compare and select the appropriate     */
/*			   functional form of a variable in continuous, categorical, and spline form        */
/*													    */
/*    Copyright (C) <2017>  <Sai Liu and Margaret R Stedman>				                    */
/*													    */
/*    This program is free software: you can redistribute it and/or modify				    */
/*    it under the terms of the GNU General Public License as published by			     	    */
/*    the Free Software Foundation, either version 3 of the License, or					    */
/*    (at your option) any later version.								    */
/*													    */
/*    This program is distributed in the hope that it will be useful,				            */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of					    */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the					    */
/*    GNU General Public License for more details.							    */
/*													    */
/*    You should have received a copy of the GNU General Public License					    */
/*    along with this program.  If not, see <http://www.gnu.org/licenses/>.				    */
/************************************************************************************************************/



options nofmterr symbolgen mlogic notes;

%MACRO SPECI;
%if &model=1 %then %do;
%DATAPREP;
%SPECI_LINEAR ;
%REPORT;
%end;
%if &model=2 %then %do;
%DATAPREP;
%SPECI_LOGISTIC ;
%REPORT;
%end;
%if &model=3 %then %do;
%DATAPREP;
%SPECI_SURVIVAL ;
%REPORT;
%end;
%MEND SPECI;


* Spline macro Reference: http://biostat.mc.vanderbilt.edu/wiki/pub/Main/SasMacros/survrisk.txt from Frank Harrell;
********************************************************************************************************************;
%MACRO RCSPLINE(x,knot1,knot2,knot3,knot4,knot5,knot6,knot7,knot8,knot9,knot10,
  norm=);

  %LOCAL j v7 k tk tk1 t k1 k2;
  %LET v7=&x; %IF %LENGTH(&v7)=8 %THEN %LET v7=%SUBSTR(&v7,1,7);
      %*Get no. knots, last knot, next to last knot;
        %DO k=1 %TO 10;
        %IF %QUOTE(&&knot&k)=  %THEN %GOTO nomorek;
        %END;
  %LET k=11;
  %nomorek: %LET k=%EVAL(&k-1); %LET k1=%EVAL(&k-1); %LET k2=%EVAL(&k-2);
  %IF &k<3 %THEN %PUT ERROR: <3 KNOTS GIVEN.  NO SPLINE VARIABLES CREATED.;
          %ELSE %DO;
          %LET tk=&&knot&k;
          %LET tk1=&&knot&k1;
          DROP _kd_; _kd_=
          %IF &norm=0 %THEN 1;
          %ELSE %IF &norm=1 %THEN &tk - &tk1;
          %ELSE (&tk - &knot1)**.666666666666; ;
                  %DO j=1 %TO &k2;
                   %LET t=&&knot&j;
  &v7&j=max((&x-&t)/_kd_,0)**3+((&tk1-&t)*max((&x-&tk)/_kd_,0)**3
                     -(&tk-&t)*max((&x-&tk1)/_kd_,0)**3)/(&tk-&tk1)%STR(;);
                  %END;
          %END;
%MEND RCSPLINE;


***********************************************************************************************
*** Part I Data set Preparation ***;
***********************************************************************************************;

%MACRO DATAPREP;

ods html close;
ods listing;

/*******Report error message in log, if none of below "MUST INPUT" variable is specified*********/
%if &dataout.= or &dataset.= or &model.= or &xvar_cont.= or &xvar_cat.= or &num_cat.= %then %do; 
    %put ERROR: any of (&dataout.,&dataset.,&model.,&xvar_cont.,&xvar_cat.,&num_cat.)variable is not specified.; %return;
%end;


* Remove "" from dataout path;
%let dataout=%qsysfunc(dequote(&dataout));

*Read in dataset from the folder you assign, or read in dataset from work library, if datain leaves blank;
%if &datain.^=" "  %then %do;
	libname lib &datain.;
	data mydata;set lib.&dataset.(where=(&xvar_cont.^=.));run;
%end;

%else %do;
	data mydata;set work.&dataset.(where=(&xvar_cont.^=.));run;
%end;


/*******Default option setting up*********/
  * knots=4, if leave knot= blank;
%if &knot.=  %then %do;
   %let knot=4;				
%end;
%else %do;
   %let knot=&knot.;
%end;

* norm=2, if leave norm= blank;
%if &norm.=  %then %do;
   %let norm=2;				
%end;
%else %do;
   %let norm=&norm.;
%end;

* ref_xvar_cat=1, if leave ref_xvar_cat= blank;
%if &ref_xvar_cat.=  %then %do;
   %let ref_xvar_cat=1;				
%end;
%else %do;
   %let ref_xvar_cat=&ref_xvar_cat.;
%end;

* reportname=Model Diagnostic Report, if leave reportname= blank;
%if &reportname.=  %then %do;
   %let reportname=Model Diagnostic Report;				
%end;
%else %do;
   %let reportname=&reportname.;
%end;


* get the percentile value for each knot;
proc univariate data=mydata noprint;
	var &xvar_cont.;
	output out=pcts pctlpts= %substr(&knot1,2) %substr(&knot2,2) %substr(&knot3,2) %substr(&knot4,2) %substr(&knot5,2)
						     %substr(&knot6,2) %substr(&knot7,2) %substr(&knot8,2) %substr(&knot9,2) %substr(&knot10,2)
	pctlpre=p;run;

* Call the %RCSPLINE to generate spline terms and save it in "Data_prep" dataset,given the percentile value from above step and method of normalization you assigned;
data Data_temp;
   if _N_ = 1 then
		set pcts;
		set mydata;
		%rcspline(&xvar_cont.,&knot1., &knot2., &knot3., &knot4.,&knot5.,
							  &knot6., &knot7., &knot8., &knot9.,&knot10., norm=&norm.);   														
temp=1;run;

* uniform order of categories to be newcat1~newcatN;
proc freq data=mydata;
	tables &xvar_cat./out=freq;
run;

data sub;
	set freq;
	keep &xvar_cat.;
run;
proc sort data=sub;by &xvar_cat;run;

data neworder;
   set sub;
   by &xvar_cat;
/*   if &xvar_cat=%eval(&ref_xvar_cat) then do;*/
/*       %let newref=newcat;*/
/*   end;*/
   newcat+1;   
run;

proc sql;
   create table addorder as
   select a.* , b.newcat
   from Data_temp a, neworder b
   where a.&xvar_cat = b.&xvar_cat;
quit;


* switch the &ref_xvar_cat. to be &newref_xvar_cat.;
data ref;
	set addorder;
	oldref=%eval(&ref_xvar_cat);
	if &xvar_cat=oldref then newref=newcat;
	temp=1;
run;
proc sort data=ref(where =(newref^=.) keep=temp oldref newref) nodupkey;by newref;run;

* generate dummy variables &xvar_cat.dum1-&xvar_cat.dum%eval(&num_cat) based on "newcat";
data data_prep(drop=d);
  merge addorder ref;
  by temp;

  array &xvar_cat.dum(%eval(&num_cat)) &xvar_cat.dum1-&xvar_cat.dum%eval(&num_cat);
  do d = 1 to %eval(&num_cat);
   &xvar_cat.dum(d)=(newcat=d);
  end;
run;	

%mend DATAPREP;

******************************************************************************************************;
****** Part II - Call specific macro for the model ***;
******************************************************************************************************;

******************************************************************************************************;
*******Linear Regression Macro ****;
******************************************************************************************************;

%MACRO SPECI_LINEAR;
* Run separate Linear regression models with exposure variable as continous, categorical and spline;
proc reg data=data_prep outest=est1Line;
  model &yvar.= &xvar_cont. &covarlist_cat. &covarlist_cont./selection=rsquare ADJRSQ rsquare mse aic bic;run;quit;
proc reg data=data_prep outest=est1cat;
  model &yvar.= &xvar_cat.dum1-&xvar_cat.dum%eval(&num_cat) &covarlist_cat. &covarlist_cont./selection=rsquare ADJRSQ rsquare mse aic bic; run; quit;	
proc reg data=data_prep outest=est1Spline;
   model &yvar. = &xvar_cont. &xvar_cont.1 -- &xvar_cont.%eval(&knot-2) &covarlist_cat. &covarlist_cont./selection=rsquare ADJRSQ rsquare mse aic bic;run; quit;

** Get estimates from model;
* continuous;
data Line(keep=int_line linpred temp);
   set est1line(rename=(intercept=int_line)) nobs=nobs;
   LinPred = &xvar_cont.;temp = 1;if _n_=nobs;run;
* categorical;
data est1cat1;set est1cat;temp=1;run;
data est1cat;
	merge est1cat1 ref;
	by temp;
  CALL SYMPUT('newref_xvar_cat', PUT(newref, 3.));
run;
PROC SQL;
	CREATE TABLE Cat_temp AS
	SELECT *, intercept as int_cat, max(_IN_) as keep
	from est1cat
	where &xvar_cat.dum%eval(&newref_xvar_cat)=.
	having _in_=calculated keep;
quit;
* only keep estimates from the model with reference group = 0;
data Cat(keep=int_cat &xvar_cat.dum1-&xvar_cat.dum%eval(&num_cat) temp); 
   set Cat_temp;
   temp = 1; if &xvar_cat.dum%eval(&newref_xvar_cat)=. then &xvar_cat.dum%eval(&newref_xvar_cat)=0;run;
* Spline;
data Spline(keep=int_sp sp sp1--sp%eval(&knot-2) temp);																		
   set est1Spline (rename=(intercept=int_sp &xvar_cont.=sp)) nobs=nobs;
   array a{%eval(&knot-2)} sp1-sp%eval(&knot-2);
   array b(%eval(&knot-2)) &xvar_cont.1-&xvar_cont.%eval(&knot-2);
   do i =1 to %eval(&knot-2);
  	a(i)=b(i);
	end; temp = 1; if _n_=nobs;run;			

data Data_est;set Line; set Cat;set Spline;run;

data Data_plot;
   merge Data_prep (drop= &xvar_cat.dum1-&xvar_cat.dum%eval(&num_cat)) Data_est;
   by temp;run;

** Calculate Predicted values for Y variable for continuous, categorical and spline models;
%macro cat;
data data_cat;
	set Data_plot;
	%do j=1 %to &num_cat.;
	* model with exposure variable as categorical;
	if newcat=&j. then cat=int_cat+&xvar_cat.dum&j.; 		/* Y-hat = a + b(&xvar_cat.&j.), if &xvar_cat.=&j. */
	%end;run;
%mend cat;
%cat;

data Data_final(drop=k sum);
   set data_cat;
   * model with exposure variable as continuous ;
   Linear = (int_line+(&xvar_cont.*LinPred));				/* Y-hat = a + b(&xvar_cont.)*x(&xvar_cont.) */
   * model with exposure variable as spline ;
   sum=0;
   array c{%eval(&knot-2)} sp1-sp%eval(&knot-2);
   array d(%eval(&knot-2)) &xvar_cont.1-&xvar_cont.%eval(&knot-2);
   do k =1 to %eval(&knot-2);
   sum=sum+c(k)*d(k);
   end;
   Spline = (int_sp+(&xvar_cont.*sp) + sum);				/* Y-hat = a + b(&xvar_cont.)*x(&xvar_cont.) + b(sp1)*x(&xvar_cont.) + b(sp2)*x(&xvar_cont.) + .... b(spn)*x(&xvar_cont.)*/			
retain;run;		

* Centralize predictions at mid-point to align predicted plots;
* Calculate the mid-point value of exposure variable of interest by ((max - min)/2 + min);
PROC SQL;
	CREATE TABLE mid AS
	SELECT temp, Linear as mid_linear, cat as mid_cat,spline as mid_spline, min(&xvar_cont.)as min, max(&xvar_cont.) as max, 
		  ((calculated max)-(calculated min))/2 + (calculated min) as mid, abs((&xvar_cont. - (calculated mid))) as abs_diff
	FROM Data_final	order by abs_diff;QUIT;
* Keep predicted values for Y from mid-point;
proc sort data=mid(firstobs=1 obs=1) out=mid1;by temp;run;
*  calculate the centralized value by subtracting midpoint Y from each predicted value Y;
data data_fig;
	merge Data_final mid1;
	by temp;
	ctr_line=Linear-mid_linear;
	ctr_cat=Cat-mid_cat;
	ctr_spline=Spline-mid_spline;
	* Assign min and max value as new macro variables so as to set up x axis range in the plot;
	call symput('minvalue',min);
	call symput('maxvalue',max);
	label ctr_line='Continuous';
	label ctr_cat='Categorical';
	label ctr_spline='Spline';

run;
* Sort dataset by exposure variable;
proc sort data=data_fig out=Data_fig;by &xvar_cont.;run;

/*Plot Figure A (Centralized Predicted Y with X) */;
ods listing gpath="&dataout.";
ods  graphics on / reset=index imagefmt=png imagename="FigA - &xvar_cont." ;
proc sgplot data =Data_fig;											
series x=&xvar_cont. y=ctr_line/LINEATTRS = (color=ROSE THICKNESS = 4);
series x=&xvar_cont. y=ctr_cat/LINEATTRS = (color=o  THICKNESS = 4);
series x=&xvar_cont. y=ctr_spline/LINEATTRS = (color=BIGB THICKNESS = 4);
YAXIS LABEL ="Centralized Predicted &yvar." ; 
XAXIS values=(&minvalue. to &maxvalue.) LABEL="&xvar_cont." ; 			
run;

ods graphics off;
ods printer close;

/*Figure B (Summary Table for Comparing Model Statistics) */;
* Integrate Statistics from models with exposure variable as continuous, categories, and spline;
data table_line;length model $ 10;set est1line (keep=_ADJRSQ_ _EDF_ _MSE_ _RSQ_ _AIC_ _BIC_) nobs=nobs;model="Continuous"; if _n_=nobs;run;
data table_cat;length model $ 10;set est1cat (keep=_ADJRSQ_ _EDF_ _MSE_ _RSQ_ _AIC_ _BIC_) nobs=nobs;model="Categorical"; if _n_=nobs;run;
data table_Spline;length model $ 10;set est1Spline (keep=_ADJRSQ_ _EDF_ _MSE_ _RSQ_ _AIC_ _BIC_) nobs=nobs;model="Spline"; if _n_=nobs;run;
data table1;set table_line table_cat table_Spline;run;

* Report Summary Table;
options printerpath=png nodate papersize=('8in','5in');
ods _all_ close;
ods printer file="&dataout.\FigB - &xvar_cont..png";
Proc report data=table1 nowd ;
column model _EDF_ _RSQ_ _ADJRSQ_ _MSE_ _AIC_ _BIC_;
define model /"Diagnostic Statistics" group order=data;
define _EDF_/ "Error Degree of Freedom (larger is better)" analysis format=10.5 ;
define _RSQ_/ "R-Squared (larger is better)" analysis format=10.5 ;
define _ADJRSQ_/ "Adjusted R-Squared (larger is better)" analysis format=10.5 ;
define _MSE_/ "Mean Squared Error (smaller is better)" analysis format=10.5 ;
define _AIC_/ "AIC (smaller is better)" analysis format=10.5 ;
define _BIC_/ "BIC (smaller is better)" analysis format=10.5 ;
run;
ods printer close;
ods listing;

/*Plot Figure C (Residual with predicted value) */;
* output residuals of model with continuous, categorical and spline term, respectively;
proc reg data=Data_prep;
  model &yvar.= &xvar_cont. &covarlist_cat. &covarlist_cont.;
  output out=rplot_c p=pred r=res;run;
proc reg data=Data_prep;
  model &yvar.= &xvar_cat. &covarlist_cat. &covarlist_cont.;
  output out=rplot_cat p=pred r=res;run;
proc reg data=Data_prep;
  model &yvar.= &xvar_cont. &xvar_cont.1 -- &xvar_cont.%eval(&knot-2) &covarlist_cat. &covarlist_cont.;
  output out=rplot_sp p=pred r=res;run;

data rplot_c;length Form $12.;set rplot_c(keep=&xvar_cont. res);Form='Continuous';run;
data rplot_cat;length Form $12.;set rplot_cat(keep=&xvar_cont. res);Form='Categorical';run;
data rplot_sp;length Form $12.;set rplot_sp(keep=&xvar_cont. res);Form='Spline';run;
proc sort data=rplot_cat;by &xvar_cont.;run;
proc sort data=rplot_c;by &xvar_cont.;run;
proc sort data=rplot_sp;by &xvar_cont.;run;
data rplot_all;	set rplot_c rplot_cat rplot_sp;	by &xvar_cont.;run;

* Plot residual with value of independent variable from linear model;
ods listing gpath="&dataout.";
ods  graphics on / reset=index imagefmt=png imagename="FigC - &xvar_cont." ;
proc sgpanel data =rplot_all;	
panelby Form/columns=3;
scatter x=&xvar_cont. y=res;
colAXIS values=(&minvalue. to &maxvalue.) LABEL="&xvar_cont." ; 			
refline 0 /transparency=0.5 axis=y;run;
ods graphics off;
ods printer close;

/* Plot Figure D (Scatter Plot of data distribution) */;
ods listing gpath="&dataout.";
ods graphics on / reset=index imagefmt=png imagename="FigD - &xvar_cont." ;
proc sgplot data = Data_fig; 
scatter x = &xvar_cont. y=&yvar.;
XAXIS values=(&minvalue. to &maxvalue.) LABEL="&xvar_cont." ; 			
run;
ods graphics off;
ods printer close;
%MEND SPECI_LINEAR;


******************************************************************************************************;
*******Logistic Regression Macro ****;
******************************************************************************************************;

%MACRO SPECI_LOGISTIC;
* Run separate logistic regression models with exposure variable as continous, categorical and spline;
ods output  Rsquare=est1line_sq  FitStatistics=est1line_fit GlobalTests=est1line_glob convergencestatus=est1line_con association=est1line_c;
proc logistic data=data_prep descending outest=est1Line;
   class &covarlist_cat. / param=glm;
  model &yvar.= &xvar_cont. &covarlist_cat. &covarlist_cont. /rsquare influence;run;quit;
data data_prep;
	set data_prep;
  CALL SYMPUT('newref_xvar_cat', PUT(newref, 3.));
run;	
%let newcat_macro=newcat;
ods output  Rsquare=est1cat_sq  FitStatistics=est1cat_fit GlobalTests=est1cat_glob convergencestatus=est1cat_con association=est1cat_c;
proc logistic data=data_prep descending outest=est1cat;
   class &newcat_macro.(ref="&newref_xvar_cat.") &covarlist_cat. / param=glm;
  model &yvar.= &newcat_macro. &covarlist_cat. &covarlist_cont./rsquare influence; run; quit;
ods output  Rsquare=est1Spline_sq  FitStatistics=est1Spline_fit GlobalTests=est1Spline_glob convergencestatus=est1spline_con association=est1spline_c;
proc logistic data = data_prep descending outest=est1Spline;
   class &covarlist_cat. / param=glm;
   model &yvar. = &xvar_cont. &xvar_cont.1 -- &xvar_cont.%eval(&knot-2) &covarlist_cat. &covarlist_cont./rsquare influence;run; quit;								

** Get estimates from modeling ;
* continuous;
data Line(keep=int_line linpred temp);
   set est1line(rename=(intercept=int_line));
   LinPred = &xvar_cont.;temp = 1;run;
* categorical;
data Cat(keep=int_cat &newcat_macro.1-&newcat_macro.%eval(&num_cat) temp); 
   set est1cat(rename=(intercept=int_cat));
   temp = 1;run;
* spline;
data Spline(keep=int_sp sp sp1--sp%eval(&knot-2) temp);																		
   set est1Spline (rename=(intercept=int_sp &xvar_cont.=sp));
   array a{%eval(&knot-2)} sp1-sp%eval(&knot-2);
   array b(%eval(&knot-2)) &xvar_cont.1-&xvar_cont.%eval(&knot-2);
   do i =1 to %eval(&knot-2);  	a(i)=b(i);	end; temp = 1;run;																								

data Data_est;set Line; set Cat;set Spline;run;

data Data_plot;
   merge Data_prep Data_est;
   by temp;run;

** Calculate Predicted values for logit(Y) for continuous, categorical and spline models;
%macro cat;
data data_logcat;
	set Data_plot;
	%do j=1 %to &num_cat.;
	/* model with exposure variable as categorical;*/
	if &newcat_macro.=&j. then Log_cat=int_cat+&newcat_macro.&j.; /* Logit(P) = a + b(&xvar_cat.&j.), if &xvar_cat.=&j. */
	%end;run;
%mend cat;
%cat;quit;

data Data_final(drop=k sum);
   set data_logcat;
   * model with exposure variable as continuous ;
   Log_Linear =  (int_line+(&xvar_cont.*LinPred));		/* Logit(P) = a + b(&xvar_cont.)*x(&xvar_cont.) */
   * model with exposure variable as spline;
   sum=0;
   array c{%eval(&knot-2)} sp1-sp%eval(&knot-2);
   array d(%eval(&knot-2)) &xvar_cont.1-&xvar_cont.%eval(&knot-2);
   do k =1 to %eval(&knot-2);
   sum=sum+c(k)*d(k);
   end;
   Log_Spline = (int_sp+(&xvar_cont.*sp) + sum);		/* Logit(P) = a + b(&xvar_cont.)*x(&xvar_cont.) + b(sp1)*x(&xvar_cont.) + b(sp2)*x(&xvar_cont.) + .... b(spn)*x(&xvar_cont.)*/											
retain;run;								

* Centralize predictions at mid-point to align predicted plots;
* Calculate the mid-point value of exposure variable of interest by ((max - min)/2 + min);
PROC SQL;
	CREATE TABLE mid AS
	SELECT temp, Log_Linear as mid_linear, Log_cat as mid_cat,Log_spline as mid_spline, min(&xvar_cont.)as min, max(&xvar_cont.) as max, 
		  ((calculated max)-(calculated min))/2 + (calculated min) as mid, abs((&xvar_cont. - (calculated mid))) as abs_diff
	FROM Data_final	order by abs_diff;QUIT;

* Keep predicted Logit(P) at mid-point;
proc sort data=mid(firstobs=1 obs=1) out=mid1;by temp;run;


*  calculate the centralized value by subtracting midpoint from each predicted value;
data data_fig;
	merge Data_final mid1;
	by temp;
	log_mid_line=Log_Linear-mid_linear;
	log_mid_cat=Log_Cat-mid_cat;
	log_mid_spline=Log_Spline-mid_spline;
	* Assign min and max value as new macro variables so as to set up x axis range in the plot;
	call symput('minvalue',min);
	call symput('maxvalue',max);

	label log_mid_line='Continuous';
	label log_mid_cat='Categorical';
	label log_mid_spline='Spline';
run;


* order dataset by exposure ;
proc sort data=data_fig out=Data_fig;by &xvar_cont.;run;

/* Plotting Figure A (Predicted logit(P) with exposure X)*/
ods listing gpath="&dataout.";
ods  graphics on /reset=index imagefmt=png imagename="FigA - &xvar_cont." ;
proc sgplot data =Data_fig;											
series x=&xvar_cont. y=log_mid_line/LINEATTRS = (color=ROSE THICKNESS = 4);
series x=&xvar_cont. y=log_mid_cat/LINEATTRS = (color=o  THICKNESS = 4);
series x=&xvar_cont. y=log_mid_spline/LINEATTRS = (color=BIGB THICKNESS = 4);
YAXIS LABEL ="Centralized Logit(P(&yvar.=1))" ; 
XAXIS values=(&minvalue. to &maxvalue.) LABEL="&xvar_cont." ; 			
run;
ods graphics off;
ods printer close;

/*Figure B (Summary Table for Comparing Model Statistics) */;
* Integrate Statistics - Model with exposure variable as continuous, categorical and spline;
%macro form(form);
data &form._sq1;set est1&form._sq (keep=label1 nvalue1 rename=(label1=label nvalue1=&form.));run;
data &form._sq2;set est1&form._sq (keep=label2 nvalue2 rename=(label2=label nvalue2=&form.));run;
data &form._fit1;set est1&form._fit (keep=Criterion InterceptAndCovariates rename=(Criterion=label InterceptAndCovariates=&form.));if label^='SC';run;
data &form._glob1;set est1&form._glob (keep=test ProbChiSq rename=(test=label ProbChiSq=&form.)); if label^='Score';run;
data &form._con1;set est1&form._con (keep=reason status rename=(reason=label status=&form.));run;
data &form._c(where=(label='c'));set Est1&form._c(keep=label2 nvalue2 rename=(label2=label nvalue2=&form.));run;
data &form._stat;set &form._sq1 &form._sq2 &form._c &form._fit1 &form._glob1 &form._con1;run;
proc transpose data=&form._stat out=&form._tran(drop=_LABEL_);run;
%mend form;
%form(line)
%form(cat)
%form(spline)

data line_tran;length _name_ $13;set line_tran; _NAME_="Continuous";run;
data cat_tran;length _name_ $13;set cat_tran; _NAME_="Categorical";run;
data spline_tran;length _name_ $13;set spline_tran; _NAME_="Spline";run;

data table;	set line_tran cat_tran spline_tran;run;

* Report Summary Table;
options printerpath=png nodate papersize=('8in','5in');
ods _all_ close;
ods printer file="&dataout.\FigB - &xvar_cont..png";
Proc report data=table nowd ;
column _name_ col1 col2 col3 col4 col5 col6 col7 col8;
define _name_ /"Diagnostic Statistics" group order=data ;
define col1/ "R-Squared (bigger is better)" analysis format=10.5 ;
define col2/ "Max-rescaled R-Squared (bigger is better)" analysis format=10.5 ;
define col3/ "C-Statistics (bigger is better)" analysis format=10.5 ;
define col4/ "AIC (smaller is better)" analysis format=10.5 ;
define col5/ "-2LogL (bigger is better)" analysis format=10.5 ;
define col6/ "Likelihood Test (P-value)" analysis format=10.4 ;
define col7/ "Wald Test (P-value)" analysis format=10.4 ;
define col8/ "Model Convergence(0=Yes, 1=No)" analysis format=1.0 ;
run;
ods printer close;
ods listing;

/*Plot Figure C (Plotting Pearson Residual with value of exposure) */;
* output Pearson residual of model with continous, categorical and spline form;
proc logistic data=Data_prep descending outest=est1Line;
  class &covarlist_cat.  / param=glm;
  model &yvar.= &xvar_cont. &covarlist_cat. &covarlist_cont.;
  output out=rplot_c prob=p reschi=pr;run;
proc logistic data=Data_prep descending outest=est1Line;
  class &xvar_cat. &covarlist_cat.  / param=glm;
  model &yvar.= &xvar_cat. &covarlist_cat. &covarlist_cont.;
  output out=rplot_cat prob=p reschi=pr;run;
proc logistic data=Data_prep descending outest=est1Line;
  class &covarlist_cat.  / param=glm;
  model &yvar.= &xvar_cont. &xvar_cont.1 -- &xvar_cont.%eval(&knot-2) &covarlist_cat. &covarlist_cont.;
  output out=rplot_sp prob=p reschi=pr;run;
  
data rplot_cat;length Form $12.;set rplot_cat(keep=&xvar_cont. pr);Form='Categorical';run;
data rplot_c;length Form $12.;set rplot_c(keep=&xvar_cont. pr);Form='Continuous';run;
data rplot_sp;length Form $12.;set rplot_sp(keep=&xvar_cont. pr);Form='Spline';run;

proc sort data=rplot_cat;by &xvar_cont.;run;
proc sort data=rplot_c;by &xvar_cont.;run;
proc sort data=rplot_sp;by &xvar_cont.;run;

data rplot_all;	set rplot_c rplot_cat rplot_sp;	by &xvar_cont.;run;

* Plotting Pearson Residual with value of exposure from continuous model;
ods listing gpath="&dataout.";
ods  graphics on /reset=index imagefmt=png imagename="FigC - &xvar_cont." ;
proc sgpanel data =rplot_all;	
panelby Form/columns=3;
scatter x=&xvar_cont. y=pr;
colAXIS values=(&minvalue. to &maxvalue.) LABEL="&xvar_cont." ; 			
refline 0 /transparency=0.2 axis=y;
refline 3 -3 /transparency=0.6 axis=y;run;
ods graphics off;
ods printer close;

/*Plot Figure D (Scatter Plot of data distribution ) */;
ods listing gpath="&dataout.";
ods graphics on /reset=index imagefmt=png imagename="FigD - &xvar_cont." ;
proc sgplot data = Data_fig; 
scatter x = &xvar_cont. y=&yvar.; 
yaxis values=(0 to 1 by 1);
XAXIS values=(&minvalue. to &maxvalue.) LABEL="&xvar_cont." ; 			
run;
ods graphics off;
ods printer close;
%MEND SPECI_LOGISTIC;



******************************************************************************************************;
*******Survival Macro ****;
******************************************************************************************************;

%MACRO SPECI_SURVIVAL;
* Run separate survival models with exposure variable as continous, categorical and spline;
ods output  FitStatistics=est1line_fit GlobalTests=est1line_glob convergencestatus=est1line_con censoredsummary=est1line_censor;
proc phreg data = data_prep outest=est1Line;
   class &covarlist_cat. / param=glm;
   model &time2event.*&EVENT.(0) =&xvar_cont. &covarlist_cat. &covarlist_cont.;
run;quit;
data data_prep;
	set data_prep;
  CALL SYMPUT('newref_xvar_cat', PUT(newref, 3.));
run;	
%let newcat_macro=newcat;
ods output  FitStatistics=est1cat_fit GlobalTests=est1cat_glob convergencestatus=est1cat_con censoredsummary=est1cat_censor;
proc phreg data=data_prep outest=est1cat;
  class &newcat_macro.(ref="&newref_xvar_cat.") &covarlist_cat. / param=glm;
  model &time2event.*&EVENT.(0)= &newcat_macro. &covarlist_cat. &covarlist_cont.;
run;quit;
ods output  FitStatistics=est1Spline_fit GlobalTests=est1Spline_glob convergencestatus=est1spline_con censoredsummary=est1spline_censor;
proc phreg data=data_prep outest=est1Spline;
   class &covarlist_cat. / param=glm;
  model &time2event.*&EVENT.(0)= &xvar_cont. &xvar_cont.1 -- &xvar_cont.%eval(&knot-2) &covarlist_cat. &covarlist_cont.;
run;quit;


** Get estimates from modeling ;
* continuous;
data Line(keep=linpred temp);
   set est1line;LinPred = &xvar_cont.;temp = 1;run;
* Categorical;
%macro rename;
data data_renamecat;
	set est1cat;
	%do j=1 %to &num_cat.;
	rename &newcat_macro.&j.=&xvar_cat.new&j. ; 
	%end;run;
%mend rename;
%rename;

data Cat (keep=&xvar_cat._newcat1-&xvar_cat._newcat%eval(&num_cat) temp); 
   set data_renamecat;
   array a{%eval(&num_cat)} &xvar_cat._newcat1-&xvar_cat._newcat%eval(&num_cat);
   array b(%eval(&num_cat)) &xvar_cat.new1-&xvar_cat.new%eval(&num_cat);
   do i =1 to %eval(&num_cat);  	a(i)=b(i);	end; temp=1;
run;

* Spline;
data Spline(keep=sp sp1--sp%eval(&knot-2) temp);																		
   set est1Spline (rename=(&xvar_cont.=sp));
   array a{%eval(&knot-2)} sp1-sp%eval(&knot-2);
   array b(%eval(&knot-2)) &xvar_cont.1-&xvar_cont.%eval(&knot-2);
   do i =1 to %eval(&knot-2);  	a(i)=b(i);	end; temp = 1;run;																								

data Data_est;set Line; set Cat;set Spline;run;

data Data_plot;
   merge Data_prep(in=a) Data_est(in=b);
   by temp;
run;

** Calculate Predicted values for log(HR) for continuous, categorical and spline models;
%macro cat;
data data_logcat;
	set Data_plot;
	%do j=1 %to &num_cat.;
    * model with exposure variable as categorical;
	if &newcat_macro.=&j. then Log_cat=&xvar_cat._newcat&j.; /* Log(HR) = b(&xvar_cat.&j.), if &xvar_cat.=&j. */
	%end;run;
%mend cat;
%cat;

data Data_final(drop=k sum);
   set data_logcat;
    * model with exposure variable as continuous ;
   Log_Linear =  (&xvar_cont.*LinPred);		/* Log(HR) = b(&xvar_cont.)*x(&xvar_cont.) */
   * model with exposure variable as spline;
   sum=0;
   array c{%eval(&knot-2)} sp1-sp%eval(&knot-2);
   array d(%eval(&knot-2)) &xvar_cont.1-&xvar_cont.%eval(&knot-2);
   do k =1 to %eval(&knot-2);
   sum=sum+c(k)*d(k);
   end;
   Log_Spline = ((&xvar_cont.*sp) + sum);		/* Log(HR) = b(&xvar_cont.)*x(&xvar_cont.) + b(sp1)*x(&xvar_cont.) + b(sp2)*x(&xvar_cont.) + .... b(spn)*x(&xvar_cont.)*/											
retain;run;		

* Centralize predictions at mid-point to align predicted plots;
* Calculate the mid-point value of exposure variable of interest by ((max - min)/2 + min);
PROC SQL;
	CREATE TABLE mid AS
	SELECT temp, Log_Linear as mid_linear, Log_cat as mid_cat,Log_spline as mid_spline, min(&xvar_cont.)as min, max(&xvar_cont.) as max, 
		  ((calculated max)-(calculated min))/2 + (calculated min) as mid, abs((&xvar_cont. - (calculated mid))) as abs_diff
	FROM Data_final	order by abs_diff;QUIT;

	* Keep predicted log(HR) at mid-point;
proc sort data=mid(firstobs=1 obs=1) out=mid1;by temp;run;

*  calculate the centralized value by subtracting midpoint from each predicted value;
data data_fig;
	merge Data_final mid1;
	by temp;
	log_mid_line=Log_Linear-mid_linear;
	log_mid_cat=Log_Cat-mid_cat;
	log_mid_spline=Log_Spline-mid_spline;
	* Assign min and max value as new macro variables so as to set up x axis range in the plot;
	call symput('minvalue',min);
	call symput('maxvalue',max);

	label log_mid_line='Continuous';
	label log_mid_cat='Categorical';
	label log_mid_spline='Spline';
run;

* Sort dataset by exposure variable;
proc sort data=data_fig out=Data_fig;by &xvar_cont.;run;

/* Plot Figure A (Predicted Log(HR) with X) */;
ods listing gpath="&dataout.";
ods  graphics on /reset=index imagefmt=png imagename="FigA - &xvar_cont." ;
proc sgplot data =Data_fig;											
series x=&xvar_cont. y=log_mid_line/LINEATTRS = (color=ROSE THICKNESS = 4);
series x=&xvar_cont. y=log_mid_cat/LINEATTRS = (color=o  THICKNESS = 4);
series x=&xvar_cont. y=log_mid_spline/LINEATTRS = (color=BIGB THICKNESS = 4);
YAXIS LABEL ="Centralized Log(HR(&event.=1))" ; 
XAXIS values=(&minvalue. to &maxvalue.) LABEL="&xvar_cont." ; 			
run;
ods graphics off;
ods printer close;

/*Figure B (Summary Table for Comparing Model Statistics) */;
* Integrate Statistics - Model with exposure variable as continuous, categorical and spline;
%macro form(form);
data &form._fit1;set est1&form._fit (keep=Criterion WithCovariates rename=(Criterion=label WithCovariates=&form.));run;
data &form._glob1;set est1&form._glob (keep=test ProbChiSq rename=(test=label ProbChiSq=&form.)); if label^='Score';run;
data &form._con1;set est1&form._con (keep=reason status rename=(reason=label status=&form.));run;
data &form._stat;set &form._fit1 &form._glob1 &form._con1;run;
proc transpose data=&form._stat out=&form._tran(drop=_LABEL_);run;
%mend form;
%form(line)
%form(cat)
%form(spline)

data line_tran;length _name_ $13;set line_tran; _NAME_="Continuous";run;
data cat_tran;length _name_ $13;set cat_tran; _NAME_="Categorical";run;
data spline_tran;length _name_ $13;set spline_tran; _NAME_="Spline";run;

data table;	set line_tran cat_tran spline_tran;run;

* Report Summary Table;
options printerpath=png nodate papersize=('8in','5in');
ods _all_ close;
ods printer file="&dataout.\FigB - &xvar_cont..png";
Proc report data=table nowd ;
column _name_ col1 col2 col3 col4 col5 col6;
define _name_ /"Covariate Specification" group order=data ;
define col1/ "-2LogL (bigger is better)" analysis format=10.5 ;
define col2/ "AIC (smaller is better)" analysis format=10.5 ;
define col3/ "SBC/BIC (smaller is better)" analysis format=10.5 ;
define col4/ "Likelihood Test (P-value)" analysis format=10.4 ;
define col5/ "Wald Test (P-value)" analysis format=10.4 ;
define col6/ "Model Convergance(0=Yes, 1=No)" analysis format=1.0 ;
run;
ods printer close;
ods listing;

/* Plot Figure C (Deviance Residual with X) */;
* output deviance residual of model with continuous, categorical and spline terms;
proc phreg data=Data_fig;
  class &covarlist_cat.  / param=glm;
  model &time2event.*&EVENT.(0)= &xvar_cont. &covarlist_cat. &covarlist_cont./ties=efron;
  output out=rplot_c resdev=dev;run;
proc phreg data=Data_fig;
  class &xvar_cat. &covarlist_cat. / param=glm;
  model &time2event.*&EVENT.(0)= &xvar_cat. &covarlist_cat. &covarlist_cont./ties=efron;
  output out=rplot_cat resdev=dev;run;
proc phreg data=Data_fig;
  class &covarlist_cat.  / param=glm;
  model &time2event.*&EVENT.(0)= &xvar_cont. &xvar_cont.1 -- &xvar_cont.%eval(&knot-2) &covarlist_cat. &covarlist_cont./ties=efron;
  output out=rplot_sp resdev=dev;run;

data rplot_cat;length Form $12.;set rplot_cat(keep=&xvar_cont. dev);Form='Categorical';run;
data rplot_c;length Form $12.;set rplot_c(keep=&xvar_cont. dev);Form='Continuous';run;
data rplot_sp;length Form $12.;set rplot_sp(keep=&xvar_cont. dev);Form='Spline';run;

proc sort data=rplot_cat;by &xvar_cont.;run;
proc sort data=rplot_c;by &xvar_cont.;run;
proc sort data=rplot_sp;by &xvar_cont.;run;

data rplot_all;	set rplot_c rplot_cat rplot_sp;	by &xvar_cont.;run;

* Plot Devance Residual with predicted value of exposure as continuous;
ods listing gpath="&dataout.";
ods  graphics on /reset=index imagefmt=png imagename="FigC - &xvar_cont." ;
proc sgpanel data =rplot_all;	
panelby Form/columns=3;
scatter x=&xvar_cont. y=dev;
colAXIS values=(&minvalue. to &maxvalue.) LABEL="&xvar_cont." ; 			
refline 0 /transparency=0.2 axis=y;
refline 3 -3 /transparency=0.6 axis=y;run;
ods graphics off;
ods printer close;

/*Plot Figure D (Kaplan Meier Curve) */;
ods listing gpath="&dataout.";
ods graphics on /reset=index imagefmt=png imagename="FigD - &xvar_cont." ;
proc lifetest data=Data_fig plots=survival(atrisk(maxlen=13 outside(0.05)));
   time &time2event.*&EVENT.(0);
   run;
ods graphics off;
ods printer close;

%MEND SPECI_SURVIVAL;



*******************************************************************************************************;
*** Part III - report output;
*******************************************************************************************************;

%MACRO REPORT;
*Generate One-Pager with four plots together;
options orientation=landscape nodate nonumber;  
goptions reset=all device=png1000 nodisplay xmax=4in ymax=3in;

goptions iback="&dataout.\figa - &xvar_cont..png" imagestyle=fit;
proc gslide;run;quit;
goptions iback="&dataout.\figb - &xvar_cont..png" imagestyle=fit;
proc gslide;run;quit;
goptions iback="&dataout.\figc - &xvar_cont..png" imagestyle=fit;
proc gslide;run;quit;
goptions iback="&dataout.\figd - &xvar_cont..png" imagestyle=fit;
proc gslide;run;quit;

goptions reset=all device=sasprtc; 
ods listing close; 
ods pdf file="&dataout.\&reportname..pdf" notoc dpi=1500;
proc greplay igout=work.gseg nofs tc=sashelp.templt template=L2R2;
    treplay 1:1 2:3 3:2 4:4;
run;quit;
ods pdf close; 

ods listing; 
* Delete all of the graphs from the WORK.GSEG catalog;
 proc catalog c=work.gseg kill;run;

 %MEND REPORT;


