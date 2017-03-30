
/********************************************************************************************************/
/*       NAME: CALL SPECI.SAS                                         				  	*/
/*      TITLE: Functional form Specification for Linear, Logistic, and Survial Models	  		*/
/*     AUTHOR: Sai Liu, MPH, Stanford University                                          		*/
/*		   OS: Windows 7 Ultimate 64-bit							*/
/*	 Software: SAS 9.4										*/
/*       DATE: 09 FEB 2017                                        					*/
/*DESCRIPTION: This program shows how to call the SPECI.sas macro					*/
/********************************************************************************************************/


/************************************************************************************************************
* Now, please set up your Parameters below
*************************************************************************************************************/
%let datain=" "; 				/*Location of permanent SAS dataset. Leave it blank, if your dataset is in the work library*/	
%let dataout=" ";				/*Location of one-pager report will be saved*/												
%let dataset=;					/*Name of the dataset*/					
%let reportname=;				/*Name of the one-pager report
						  If you leave it blank, this program will give a name as "Model Diagnostic Report" */

%let model=;	           			/*1=linear, 2=logistic, 3=survival*/

%let yvar=;					/*Dependent variable of interest for linear and logistic regression model only, otherwise leave it blank
						  e.g. %let yvar= heartfail; or %let yvar= ;(if this is not for linear nor logistic regression model)*/

%let event=;					/*Dependent variable of interest for survival model only, otherwise leave it blank 
						  e.g. %let event= death; or %let event= ;(if this is not for survival model)*/
%let time2event=;	        		/*Dependent variable time component for survival models only, otherwise leave it blank 
						  e.g. %let time2event= time2death; or %let time2event= ;(if this is not for survival model)*/

%let xvar_cont=;				/*Independent variable of interest (continuous). e.g. %let xvar_cont= age; */								
%let xvar_cat=;					/*Independent variable of interest (categorical). 
						  If you don't have an exist categorical variable in the dataset, 
						  please create a categorical variable and entry the name of created variable here and leave datain=   blank,
                              			  because the main dataset is already in work library, this program won't read a permanent SAS dataset. 
						  e.g. %let xvar_cat= bmi_cat; */	
	
%let num_cat=;					/*# of categories of above categorical variable. MUST enter a numeric number, can't leave it blank.
						         e.g. bmi_cat has 4 categories, then %let num_cat= 4;*/
	
%let ref_xvar_cat=;				/*Specify the reference group of the xvar_cat variable. 
							  e.g. the 2nd category of bmi_cat is the reference, then %let ref_xvar_cat= 2; 
							       If you leave it blank, this program will set 1st category as the reference */
	
%let covarlist_cont=;				/*continuous covariates for model adjustment. 
							  OPTIONS: 1) Leave it blank if no continuous covariate to be included in the model or
								   2) add continuous covariates as needed and add one space between multicovariates. 
							  			  e.g. %let covarlist_cont= age height weight; */		
	
%let covarlist_cat=;				/*categorical covariates for model adjustment. Need one space between multicovariates. 
							  OPTIONS: 1) Leave it blank if no categorical covariate to be included in the model or
								   2) add categorical covariates as needed and add one space between multicovariates. 
							  			  e.g. %let covarlist_cat= race year; */		
	
%let knot=;					/*# of Knots for Spline. MUST enter a number from 4 to 10. 
						  Because NO SPLINE VARIABLES CREATED if number of knots <=3 in this program
						  Default is 4*/	
%let norm=;					/*Normalization method. Options: 0, 1, or 2. Default is 2*/
%let knot1=;					/*Placement of percentile Cutoff at 1st knot e.g. %let knot1= p5 (p5 for 5th percentile)*/												
%let knot2=;					/*Placement of percentile Cutoff at 2nd knot e.g. %let knot1= p35 (p35 for 35th percentile)*/												
%let knot3=;					/*Placement of percentile Cutoff at 3rd knot e.g. %let knot1= p65 (p65 for 65th percentile)*/												
%let knot4=;					/*Placement of percentile Cutoff at 4th knot e.g. %let knot1= p95 (p95 for 95th percentile)*/												
%let knot5=;					/*Placement of percentile Cutoff at 5th knot. Leave it blank even though you don't specify*/												
%let knot6=;					/*Placement of percentile Cutoff at 6th knot. Leave it blank even though you don't specify*/												
%let knot7=;					/*Placement of percentile Cutoff at 7th knot. Leave it blank even though you don't specify*/												
%let knot8=;					/*Placement of percentile Cutoff at 8th knot. Leave it blank even though you don't specify*/												
%let knot9=;					/*Placement of percentile Cutoff at 9th knot. Leave it blank even though you don't specify*/												
%let knot10=;					/*Placement of percentile Cutoff at 10th knot. Leave it blank even though you don't specify*/
/************************************************************************************************************/

%include "Directory\speci.sas";

/*
libname lib "";
data data;
	set lib.&dataset.;

	* creating categorical variable;

run;
*/

%SPECI;
Quit;

