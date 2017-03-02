# A SAS® Macro for Covariate Specification in Linear, Logistic, or Survival Regression

This website will introduce a customizable user-friendly SAS® macro %SPECI to quickly produce a one page report that organizes multiple commonly-used statistics to help you compare and select the appropriate functional form from continuous, categorical, and spline terms in linear regression, logistic regression, and survival analysis models.  

The statistics in the final report include:

•	Plot showing an overlay of predicted values from the three functional forms. 

•	Summary table of model statistics. (See complete list and descriptions for each model in Appendix A) 

•	Plot of the residual values from the model where the covariate is continuous. 

•	Plot of the observed values of the covariate and the outcome variable (linear and logistic regression models only)

•	Kaplan Meier plot (survival model only). 


# INSTRUCTIONS FOR USING MACRO %SPECI 

There are two SAS® editor programs: the main macro (SPECI.sas) and the program to call the macro (CALL SPECI.sas). The call program is provided in the Appendix B and both the main macro and the call program are available upon request from the author (Sai Liu) and are posted to the GitHub website (https://github.com/SaiLMainpage/ModelSpecification).

• First, save the CALL SPECI.sas and SPECI.sas programs to your computer. Next, open “CALL SPECI.sas” and update the include statement to the directory where the “SPECI.sas” macro stored 

   %include "Directory/speci.sas";

• Next, specify the parameters for the macro program (for example %let dataset= mydata) see CALL SPECI.sas program.

     
# CONTACT INFORMATION 
  Your comments and questions are valued and encouraged. Contact the author at:

  # Sai Liu

  Division of Nephrology, Department of Medicine 
  Stanford University School of Medicine
  1070 Arastradero Rd., Suite 100
  Palo Alto, CA. 94304
  Phone: 213-793-1055
  Email: sailiu.tian@gmail.com
