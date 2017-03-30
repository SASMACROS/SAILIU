# A SAS® Macro for Covariate Specification in Linear, Logistic, and Survival Regression

This website will introduce a customizable user-friendly SAS® macro %SPECI to quickly produce a one page report that organizes multiple commonly-used statistics to help you compare and select the appropriate functional form from continuous, categorical, and spline terms in linear regression, logistic regression, and survival analysis models.  


The statistics in the final report include:

•	Predicted Plot showing predicted values for each of the three functional forms. 

•	Summary Diagnostic Statistics Table (See complete list and descriptions for each model in Appendix A) 

•	Residual Plot - for model with the covariate in continuous, categorical and spline terms. 

•	Correlation Plot - Distribution of outcome variable by the covariate in continuous form (linear and logistic regression models only)

•	Kaplan Meier plot (survival model only). 


# INSTRUCTIONS FOR USING MACRO %SPECI 

There are two SAS® editor programs: the main macro (SPECI.sas) and the program to call the macro (CALL SPECI.sas). The SAS proceeding paper of this project with more instructions and details are available upon request from the author (Sai Liu) and will be posted on the website (https://www.sas.com/en_us/events/sas-global-forum/sas-global-forum-2017.html, search "A SAS® Macro for Covariate Specification in Linear, Logistic, and Survival Regression" or "Paper 1223").

• First, save the CALL SPECI.sas and SPECI.sas programs to your computer; 

• Next, open “CALL SPECI.sas” and update the include statement (%include "Directory/speci.sas") to the directory where the “SPECI.sas” macro stored;

• Lastly, specify the parameters for the macro program (for example %let dataset= mydata) see CALL SPECI.sas program.

     
# CONTACT INFORMATION 
  Your comments and questions are valued and encouraged. Contact the author at:

  Sai Liu

  Division of Nephrology, Department of Medicine 
  Stanford University School of Medicine
  1070 Arastradero Rd., Suite 100
  Palo Alto, CA. 94304
  Email: sailiu.tian@gmail.com
