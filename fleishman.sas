
*****   Program CHAPTER1_1_NONNORMAL_FLEISHMAN.SAS  Deriving Fleishman Coefficients for Desired Skewness and Kurtosis  *****;

%macro fleishman;
   /*	This program calculates the coefficients for Fleishman's power transformation in order     */
   /*   to obtain univariate non-normal variables.  For references, see Allen I. Fleishman, (1978).*/
   /*   A method for simulating non-normal distributions, Psychometrika, 43, 521-532.  Also see    */
   /*   Vale, C. David and Maurelli, Vincent A.  (1983).  Simulating multivariate non-normal       */
   /*   distributions, Psychometrika, 48, 465-471.                                                 */                   

PROC IML;

 use skewkurt; 
 read all var{skewness kurtosis} into skewkurt; 

START NEWTON;
  RUN FUN;
  DO ITER = 1 TO MAXITER
  WHILE(MAX(ABS(F))>CONVERGE);
        RUN DERIV;
        DELTA=-SOLVE(J,F);
        COEF=COEF+DELTA;
        RUN FUN;
  END;
FINISH NEWTON;
MAXITER=25;
CONVERGE=.000001;
START FUN;
  X1=COEF[1];
  X2=COEF[2];
  X3=COEF[3];
  F=(X1**2+6*X1*X3+2*X2**2+15*X3**2-1)//
    (2*X2*(X1**2+24*X1*X3+105*X3**2+2)-SKEWNESS)//
    (24*(X1*X3+X2**2*(1+X1**2+28*X1*X3)+X3**2*
      (12+48*X1*X3+141*X2**2+225*X3**2))-KURTOSIS);
FINISH FUN;
START DERIV;
  J=((2*X1+6*X3)||(4*X2)||(6*X1+30*X3))//
    ((4*X2*(X1+12*X3))||(2*(X1**2+24*X1*X3+105*X3**2+2))
     ||(4*X2*(12*X1+105*X3)))//
    ((24*(X3+X2**2*(2*X1+28*X3)+48*X3**3))||
     (48*X2*(1+X1**2+28*X1*X3+141*X3**2))||
     (24*(X1+28*X1*X2**2+2*X3*(12+48*X1*X3+141*X2**2+225*X3**2)
 
     +X3**2*(48*X1+450*X3))));
FINISH DERIV;
DO;
NUM = NROW(SKEWKURT);
DO VAR=1 TO NUM;
  SKEWNESS=SKEWKURT[VAR,1];
  KURTOSIS=SKEWKURT[VAR,2];
  COEF={1.0, 0.0, 0.0};
  RUN NEWTON;
  COEF=COEF`;
  SK_KUR=SKEWKURT[VAR,];
  COMBINE=SK_KUR || COEF;
  IF VAR=1 THEN RESULT=COMBINE;
  ELSE IF VAR>1 THEN RESULT=RESULT // COMBINE;
END;
  PRINT "COEFFICEINTS OF B, C, D FOR FLEISHMAN'S POWER TRANSFORMATION";
  PRINT "Y = A + BX + CX^2 + DX^3";
  PRINT " A = -C";
  MATTRIB RESULT COLNAME=({SKEWNESS KURTOSIS B C D})
                 FORMAT=12.9;
  PRINT RESULT;
END;
 create fleishman from result[colname={SKEWNESS KURTOSIS B C D}];
   append from result;

QUIT;
 
%mend fleishman;

***************************************************************************;

/*  EXAMPLE 1: obtaining coefficients for 4 variables */

/*	In the following matrix 'SKEWKURT', specify the skewness and kurtosis for each variable.  */
/*  Each row represents one variable. In each row, the 1st number is the skewness, the 2nd is */
/*  the kurtosis of the variable;                                                             */

data skewkurt;
 input skewness kurtosis;
 datalines;
 -2 6
 -2 11
 0 0
 0 6
 0 11
 2 6
 2 11
;
run;

%fleishman; 

