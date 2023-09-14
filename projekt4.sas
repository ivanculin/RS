/*	MACRO koji generira replikacije uzoraka		*/
/*	i racuna snagu jakosti testa za uzorke koji	*/
/*	dolaze iz populacija jednakih varijanci		*/

/*	za varijable macro uzima nizove duljina prvog	*/
/*	i drugog uzorka (nn1 i nn2), uzima standardnu 	*/
/*	devijaciju uzoraka (SIG), pomocnu varijablu 	*/
/*	za oznake data setova (skewness i kurtosis) */
%MACRO GENERATOR(N1,N2,SIG_1,SIG_2,SKEW_1,KURT_1);

%LET K1=%SYSEVALF(&SIG_1**2);
%LET K2=%SYSEVALF(&SIG_2**2);




%fleishman;
DATA GEN (KEEP= SAMPLEID Cc X DELTA);
    set fleishman(where=(skewness=&skew_1 and kurtosis=&kurt_1));
	a = -c;		
	DO DELTA=-6 TO 6 BY 2;
		DO SAMPLEID=1 TO &NUMSAMPLES;
			Cc=1;
			DO I=1 TO &N1;
				CALL STREAMINIT(&SEED1);
				X=RAND("NORMAL");
				X = a+b*x+c*x**2+d*x**3;
				x = delta + &sig_1*x;
				OUTPUT;
					
			END;
			Cc=2;
			DO I=1 TO &N2;
				CALL STREAMINIT(&SEED2);
				x=RAND("NORMAL");
				x = a+b*x+c*x**2+d*x**3;
				x = &sig_2*x;
				OUTPUT;
			END;
		END;
	END;
RUN;
%MEND GENERATOR;
/*	*	*	*	*	*	*	*/
/*	*	*	*	*	*	*	*/
%MACRO TVALUE(N1,N2);
%LET DF=%EVAL(&N1+&N2-2);		/* stupnjevi slobode t statistike	*/

 

/*	T statistiku racunamo u dva koraka	*/

/*	Racunamo statistike koja su nam potrebne za racunanje t statistike	*/
/*	Racunamo statistike odvojeno za uzorke po replikacijama */ 

PROC MEANS DATA=GEN(WHERE=(Cc=1)) NOPRINT;
	VAR X;
	BY DELTA SAMPLEID;
	OUTPUT OUT=T_C1	MEAN=M1_X STD=STD1_X;
RUN;

PROC MEANS DATA=GEN(WHERE=(Cc=2)) NOPRINT;
	VAR X;
	BY DELTA SAMPLEID;
	OUTPUT OUT=T_C2	MEAN=M2_X STD=STD2_X;
RUN;

/*	Spajamo outpute i racunamo t statistiku testa	i usoredujemo sa 	*/
/*	kriticnom vrijednosti, za alfa=0.01. U varijablu rejecth0_n spremamo	*/
/*	podatak (oznaka 1) ukoliko odbacujemo nultu hipotezu 							*/

DATA TALL;
	SET T_C1;
	SET T_C2;
	
	T1=TINV(0.995,&DF);
	T2=TINV(0.005,&DF);/* Kriticna vrijednost	*/
		TN=(M1_X-M2_X)/(SQRT((((&N1-1)*STD1_X**2+(&N2-1)*STD2_X**2)*(&N1+&N2))/((&N1+&N2-2)*&N1*&N2)));
		REJECTH0=((TN GE T1) OR (TN LE T2)) ;	/*	usporedba	*/
RUN;
%MEND TVALUE;
/*	*	*	*	*	*	*	*/
/*	*	*	*	*	*	*	*/

%MACRO USPOREDBA();
/*	pomocu power procedure usporedujemo jakost testa	*/
/*	sa jakosti testa prema nasim podacima			*/
	PROC POWER;
		TWOSAMPLEMEANS
		TEST=DIFF
		DIST=NORMAL
		SIDES=2
		ALPHA=0.01
		MEANDIFF=-6 TO 6 BY 2
		STDDEV=&SIG_1
		GROUPNS= (&N1 &N2)
		POWER=.;
		ODS OUTPUT OUTPUT=POWER_PROC;
	RUN;
	
	DATA POWER_PROC;
		SET POWER_PROC;
		P=PERCENT/100;
		LABEL P="POWWER";
		DROP PERCENT;
	RUN;
/*	spremamo nase podatke sa podacima iz power procedure u isti data set 	*/	
	DATA POM;
		SET POWER_PROC SIMPOWER;
	RUN;	
	ods graphics off;
	goptions reset=all ftext=swiss cback=white ctext=BLACK;
	legend label=none value=(font=swiss color=black 'SNAGA TESTA ZA SIMULIRANI UZORAK' 'SNAGA TESTA POMOCU PROC POWER')
       position=(top right inside) mode=share cborder=black;;
symbol1 value="dot" interpol=line color=red;
symbol2 value="dot" interpol=line color=blue;

axis1 label=("Razlika u sredinama uzorka")
      offset=(5)
      width=3;

axis2 label=('Snaga testa')
      order=(-0.2 to 1.2 by 0.2)
      width=3;
/*	graf usporedbe power procedure i nasih podataka  	*/
	TITLE "SNAGA TESTA ZA N(0,&SIG_1)";
	title2 "duljina prvog uzorka &N1, duljina drugog uzorka &N2";
	PROC GPLOT DATA=POM; 
		plot P*DELTA POWER*MEANDIFF/ overlay
				haxis=axis1
                  		vaxis=axis2
				frame
				legend=legend1; 
	RUN;
	TITLE;
%MEND USPOREDBA;
/*	*	*	*	*	*	*	*/
/*	*	*	*	*	*	*	*/		
%macro program(NN1,NN2,SIG_1,SIG_2,SKEW_1,KURT_1);

DATA SNAGA;
	LENGTH REJECTH0 DELTA COUNT P N1 N2 8 ;
	STOP;
RUN;

%LET KK=1;
%LET N1=%SCAN(&NN1,&KK);
%LET N2=%SCAN(&NN2,&KK);

%IF &N1=&N2 %THEN %DO; %LET UZORAK="Uzorci jednakih duljina";%END;
%ELSE %DO; %LET UZORAK="Uzorci različitih duljina";%END;
%IF &SIG_1=&SIG_2 %THEN %DO; %LET VARIJANCA="Uzorci jednakih varijanci";%END;
%ELSE %DO; %LET VARIJANCA="Uzorci različitih varijanci";%END;

%DO %WHILE(&N1 NE);			/* pocetak petlje po duljini uzoraka	*/

%GENERATOR(&N1,&N2,&SIG_1,&SIG_2,&SKEW_1,&KURT_1);

%TVALUE(&N1,&N2);
/*	Sa proc freq, racunamo koliko smo pota odbacili H_0	*/

PROC FREQ DATA=TALL NOPRINT;
	TABLES REJECTH0*delta/CROSSLIST SPARSE OUT=SIMPOWER(WHERE=(REJECTH0=1));
RUN;

/*	*/
DATA SIMPOWER;
	SET SIMPOWER;
	P=COUNT/&NUMSAMPLES;
	N1=&N1;
	N2=&N2;
	LABEL P="POWER";
	DROP PERCENT;
RUN;

/*	bitne podatke cuvamo u datasetu snaga_normal za */
/*	normalnu distribuciju.				*/

DATA SNAGA;
	SET SNAGA SIMPOWER;
RUN;

%IF &SKEW_1=0 and &KURT_1=0 AND &SIG_1=&SIG_2 %then
		%USPOREDBA();





%LET KK=%EVAL(&KK+1); 
%LET N1=%SCAN(&NN1,&KK);
%LET N2=%SCAN(&NN2,&KK);
/*	prelazimo na iduci par duljine uzoraka */

%END; /*	zatvara se pocetna petlja	*/




PROC SORT DATA=SNAGA;
	BY DELTA N1;
RUN;
TITLE "Usporedba uzoraka sa Skewness = &SKEW_1  i Kurtosis = &KURT_1";
title2 "&varijanca";
title3 "&uzorak";
PROC PRINT;RUN;


PROC SGPLOT DATA=SNAGA; 
	SERIES X=DELTA Y=P /GROUP=N1 GROUPORDER=ASCENDING;
RUN;
title;
TITLE;

%MEND program;


%LET SEED1=22334;
%LET SEED2=45678;
%LET UZORAK0=10 20 40 100;
%LET UZORAK1=5 10 20 50;
%LET UZORAK2=15 30 60 150;
%LET GOPT=HBY=0; /**SUPPRESES THE BY LINE**/
%LET SIGMA0=2.24;
%LET SIGMA1=3;
%LET SIGMA2=1;
%LET NUMSAMPLES=500;
%let sk1 = -2;
%let sk2 = 0;
%let sk3 = 2;
%let kurt1 = 0;
%let kurt2 = 6;
%let kurt3 = 11;

/* DULJINE ISTE, VARIJANCA ISTA */
%program(&UZORAK0,&UZORAK0,&SIGMA0,&SIGMA0, &sk1, &kurt2);
%program(&UZORAK0,&UZORAK0,&SIGMA0,&SIGMA0, &sk1, &kurt3);
%program(&UZORAK0,&UZORAK0,&SIGMA0,&SIGMA0, &sk2, &kurt1);
%program(&UZORAK0,&UZORAK0,&SIGMA0,&SIGMA0, &sk2, &kurt2);
%program(&UZORAK0,&UZORAK0,&SIGMA0,&SIGMA0, &sk2, &kurt3);
%program(&UZORAK0,&UZORAK0,&SIGMA0,&SIGMA0, &sk3, &kurt2);
%program(&UZORAK0,&UZORAK0,&SIGMA0,&SIGMA0, &sk3, &kurt3);

/* DULJINE ISTE, VARIJANCA RAZLIČITA*/
%program(&UZORAK0,&UZORAK0,&SIGMA1,&SIGMA2, &sk1, &kurt2);
%program(&UZORAK0,&UZORAK0,&SIGMA1,&SIGMA2, &sk1, &kurt3);
%program(&UZORAK0,&UZORAK0,&SIGMA1,&SIGMA2, &sk2, &kurt1);
%program(&UZORAK0,&UZORAK0,&SIGMA1,&SIGMA2, &sk2, &kurt2);
%program(&UZORAK0,&UZORAK0,&SIGMA1,&SIGMA2, &sk2, &kurt3);
%program(&UZORAK0,&UZORAK0,&SIGMA1,&SIGMA2, &sk3, &kurt2);
%program(&UZORAK0,&UZORAK0,&SIGMA1,&SIGMA2, &sk3, &kurt3);

/* DULJINE RAZLIČITE, VARIJANCA ISTA*/
%program(&UZORAK1,&UZORAK2,&SIGMA0,&SIGMA0, &sk1, &kurt2);
%program(&UZORAK1,&UZORAK2,&SIGMA0,&SIGMA0, &sk1, &kurt3);
%program(&UZORAK1,&UZORAK2,&SIGMA0,&SIGMA0, &sk2, &kurt1);
%program(&UZORAK1,&UZORAK2,&SIGMA0,&SIGMA0, &sk2, &kurt2);
%program(&UZORAK1,&UZORAK2,&SIGMA0,&SIGMA0, &sk2, &kurt3);
%program(&UZORAK1,&UZORAK2,&SIGMA0,&SIGMA0, &sk3, &kurt2);
%program(&UZORAK1,&UZORAK2,&SIGMA0,&SIGMA0, &sk3, &kurt3);

/* DULJINE RAZLIČITE, VARIJANCA RAZLIČITA */
%program(&UZORAK1,&UZORAK2,&SIGMA1,&SIGMA2, &sk1, &kurt2);
%program(&UZORAK1,&UZORAK2,&SIGMA1,&SIGMA2, &sk1, &kurt3);
%program(&UZORAK1,&UZORAK2,&SIGMA1,&SIGMA2, &sk2, &kurt1);
%program(&UZORAK1,&UZORAK2,&SIGMA1,&SIGMA2, &sk2, &kurt2);
%program(&UZORAK1,&UZORAK2,&SIGMA1,&SIGMA2, &sk2, &kurt3);
%program(&UZORAK1,&UZORAK2,&SIGMA1,&SIGMA2, &sk3, &kurt2);
%program(&UZORAK1,&UZORAK2,&SIGMA1,&SIGMA2, &sk3, &kurt3);


PROC POWER;
	TWOSAMPLEMEANS
	TEST=DIFF
	DIST=NORMAL
	SIDES=2
	ALPHA=0.01
	MEANDIFF=1
	STDDEV=1
	ntotal=.
	POWER=0.8;
RUN;	