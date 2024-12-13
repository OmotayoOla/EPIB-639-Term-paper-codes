
/*............................Primary Analyses*................................../

/*PART 1: Cohort Definition*/

/*NOTE: Replace filepaths with appropriate filepaths*/

proc import datafile="/home/u64009092/group_assignment/prescriptions.csv" out=rx dbms=csv replace;
    getnames=yes;
run;

proc import datafile="/home/u64009092/group_assignment/event.csv" out=event dbms=csv replace;
    getnames=yes;
run;

proc import datafile="/home/u64009092/group_assignment/patient.csv" out=patient dbms=csv replace;
    getnames=yes;
run;

proc import datafile="/home/u64009092/group_assignment/covariates.csv" out=covariates dbms=csv replace;
    getnames=yes;
run;

/*Restrict rx dataset to 1995 - 2015*/

/* Identify IDs with pre-1995 prescriptions. Total = 903 patients*/

proc sql;
    create table pre1995_ids as
    select distinct id
    from rx
    where date < '1jan1995'd;

    /* Exclude these IDs from the cohort. 9097 patients remain in the cohort*/
   
proc sql;
    create table first_rx as
    select id, min(date) as first_rxdate format=date9.
    from rx
    where id not in (select id from pre1995_ids)
    group by id;
quit;

/*Apply this 1995 restriction to the other datasets: patients, events and covariates*/

/*Patients dataset : We end up with 9097 patients*/

proc sql;
	create table fil_patients as
	select a.*        				/*Select all columns from the patients table aliased as a*/
	from patient as a
	inner join first_rx as b
	on a.id = b.id;
quit;

/*Events dataset: We end up with 3312 events*/

proc sql;
    create table fil_event as
    select a.*
    from event as a
    inner join first_rx as b
    on a.id = b.id;
quit;

/*Covariates dataset: We end up with covariates for 9097 patients, as expected*/

proc sql;
    create table fil_covariates as
    select a.*
    from covariates as a
    inner join first_rx as b
    on a.id = b.id;
quit;

/*Remember you took only the first rx row to identify eligible patients. Now, restore all rows*/
/*We end up with 143303 rows*/

proc sql;
    create table fil_rx as
    select a.*
    from rx as a
    inner join first_rx as b
    on a.id = b.id;
quit;

/*Now, our 4 datasets modified to 1995 - 2015 are: fil_rx, fil_event, fil_patient, and fil_covariates*/

/*Step 2: Write codes for exclusion criteria*/

/*Exclude people with diagnosis before start of follow up. Total = 0*/

proc sql;
	create table exc_event as
	select *						/*Asterisk for consistency because it selects all columns*/
	from fil_event
	where eventdate < '1jan1995'd;
quit;

/*Exclude people with prescriptions before 1995. Total 0 */

proc sql;
	create table exc_rx1 as
	select *
	from fil_rx
	where date < '1jan1995'd;
quit;

/*Exclude people with concominant use of 2 rx generations at start of follow up. Total = 27 */

proc sql;
    create table exc_rx2 as
    select a.*
    from fil_rx as a
    inner join first_rx as b
    on a.id = b.id and a.date = b.first_rxdate
    group by a.id
    having count(distinct a.Generation) > 1; /* Multiple distinct generations */
quit;

/*To see count for number of distinct IDs. Total = 10 people*/

proc sql;
    select count(distinct id) as distinct_id_count
    from exc_rx2;
quit;

/*Exclude people with concomitant use of 2 rx from same generation at start of follow up*/

proc sql;
    create table exc_rx3 as
    select a.*
    from fil_rx as a
    inner join first_rx as b
    on a.id = b.id and a.date = b.first_rxdate
    where a.id not in (select id from exc_rx2) 
    group by a.id, a.Generation
    having count(*) > 1; /* Multiple prescriptions of the same generation */
quit;

/*To see count for number of distinct IDs. Total = 385 people*/

proc sql;
    select count(distinct id) as distinct_id_count
    from exc_rx3;
quit;

/*Exclude people less than 18, that were not excluded in previous steps. Total = 187 people*/

proc sql;
    create table exc_age as
    select *
    from fil_covariates
    where age < 18
      and id not in (select id from exc_rx2)
      and id not in (select id from exc_rx3);
quit;

/*Combine all exclusion tables. Total = 582*/

proc sql;
    create table exclusions as
    select id from exc_event
    union
    select id from exc_rx1 
    union
    select id from exc_rx2 
    union
    select id from exc_rx3
    union
    select id from exc_age;
quit;

/*Apply exclusions to each dataset*/

/*Event; remaining = 3175 events*/

data final_event;
    merge fil_event (in=a) exclusions (in=b);
    by id;
    if a and not b;
run;

/*Prescriptions (rx); remaining = 127792 rows */

data final_rx;
    merge fil_rx (in=a) exclusions (in=b);
    by id;
    if a and not b;
run;

/*Covariates; remaining = 8515 patients*/

data final_covariates;
    merge fil_covariates (in=a) exclusions (in=b);
    by id;
    if a and not b;
run;

/*Patients; remaining = 8515 patients*/

data final_patients;
    merge fil_patients (in=a) exclusions (in=b);
    by id;
    if a and not b;
run;

/*Create cohort merging patients, prescription and covariates*/

/*First, sort the 3 datasets*/

proc sort data = final_rx; by id; run;

proc sort data = final_covariates; by id; run;

proc sort data = final_patients; by id; run;

/*Then merge these 3 datasets. Total = 127792 rows*/

data final;
 merge final_rx final_covariates final_patients;
 by id;
run;


/*PART 2: Define the Exposure, Lag, Grace periods, Censoring dates etc*/

data cohort1;
 set final_rx;
 by id date;
 if first.id;		/*We select only the first id and assign the dates time 0*/
 rename date=t0;
run;

/*Identify switchers and determine the switch date*/

data switch;
	merge cohort1(in=a) final_rx(where=(date <'31dec2015'd) rename=(Generation=gen));
	by id;
	if a;
	if Generation ne gen;		/*Only the people switch make it into this dataset*/
run;

/*Create a dataset with only the first patient ID and the date of censoring*/

data switch_data;
	set switch;
	by id date;
	if first.id;			/*Keep only first id. Remember you merged with longform final_rx above*/
	keep id date;			/*Keep these 2 variables and drop everything else*/
	rename date=censor1;	
run;

/*Merge 3 datasets. So we now have t0 from 'cohort1' and censoring date (censor1) from 'switch_data' 
added to the long form (final_rx). T0 is fixed for each person while date varies (date of every refill rx)*/

data continuous_data;
	merge cohort1(in=a) final_rx switch_data;
	by id;
	if a;
	if date < censor1 or censor1=.;
run;

/*Apply rxdate + 30 + grace period (365 days)*/

data c1;
	set continuous_data;
	by id date;
	rxdur=date+30+365;				/*rxdur id prescription duration plus grace period*/
	format rxdur lagrxdur date9.;
	lagrxdur=lag(rxdur);			/*lagrxdur is the end date for the previous rx*/
	if first.id then lagrxdur=date;
	retain discon;
	if first.id then discon=0;
	if date > lagrxdur then discon=1;
	if discon=1 and id=lag(id) and discon=lag(discon) then delete; 
run;

/*Calculate the final censoring date (censor2) for each patient. 'Censor1' only identified 
switchers by capturing their first switching date. Censor two accounts for grace periods in both switchers and non-switchers.*/

data c2;
 	set c1;
 	by id date;
 	format censor2 date9.;
 	if discon = 1 then censor2 = lagrxdur + 365; /*For switchers: Added the 365 for the grace period post rx cessation*/
 	else if last.id then censor2 = rxdur + 365;	 /*For those who stop completely: Added the 365 for the grace period post rx cessation*/
 	if last.id;
run;

/*Now, define and merge the events*/

data merged_data;
    merge final_event (in=a) c2(in=b) final_patients(keep=id end); /* Add patients to get 'end' date*/
    by id;
    if b; 
    end_followup = min(eventdate,end,'31dec2015'd,censor2);
    format end date9.;
    if a and t0 < eventdate <= end_followup then event=1; 
	else event = 0;
	time = end_followup - t0 + 1;
	time_years = time/365.25;
	if time < 365.25 then event=0; /*Apply 1 year lag where we don't count events*/
run;


/*PART 3: Estimate PS and deal with confounding*/

/*First, merge the covariates dataset with the current dataset (but sort first)*/

proc sort data=merged_data; by id; run;
proc sort data=final_covariates; by id; run;

/* Step 2: Merge datasets */

data merged_data2;
    merge merged_data (in=a) final_covariates (in=b);
    by id;
    if a; 	/* Keeps only rows that exist in the main dataset*/
run;

/*Extract the year variable*/

data merged_data3;
    set merged_data2;
    t0_year = year(t0); /* Extract the year from the t0 date */
run;

/* Print the first 10 observations to verify */

proc print data=merged_data3(obs=10); 
    var t0 t0_year;
run;

/*We will create 5-year calendar time bands*/

data merged_data_bands;
    set merged_data3;
    
    /* Create calendar time bands */
    if t0_year >= 1995 and t0_year < 2000 then time_band = "1995-1999";
    else if t0_year >= 2000 and t0_year < 2005 then time_band = "2000-2004";
    else if t0_year >= 2005 and t0_year < 2010 then time_band = "2005-2009";
    else time_band = "2010-2015";

    /* Combine age and calendar time bands */
    combined_band = time_band;
run;

/*Print a few rows to confirm*/

proc print data=merged_data_bands(obs=10); 
    var time_band t0_year combined_band;
run;

proc sort data=merged_data_bands; by combined_band; run;

/*Estimate the PS*/

proc logistic data=merged_data_bands;
	
	/*Binary and categorical variables*/
    class obesity(ref="0") smoking(ref="0") male(ref="0") alcohol(ref="0") cancer(ref="0") 
          dementia(ref="0") nsaids(ref="0") aspirin(ref="0") pvd(ref="0") cad(ref="0") revasc(ref="0") 
          nitrates(ref="0") chf(ref="0") cerebrovascular(ref="0") thyroid(ref="0") statins(ref="0") 
          renal(ref="0") digoxin(ref="0") fibrates(ref="0") clopidogrel(ref="0") warfarin(ref="0") 
          opioids(ref="0") af(ref="0") copd(ref="0") asthma(ref="0") dyslipidemia(ref="0") anaemia(ref="0") 
          angina(ref="0") htn(ref="0") hypotension(ref="0") diabetes(ref="0") mi(ref="0") stroke(ref="0") / param=ref;

    /* Include quadratic and cubic terms for continuous variables */
    model Generation(event="1") = age age*age age*age*age visit visit*visit visit*visit*visit
                                  obesity smoking male alcohol cancer dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
    by combined_band;
    output out=ps_scores p=pscore;
run;

/*Print 30 rows to see propensity scores*/

proc print data=ps_scores(obs=30); run;

/*View overlap*/

proc sgplot data=ps_scores;
    histogram pscore / group=Generation transparency=0.5;
run;


/*Now, we will conduct fine stratification weighting*/

/*Create 50 strata*/

proc rank data=ps_scores out=ps_scores_strata groups=50;
    var pscore; /* The propensity score variable */
    ranks stratum; /* Creates the stratum variable (0 to 49) */
run;

/* Adjust stratum to 1-based index for easier interpretation*/

data ps_scores_strata;
    set ps_scores_strata;
    stratum = stratum + 1;
run;

/*Calculate population counts for each stratum*/

proc sql;
    create table stratum_counts as
    select stratum,
           count(*) as Ntotal,
           sum(Generation = 1) as Nexposed, /* Assuming Generation=1 is the exposed group */
           sum(Generation = 2) as Nunexposed, /* Assuming Generation=2 is the unexposed group */
           (select count(*) from ps_scores_strata where Generation = 1) as Ntotal_exposed,
           (select count(*) from ps_scores_strata where Generation = 2) as Ntotal_unexposed
    from ps_scores_strata
    group by stratum;
quit;

/*Calculate weights*/

/*First, calculate the total number of obs (sum(Ntotal)) across all strata and save it as a macro variable*/

proc sql noprint;
    select sum(Ntotal) into :Total_N from stratum_counts;
quit;

/* Check the macro variable value */

%put &=Total_N;

/*Use precomputed &Total_N to calculate the weights (but sort first)*/

proc sort data=ps_scores_strata; by stratum; run;

proc sort data=stratum_counts; by stratum; run;

/*Calculate PS Weights*/

data ps_scores_weighted;
    merge ps_scores_strata (in=a) stratum_counts;
    by stratum;
    if a then do;
        
/* Calculate weights for treated and untreated */
    if Generation = 1 then 
         Weight = (Ntotal / &Total_N) / (Nexposed / Ntotal_exposed);
    else if Generation = 2 then 
         Weight = (Ntotal / &Total_N) / (Nunexposed / Ntotal_unexposed);
    
    /* Truncate weights above 10 */
        if Weight > 10 then Weight = 10;
    end;
run;

/*Check out the first 20 weights*/

proc print data=ps_scores_weighted(obs=20);
run;

/*Check summary statistics of the weight*/

proc means data=ps_scores_weighted n mean median min max stddev;
    var Weight;
    class Generation; /* Separate summary by exposure group */
run;

/*Check histogram of the weights*/

proc sgplot data=ps_scores_weighted;
    histogram Weight / group=Generation transparency=0.5;
    xaxis label="Weights";
    yaxis label="Frequency";
run;

proc sgplot data=ps_scores_weighted;
    density Weight / group=Generation type=kernel transparency=0.5;
    keylegend / position=topright; /* Legend for clarity */
    xaxis label="Weights";
    yaxis label="Density";
    title "Weight Distribution by Group (Area Plot)";
run;

/*Check how many weights are >10*/

proc sql;
    select count(*) as WeightsAbove10
    from ps_scores_weighted
    where Weight > 10;
quit;

/*Sort and Merge the dataset*/

proc sort data=ps_scores_weighted; by id; run;
proc sort data=final_event; by id; run;

/*Merge*/

data finally;
    merge ps_scores_weighted (in=a) final_event;
    by id;
    if a;
run;

/*Finally, run the cox model*/

proc phreg data=finally;
    class Generation(ref="2") / param=ref;
    model time*event(0) = Generation / ties=efron;
    strata combined_band; /* Stratify by combined age and calendar time bands */
    weight Weight; /* Apply fine stratification weights */
    hazardratio Generation / cl=wald; /* Include confidence intervals */
run;

/*Get incidence rates, person times and CIs (these are unweighted and we didnt use them)*/

/*............ Weighted incidence rate (what we used)..........*/

proc means data=finally sum; 
var event time_years;
class Generation;
weight Weight; /* Apply weights */
output out=Inc_rate sum=event py;
run;

data Inc_rate2;													/*Create a new dataset named 'IR'*/
set Inc_rate;													/*Tell SAS to use data from IRprep*/
Inc_rate2 = event/(py/1000);									/*Calculate the incidence rate*/
LCI=quantile('chisq',0.025, event*2)/((py/100000)*2);			/*Calculate lower confidence interval*/
UCI=quantile('chisq',0.975, (event+1)*2)/((py/100000)*2);run;	/*calculate the upper confidence interval*/
run;

proc print data=Inc_rate2;
run;

/*Some diagnostics since event numbers are low*/

/*Total number of events that occured after censoring = 117 events */

proc sql;
    create table events_after_censor as
    select * from final_event as fe
    left join finally as f
    on fe.id = f.id
    where fe.eventdate > f.end_followup;
quit;

proc sql;
    select count(*) as events_after_censor from events_after_censor;
quit;

/*Total number of events that occured within the one year lag = 2884 events */

proc sql;
    create table first_year_events as
    select * from final_event as fe
    left join finally as f
    on fe.id = f.id
    where f.time < 365.25;
quit;

proc sql;
    select count(*) as excluded_first_year_events from first_year_events;
quit;



/*...................................Secondary Amalyses......................................*/

/*..................................Secondary Analyses.......................................*/


/*..............................SECONDARY ANALYSES.........................................*/

/*Split data into males and female and basically repeat eveything for each group*/

/*First, merge this dataset with covariates*/

proc sort data=merged_data; by id; run;
proc sort data=final_covariates; by id; run;

data merged_data_cov;
    merge merged_data(in=a) final_covariates(in=b);
    by id;
    if a; /* Keep only rows present in the main dataset */
run;

data merged_data_male merged_data_female;
    set merged_data_cov;
    if male = 1 then output merged_data_male;
    else if male = 0 then output merged_data_female;
run;

/*PART 3: Estimate PS and deal with confounding*/

/*Extract the year variable*/

data merged_data_new;
    set merged_data_male;
    t0_year = year(t0); /* Extract the year from the t0 date */
run;

/* Print the first 10 observations to verify */

proc print data=merged_data_new(obs=10); 
    var t0 t0_year;
run;

/*We will create 5-year calendar time bands*/

data merged_data_bands;
    set merged_data_new;
    
    /* Create calendar time bands */
    if t0_year >= 1995 and t0_year < 2000 then time_band = "1995-1999";
    else if t0_year >= 2000 and t0_year < 2005 then time_band = "2000-2004";
    else if t0_year >= 2005 and t0_year < 2010 then time_band = "2005-2009";
    else time_band = "2010-2015";

    /* Combine age and calendar time bands */
    combined_band = time_band;
run;

/*Print a few rows to confirm*/

proc print data=merged_data_bands(obs=10); 
    var time_band t0_year combined_band;
run;

proc sort data=merged_data_bands; by combined_band; run;

/*Estimate the PS*/

proc logistic data=merged_data_bands;
	
	/*Binary and categorical variables*/
    class obesity(ref="0") smoking(ref="0") alcohol(ref="0") cancer(ref="0") 
          dementia(ref="0") nsaids(ref="0") aspirin(ref="0") pvd(ref="0") cad(ref="0") revasc(ref="0") 
          nitrates(ref="0") chf(ref="0") cerebrovascular(ref="0") thyroid(ref="0") statins(ref="0") 
          renal(ref="0") digoxin(ref="0") fibrates(ref="0") clopidogrel(ref="0") warfarin(ref="0") 
          opioids(ref="0") af(ref="0") copd(ref="0") asthma(ref="0") dyslipidemia(ref="0") anaemia(ref="0") 
          angina(ref="0") htn(ref="0") hypotension(ref="0") diabetes(ref="0") mi(ref="0") stroke(ref="0") / param=ref;

    /* Include quadratic and cubic terms for continuous variables */
    model Generation(event="1") = age age*age age*age*age visit visit*visit visit*visit*visit
                                  obesity smoking alcohol cancer dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
    by combined_band;
    output out=ps_scores p=pscore;
run;

/*Print 30 rows to see propensity scores*/

proc print data=ps_scores(obs=30); run;

/*View overlap*/

proc sgplot data=ps_scores;
    histogram pscore / group=Generation transparency=0.5;
run;


/*Now, we will conduct fine stratification weighting*/

/*Create 50 strata*/

proc rank data=ps_scores out=ps_scores_strata groups=50;
    var pscore; /* The propensity score variable */
    ranks stratum; /* Creates the stratum variable (0 to 49) */
run;

/* Adjust stratum to 1-based index for easier interpretation*/

data ps_scores_strata;
    set ps_scores_strata;
    stratum = stratum + 1;
run;

/*Calculate population counts for each stratum*/

proc sql;
    create table stratum_counts as
    select stratum,
           count(*) as Ntotal,
           sum(Generation = 1) as Nexposed, /* Assuming Generation=1 is the exposed group */
           sum(Generation = 2) as Nunexposed, /* Assuming Generation=2 is the unexposed group */
           (select count(*) from ps_scores_strata where Generation = 1) as Ntotal_exposed,
           (select count(*) from ps_scores_strata where Generation = 2) as Ntotal_unexposed
    from ps_scores_strata
    group by stratum;
quit;

/*Calculate weights*/

/*First, calculate the total number of obs (sum(Ntotal)) across all strata and save it as a macro variable*/

proc sql noprint;
    select sum(Ntotal) into :Total_N from stratum_counts;
quit;

/* Check the macro variable value */

%put &=Total_N;

/*Use precomputed &Total_N to calculate the weights (but sort first)*/

proc sort data=ps_scores_strata; by stratum; run;

proc sort data=stratum_counts; by stratum; run;

/*Calculate PS Weights*/

data ps_scores_weighted;
    merge ps_scores_strata (in=a) stratum_counts;
    by stratum;
    if a then do;
        
/* Calculate weights for treated and untreated */
    if Generation = 1 then 
         Weight = (Ntotal / &Total_N) / (Nexposed / Ntotal_exposed);
    else if Generation = 2 then 
         Weight = (Ntotal / &Total_N) / (Nunexposed / Ntotal_unexposed);
    
    /* Truncate weights above 10 */
        if Weight > 10 then Weight = 10;
    end;
run;

/*Check out the first 20 weights*/

proc print data=ps_scores_weighted(obs=20);
run;

/*Check summary statistics of the weight*/

proc means data=ps_scores_weighted n mean median min max stddev;
    var Weight;
    class Generation; /* Separate summary by exposure group */
run;

/*Check histogram of the weights*/

proc sgplot data=ps_scores_weighted;
    histogram Weight / group=Generation transparency=0.5;
    xaxis label="Weights";
    yaxis label="Frequency";
run;

proc sgplot data=ps_scores_weighted;
    density Weight / group=Generation type=kernel transparency=0.5;
    keylegend / position=topright; /* Legend for clarity */
    xaxis label="Weights";
    yaxis label="Density";
    title "Weight Distribution by Group (Area Plot)";
run;

/*Check how many weights are >10*/

proc sql;
    select count(*) as WeightsAbove10
    from ps_scores_weighted
    where Weight > 10;
quit;

/*Sort and Merge the dataset*/

proc sort data=ps_scores_weighted; by id; run;
proc sort data=final_event; by id; run;

/*Merge*/

data finally;
    merge ps_scores_weighted (in=a) final_event;
    by id;
    if a;
run;

/*Finally, run the cox model*/

proc phreg data=finally;
    class Generation(ref="2") / param=ref;
    model time*event(0) = Generation / ties=efron;
    strata combined_band; /* Stratify by combined age and calendar time bands */
    weight Weight; /* Apply fine stratification weights */
    hazardratio Generation / cl=wald; /* Include confidence intervals */
run;

/*Get incidence rates, person times and CIs*/

/*Now, calculate your incidence rate*/
proc means data = finally;					/*Calls on the proc means procedure for basic descriptive statistics*/
var event time_years;						/*Specifies the variables we are intetested in to clculate the incidence rate*/
class Generation;							/*Tells SAS to do the cuaculation for each treatment group*/
weight Weight; /* Apply weights */
output out = Inc_rate sum = event py;		/*Creates a new dataset (Irate) to store the results we are about to generate*/
											/*The sum part provides names for the sum of the variables listed above (events and follow-up)*/
quit;

data Inc_rate2;													/*Create a new dataset named 'IR'*/
set Inc_rate;													/*Tell SAS to use data from IRprep*/
Inc_rate2 = event/(py/1000);									/*Calculate the incidence rate*/
LCI=quantile('chisq',0.025, event*2)/((py/100000)*2);			/*Calculate lower confidence interval*/
UCI=quantile('chisq',0.975, (event+1)*2)/((py/100000)*2);run;	/*calculate the upper confidence interval*/
run;

proc print data=Inc_rate2;
run;

/*********************************Females****************************************/

/*PART 3: Estimate PS and deal with confounding*/

/*Extract the year variable*/

data merged_data_new;
    set merged_data_female;
    t0_year = year(t0); /* Extract the year from the t0 date */
run;

/* Print the first 10 observations to verify */

proc print data=merged_data_new(obs=10); 
    var t0 t0_year;
run;

/*We will create 5-year calendar time bands*/

data merged_data_bands;
    set merged_data_new;
    
    /* Create calendar time bands */
    if t0_year >= 1995 and t0_year < 2000 then time_band = "1995-1999";
    else if t0_year >= 2000 and t0_year < 2005 then time_band = "2000-2004";
    else if t0_year >= 2005 and t0_year < 2010 then time_band = "2005-2009";
    else time_band = "2010-2015";

    /* Combine age and calendar time bands */
    combined_band = time_band;
run;

/*Print a few rows to confirm*/

proc print data=merged_data_bands(obs=10); 
    var time_band t0_year combined_band;
run;

proc sort data=merged_data_bands; by combined_band; run;

/*Estimate the PS*/

proc logistic data=merged_data_bands;
	
	/*Binary and categorical variables*/
    class obesity(ref="0") smoking(ref="0") alcohol(ref="0") cancer(ref="0") 
          dementia(ref="0") nsaids(ref="0") aspirin(ref="0") pvd(ref="0") cad(ref="0") revasc(ref="0") 
          nitrates(ref="0") chf(ref="0") cerebrovascular(ref="0") thyroid(ref="0") statins(ref="0") 
          renal(ref="0") digoxin(ref="0") fibrates(ref="0") clopidogrel(ref="0") warfarin(ref="0") 
          opioids(ref="0") af(ref="0") copd(ref="0") asthma(ref="0") dyslipidemia(ref="0") anaemia(ref="0") 
          angina(ref="0") htn(ref="0") hypotension(ref="0") diabetes(ref="0") mi(ref="0") stroke(ref="0") / param=ref;

    /* Include quadratic and cubic terms for continuous variables */
    model Generation(event="1") = age age*age age*age*age visit visit*visit visit*visit*visit
                                  obesity smoking alcohol cancer dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
    by combined_band;
    output out=ps_scores p=pscore;
run;

/*Print 30 rows to see propensity scores*/

proc print data=ps_scores(obs=30); run;

/*View overlap*/

proc sgplot data=ps_scores;
    histogram pscore / group=Generation transparency=0.5;
run;


/*Now, we will conduct fine stratification weighting*/

/*Create 50 strata*/

proc rank data=ps_scores out=ps_scores_strata groups=50;
    var pscore; /* The propensity score variable */
    ranks stratum; /* Creates the stratum variable (0 to 49) */
run;

/* Adjust stratum to 1-based index for easier interpretation*/

data ps_scores_strata;
    set ps_scores_strata;
    stratum = stratum + 1;
run;

/*Calculate population counts for each stratum*/

proc sql;
    create table stratum_counts as
    select stratum,
           count(*) as Ntotal,
           sum(Generation = 1) as Nexposed, /* Assuming Generation=1 is the exposed group */
           sum(Generation = 2) as Nunexposed, /* Assuming Generation=2 is the unexposed group */
           (select count(*) from ps_scores_strata where Generation = 1) as Ntotal_exposed,
           (select count(*) from ps_scores_strata where Generation = 2) as Ntotal_unexposed
    from ps_scores_strata
    group by stratum;
quit;

/*Calculate weights*/

/*First, calculate the total number of obs (sum(Ntotal)) across all strata and save it as a macro variable*/

proc sql noprint;
    select sum(Ntotal) into :Total_N from stratum_counts;
quit;

/* Check the macro variable value */

%put &=Total_N;

/*Use precomputed &Total_N to calculate the weights (but sort first)*/

proc sort data=ps_scores_strata; by stratum; run;

proc sort data=stratum_counts; by stratum; run;

/*Calculate PS Weights*/

data ps_scores_weighted;
    merge ps_scores_strata (in=a) stratum_counts;
    by stratum;
    if a then do;
        
/* Calculate weights for treated and untreated */
    if Generation = 1 then 
         Weight = (Ntotal / &Total_N) / (Nexposed / Ntotal_exposed);
    else if Generation = 2 then 
         Weight = (Ntotal / &Total_N) / (Nunexposed / Ntotal_unexposed);
    
    /* Truncate weights above 10 */
        if Weight > 10 then Weight = 10;
    end;
run;

/*Check out the first 20 weights*/

proc print data=ps_scores_weighted(obs=20);
run;

/*Check summary statistics of the weight*/

proc means data=ps_scores_weighted n mean median min max stddev;
    var Weight;
    class Generation; /* Separate summary by exposure group */
run;

/*Check histogram of the weights*/

proc sgplot data=ps_scores_weighted;
    histogram Weight / group=Generation transparency=0.5;
    xaxis label="Weights";
    yaxis label="Frequency";
run;

proc sgplot data=ps_scores_weighted;
    density Weight / group=Generation type=kernel transparency=0.5;
    keylegend / position=topright; /* Legend for clarity */
    xaxis label="Weights";
    yaxis label="Density";
    title "Weight Distribution by Group (Area Plot)";
run;

/*Check how many weights are >10*/

proc sql;
    select count(*) as WeightsAbove10
    from ps_scores_weighted
    where Weight > 10;
quit;

/*Sort and Merge the dataset*/

proc sort data=ps_scores_weighted; by id; run;
proc sort data=final_event; by id; run;

/*Merge*/

data finally;
    merge ps_scores_weighted (in=a) final_event;
    by id;
    if a;
run;

/*Finally, run the cox model*/

proc phreg data=finally;
    class Generation(ref="2") / param=ref;
    model time*event(0) = Generation / ties=efron;
    strata combined_band; /* Stratify by combined age and calendar time bands */
    weight Weight; /* Apply fine stratification weights */
    hazardratio Generation / cl=wald; /* Include confidence intervals */
run;

/*Get incidence rates, person times and CIs*/

/*Now, calculate your incidence rate*/
proc means data = finally;					/*Calls on the proc means procedure for basic descriptive statistics*/
var event time_years;						/*Specifies the variables we are intetested in to clculate the incidence rate*/
class Generation;							/*Tells SAS to do the cuaculation for each treatment group*/
weight Weight; /* Apply weights */
output out = Inc_rate sum = event py;		/*Creates a new dataset (Irate) to store the results we are about to generate*/
											/*The sum part provides names for the sum of the variables listed above (events and follow-up)*/
quit;

data Inc_rate2;													/*Create a new dataset named 'IR'*/
set Inc_rate;													/*Tell SAS to use data from IRprep*/
Inc_rate2 = event/(py/1000);									/*Calculate the incidence rate*/
LCI=quantile('chisq',0.025, event*2)/((py/100000)*2);			/*Calculate lower confidence interval*/
UCI=quantile('chisq',0.975, (event+1)*2)/((py/100000)*2);run;	/*calculate the upper confidence interval*/
run;

proc print data=Inc_rate2;
run;

/*......................................Cancer...........................................*/


/*Split data into males and female and basically repeat eveything for each group*/

data merged_data_cancer merged_data_nocancer;
    set merged_data_cov;
    if cancer = 1 then output merged_data_cancer;
    else if cancer = 0 then output merged_data_nocancer;
run;

/*PART 3: Estimate PS and deal with confounding*/

/*..........................................Cancer.............................................*/

/*Extract the year variable*/

data merged_data_new;
    set merged_data_cancer;
    t0_year = year(t0); /* Extract the year from the t0 date */
run;

/* Print the first 10 observations to verify */

proc print data=merged_data_new(obs=10); 
    var t0 t0_year;
run;

/*We will create 5-year calendar time bands*/

data merged_data_bands;
    set merged_data_new;
    
    /* Create calendar time bands */
    if t0_year >= 1995 and t0_year < 2000 then time_band = "1995-1999";
    else if t0_year >= 2000 and t0_year < 2005 then time_band = "2000-2004";
    else if t0_year >= 2005 and t0_year < 2010 then time_band = "2005-2009";
    else time_band = "2010-2015";

    /* Combine age and calendar time bands */
    combined_band = time_band;
run;

/*Print a few rows to confirm*/

proc print data=merged_data_bands(obs=10); 
    var time_band t0_year combined_band;
run;

proc sort data=merged_data_bands; by combined_band; run;

/*Estimate the PS*/

proc logistic data=merged_data_bands;
	
	/*Binary and categorical variables*/
    class obesity(ref="0") smoking(ref="0") alcohol(ref="0") male(ref="0") 
          dementia(ref="0") nsaids(ref="0") aspirin(ref="0") pvd(ref="0") cad(ref="0") revasc(ref="0") 
          nitrates(ref="0") chf(ref="0") cerebrovascular(ref="0") thyroid(ref="0") statins(ref="0") 
          renal(ref="0") digoxin(ref="0") fibrates(ref="0") clopidogrel(ref="0") warfarin(ref="0") 
          opioids(ref="0") af(ref="0") copd(ref="0") asthma(ref="0") dyslipidemia(ref="0") anaemia(ref="0") 
          angina(ref="0") htn(ref="0") hypotension(ref="0") diabetes(ref="0") mi(ref="0") stroke(ref="0") / param=ref;

    /* Include quadratic and cubic terms for continuous variables */
    model Generation(event="1") = age age*age age*age*age visit visit*visit visit*visit*visit
                                  obesity smoking alcohol male dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
    by combined_band;
    output out=ps_scores p=pscore;
run;

/*Print 30 rows to see propensity scores*/

proc print data=ps_scores(obs=30); run;

/*View overlap*/

proc sgplot data=ps_scores;
    histogram pscore / group=Generation transparency=0.5;
run;


/*Now, we will conduct fine stratification weighting*/

/*Create 50 strata*/

proc rank data=ps_scores out=ps_scores_strata groups=50;
    var pscore; /* The propensity score variable */
    ranks stratum; /* Creates the stratum variable (0 to 49) */
run;

/* Adjust stratum to 1-based index for easier interpretation*/

data ps_scores_strata;
    set ps_scores_strata;
    stratum = stratum + 1;
run;

/*Calculate population counts for each stratum*/

proc sql;
    create table stratum_counts as
    select stratum,
           count(*) as Ntotal,
           sum(Generation = 1) as Nexposed, /* Assuming Generation=1 is the exposed group */
           sum(Generation = 2) as Nunexposed, /* Assuming Generation=2 is the unexposed group */
           (select count(*) from ps_scores_strata where Generation = 1) as Ntotal_exposed,
           (select count(*) from ps_scores_strata where Generation = 2) as Ntotal_unexposed
    from ps_scores_strata
    group by stratum;
quit;

/*Calculate weights*/

/*First, calculate the total number of obs (sum(Ntotal)) across all strata and save it as a macro variable*/

proc sql noprint;
    select sum(Ntotal) into :Total_N from stratum_counts;
quit;

/* Check the macro variable value */

%put &=Total_N;

/*Use precomputed &Total_N to calculate the weights (but sort first)*/

proc sort data=ps_scores_strata; by stratum; run;

proc sort data=stratum_counts; by stratum; run;

/*Calculate PS Weights*/

data ps_scores_weighted;
    merge ps_scores_strata (in=a) stratum_counts;
    by stratum;
    if a then do;
        
/* Calculate weights for treated and untreated */
    if Generation = 1 then 
         Weight = (Ntotal / &Total_N) / (Nexposed / Ntotal_exposed);
    else if Generation = 2 then 
         Weight = (Ntotal / &Total_N) / (Nunexposed / Ntotal_unexposed);
    
    /* Truncate weights above 10 */
        if Weight > 10 then Weight = 10;
    end;
run;

/*Check out the first 20 weights*/

proc print data=ps_scores_weighted(obs=20);
run;

/*Check summary statistics of the weight*/

proc means data=ps_scores_weighted n mean median min max stddev;
    var Weight;
    class Generation; /* Separate summary by exposure group */
run;

/*Check histogram of the weights*/

proc sgplot data=ps_scores_weighted;
    histogram Weight / group=Generation transparency=0.5;
    xaxis label="Weights";
    yaxis label="Frequency";
run;

proc sgplot data=ps_scores_weighted;
    density Weight / group=Generation type=kernel transparency=0.5;
    keylegend / position=topright; /* Legend for clarity */
    xaxis label="Weights";
    yaxis label="Density";
    title "Weight Distribution by Group (Area Plot)";
run;

/*Check how many weights are >10*/

proc sql;
    select count(*) as WeightsAbove10
    from ps_scores_weighted
    where Weight > 10;
quit;

/*Sort and Merge the dataset*/

proc sort data=ps_scores_weighted; by id; run;
proc sort data=final_event; by id; run;

/*Merge*/

data finally;
    merge ps_scores_weighted (in=a) final_event;
    by id;
    if a;
run;

/*Finally, run the cox model*/

proc phreg data=finally;
    class Generation(ref="2") / param=ref;
    model time*event(0) = Generation / ties=efron;
    strata combined_band; /* Stratify by combined age and calendar time bands */
    weight Weight; /* Apply fine stratification weights */
    hazardratio Generation / cl=wald; /* Include confidence intervals */
run;

/*Get incidence rates, person times and CIs*/

/*Now, calculate your incidence rate*/
proc means data = finally;					/*Calls on the proc means procedure for basic descriptive statistics*/
var event time_years;						/*Specifies the variables we are intetested in to clculate the incidence rate*/
class Generation;							/*Tells SAS to do the cuaculation for each treatment group*/
weight Weight; /* Apply weights */
output out = Inc_rate sum = event py;		/*Creates a new dataset (Irate) to store the results we are about to generate*/
											/*The sum part provides names for the sum of the variables listed above (events and follow-up)*/
quit;

data Inc_rate2;													/*Create a new dataset named 'IR'*/
set Inc_rate;													/*Tell SAS to use data from IRprep*/
Inc_rate2 = event/(py/1000);									/*Calculate the incidence rate*/
LCI=quantile('chisq',0.025, event*2)/((py/100000)*2);			/*Calculate lower confidence interval*/
UCI=quantile('chisq',0.975, (event+1)*2)/((py/100000)*2);run;	/*calculate the upper confidence interval*/
run;

proc print data=Inc_rate2;
run;

/*........................................No cancer............................................*

/*Split data into males and female and basically repeat eveything for each group*/

data merged_data_cancer merged_data_nocancer;
    set merged_data_cov;
    if cancer = 1 then output merged_data_cancer;
    else if cancer = 0 then output merged_data_nocancer;
run;

/*PART 3: Estimate PS and deal with confounding*/

/*Cancer*/

/*Extract the year variable*/

data merged_data_new;
    set merged_data_nocancer;
    t0_year = year(t0); /* Extract the year from the t0 date */
run;

/* Print the first 10 observations to verify */

proc print data=merged_data_new(obs=10); 
    var t0 t0_year;
run;

/*We will create 5-year calendar time bands*/

data merged_data_bands;
    set merged_data_new;
    
    /* Create calendar time bands */
    if t0_year >= 1995 and t0_year < 2000 then time_band = "1995-1999";
    else if t0_year >= 2000 and t0_year < 2005 then time_band = "2000-2004";
    else if t0_year >= 2005 and t0_year < 2010 then time_band = "2005-2009";
    else time_band = "2010-2015";

    /* Combine age and calendar time bands */
    combined_band = time_band;
run;

/*Print a few rows to confirm*/

proc print data=merged_data_bands(obs=10); 
    var time_band t0_year combined_band;
run;

proc sort data=merged_data_bands; by combined_band; run;

/*Estimate the PS*/

proc logistic data=merged_data_bands;
	
	/*Binary and categorical variables*/
    class obesity(ref="0") smoking(ref="0") alcohol(ref="0") male(ref="0") 
          dementia(ref="0") nsaids(ref="0") aspirin(ref="0") pvd(ref="0") cad(ref="0") revasc(ref="0") 
          nitrates(ref="0") chf(ref="0") cerebrovascular(ref="0") thyroid(ref="0") statins(ref="0") 
          renal(ref="0") digoxin(ref="0") fibrates(ref="0") clopidogrel(ref="0") warfarin(ref="0") 
          opioids(ref="0") af(ref="0") copd(ref="0") asthma(ref="0") dyslipidemia(ref="0") anaemia(ref="0") 
          angina(ref="0") htn(ref="0") hypotension(ref="0") diabetes(ref="0") mi(ref="0") stroke(ref="0") / param=ref;

    /* Include quadratic and cubic terms for continuous variables */
    model Generation(event="1") = age age*age age*age*age visit visit*visit visit*visit*visit
                                  obesity smoking alcohol male dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
    by combined_band;
    output out=ps_scores p=pscore;
run;

/*Print 30 rows to see propensity scores*/

proc print data=ps_scores(obs=30); run;

/*View overlap*/

proc sgplot data=ps_scores;
    histogram pscore / group=Generation transparency=0.5;
run;


/*Now, we will conduct fine stratification weighting*/

/*Create 50 strata*/

proc rank data=ps_scores out=ps_scores_strata groups=50;
    var pscore; /* The propensity score variable */
    ranks stratum; /* Creates the stratum variable (0 to 49) */
run;

/* Adjust stratum to 1-based index for easier interpretation*/

data ps_scores_strata;
    set ps_scores_strata;
    stratum = stratum + 1;
run;

/*Calculate population counts for each stratum*/

proc sql;
    create table stratum_counts as
    select stratum,
           count(*) as Ntotal,
           sum(Generation = 1) as Nexposed, /* Assuming Generation=1 is the exposed group */
           sum(Generation = 2) as Nunexposed, /* Assuming Generation=2 is the unexposed group */
           (select count(*) from ps_scores_strata where Generation = 1) as Ntotal_exposed,
           (select count(*) from ps_scores_strata where Generation = 2) as Ntotal_unexposed
    from ps_scores_strata
    group by stratum;
quit;

/*Calculate weights*/

/*First, calculate the total number of obs (sum(Ntotal)) across all strata and save it as a macro variable*/

proc sql noprint;
    select sum(Ntotal) into :Total_N from stratum_counts;
quit;

/* Check the macro variable value */

%put &=Total_N;

/*Use precomputed &Total_N to calculate the weights (but sort first)*/

proc sort data=ps_scores_strata; by stratum; run;

proc sort data=stratum_counts; by stratum; run;

/*Calculate PS Weights*/

data ps_scores_weighted;
    merge ps_scores_strata (in=a) stratum_counts;
    by stratum;
    if a then do;
        
/* Calculate weights for treated and untreated */
    if Generation = 1 then 
         Weight = (Ntotal / &Total_N) / (Nexposed / Ntotal_exposed);
    else if Generation = 2 then 
         Weight = (Ntotal / &Total_N) / (Nunexposed / Ntotal_unexposed);
    
    /* Truncate weights above 10 */
        if Weight > 10 then Weight = 10;
    end;
run;

/*Check out the first 20 weights*/

proc print data=ps_scores_weighted(obs=20);
run;

/*Check summary statistics of the weight*/

proc means data=ps_scores_weighted n mean median min max stddev;
    var Weight;
    class Generation; /* Separate summary by exposure group */
run;

/*Check histogram of the weights*/

proc sgplot data=ps_scores_weighted;
    histogram Weight / group=Generation transparency=0.5;
    xaxis label="Weights";
    yaxis label="Frequency";
run;

proc sgplot data=ps_scores_weighted;
    density Weight / group=Generation type=kernel transparency=0.5;
    keylegend / position=topright; /* Legend for clarity */
    xaxis label="Weights";
    yaxis label="Density";
    title "Weight Distribution by Group (Area Plot)";
run;

/*Check how many weights are >10*/

proc sql;
    select count(*) as WeightsAbove10
    from ps_scores_weighted
    where Weight > 10;
quit;

/*Sort and Merge the dataset*/

proc sort data=ps_scores_weighted; by id; run;
proc sort data=final_event; by id; run;

/*Merge*/

data finally;
    merge ps_scores_weighted (in=a) final_event;
    by id;
    if a;
run;

/*Finally, run the cox model*/

proc phreg data=finally;
    class Generation(ref="2") / param=ref;
    model time*event(0) = Generation / ties=efron;
    strata combined_band; /* Stratify by combined age and calendar time bands */
    weight Weight; /* Apply fine stratification weights */
    hazardratio Generation / cl=wald; /* Include confidence intervals */
run;

/*Get incidence rates, person times and CIs*/

/*Now, calculate your incidence rate*/
proc means data = finally;					/*Calls on the proc means procedure for basic descriptive statistics*/
var event time_years;						/*Specifies the variables we are intetested in to clculate the incidence rate*/
class Generation;							/*Tells SAS to do the cuaculation for each treatment group*/
weight Weight; 								/* Apply weights */
output out = Inc_rate sum = event py;		/*Creates a new dataset (Irate) to store the results we are about to generate*/
											/*The sum part provides names for the sum of the variables listed above (events and follow-up)*/
quit;

data Inc_rate2;													/*Create a new dataset named 'IR'*/
set Inc_rate;													/*Tell SAS to use data from IRprep*/
Inc_rate2 = event/(py/1000);									/*Calculate the incidence rate*/
LCI=quantile('chisq',0.025, event*2)/((py/100000)*2);			/*Calculate lower confidence interval*/
UCI=quantile('chisq',0.975, (event+1)*2)/((py/100000)*2);run;	/*calculate the upper confidence interval*/
run;

proc print data=Inc_rate2;
run;


/*...................................Age.........................................*

/*Split data into males and female and basically repeat eveything for each group*/

data merged_data_age1 merged_data_age2;
    set merged_data_cov;
    if 18 <= age <= 69 then output merged_data_age1; /* Age 50-69 */
    else if age >= 70 then output merged_data_age2; /* Age 70+ */
run;

/*PART 3: Estimate PS and deal with confounding*/

/*60-69*/

/*Extract the year variable*/

data merged_data_new;
    set merged_data_age1;
    t0_year = year(t0); /* Extract the year from the t0 date */
run;

/* Print the first 10 observations to verify */

proc print data=merged_data_new(obs=10); 
    var t0 t0_year;
run;

/*We will create 5-year calendar time bands*/

data merged_data_bands;
    set merged_data_new;
    
    /* Create calendar time bands */
    if t0_year >= 1995 and t0_year < 2000 then time_band = "1995-1999";
    else if t0_year >= 2000 and t0_year < 2005 then time_band = "2000-2004";
    else if t0_year >= 2005 and t0_year < 2010 then time_band = "2005-2009";
    else time_band = "2010-2015";

    /* Combine age and calendar time bands */
    combined_band = time_band;
run;

/*Print a few rows to confirm*/

proc print data=merged_data_bands(obs=10); 
    var time_band t0_year combined_band;
run;

proc sort data=merged_data_bands; by combined_band; run;

/*Estimate the PS*/

proc logistic data=merged_data_bands;
	
	/*Binary and categorical variables*/
    class obesity(ref="0") smoking(ref="0") alcohol(ref="0") male(ref="0") cancer(ref="0") 
          dementia(ref="0") nsaids(ref="0") aspirin(ref="0") pvd(ref="0") cad(ref="0") revasc(ref="0") 
          nitrates(ref="0") chf(ref="0") cerebrovascular(ref="0") thyroid(ref="0") statins(ref="0") 
          renal(ref="0") digoxin(ref="0") fibrates(ref="0") clopidogrel(ref="0") warfarin(ref="0") 
          opioids(ref="0") af(ref="0") copd(ref="0") asthma(ref="0") dyslipidemia(ref="0") anaemia(ref="0") 
          angina(ref="0") htn(ref="0") hypotension(ref="0") diabetes(ref="0") mi(ref="0") stroke(ref="0") / param=ref;

    /* Include quadratic and cubic terms for continuous variables */
    model Generation(event="1") = visit visit*visit visit*visit*visit
                                  obesity smoking alcohol male dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
    by combined_band;
    output out=ps_scores p=pscore;
run;

/*Print 30 rows to see propensity scores*/

proc print data=ps_scores(obs=30); run;

/*View overlap*/

proc sgplot data=ps_scores;
    histogram pscore / group=Generation transparency=0.5;
run;


/*Now, we will conduct fine stratification weighting*/

/*Create 50 strata*/

proc rank data=ps_scores out=ps_scores_strata groups=50;
    var pscore; /* The propensity score variable */
    ranks stratum; /* Creates the stratum variable (0 to 49) */
run;

/* Adjust stratum to 1-based index for easier interpretation*/

data ps_scores_strata;
    set ps_scores_strata;
    stratum = stratum + 1;
run;

/*Calculate population counts for each stratum*/

proc sql;
    create table stratum_counts as
    select stratum,
           count(*) as Ntotal,
           sum(Generation = 1) as Nexposed, /* Assuming Generation=1 is the exposed group */
           sum(Generation = 2) as Nunexposed, /* Assuming Generation=2 is the unexposed group */
           (select count(*) from ps_scores_strata where Generation = 1) as Ntotal_exposed,
           (select count(*) from ps_scores_strata where Generation = 2) as Ntotal_unexposed
    from ps_scores_strata
    group by stratum;
quit;

/*Calculate weights*/

/*First, calculate the total number of obs (sum(Ntotal)) across all strata and save it as a macro variable*/

proc sql noprint;
    select sum(Ntotal) into :Total_N from stratum_counts;
quit;

/* Check the macro variable value */

%put &=Total_N;

/*Use precomputed &Total_N to calculate the weights (but sort first)*/

proc sort data=ps_scores_strata; by stratum; run;

proc sort data=stratum_counts; by stratum; run;

/*Calculate PS Weights*/

data ps_scores_weighted;
    merge ps_scores_strata (in=a) stratum_counts;
    by stratum;
    if a then do;
        
/* Calculate weights for treated and untreated */
    if Generation = 1 then 
         Weight = (Ntotal / &Total_N) / (Nexposed / Ntotal_exposed);
    else if Generation = 2 then 
         Weight = (Ntotal / &Total_N) / (Nunexposed / Ntotal_unexposed);
    
    /* Truncate weights above 10 */
        if Weight > 10 then Weight = 10;
    end;
run;

/*Check out the first 20 weights*/

proc print data=ps_scores_weighted(obs=20);
run;

/*Check summary statistics of the weight*/

proc means data=ps_scores_weighted n mean median min max stddev;
    var Weight;
    class Generation; /* Separate summary by exposure group */
run;

/*Check histogram of the weights*/

proc sgplot data=ps_scores_weighted;
    histogram Weight / group=Generation transparency=0.5;
    xaxis label="Weights";
    yaxis label="Frequency";
run;

proc sgplot data=ps_scores_weighted;
    density Weight / group=Generation type=kernel transparency=0.5;
    keylegend / position=topright; /* Legend for clarity */
    xaxis label="Weights";
    yaxis label="Density";
    title "Weight Distribution by Group (Area Plot)";
run;

/*Check how many weights are >10*/

proc sql;
    select count(*) as WeightsAbove10
    from ps_scores_weighted
    where Weight > 10;
quit;

/*Sort and Merge the dataset*/

proc sort data=ps_scores_weighted; by id; run;
proc sort data=final_event; by id; run;

/*Merge*/

data finally;
    merge ps_scores_weighted (in=a) final_event;
    by id;
    if a;
run;

/*Finally, run the cox model*/

proc phreg data=finally;
    class Generation(ref="2") / param=ref;
    model time*event(0) = Generation / ties=efron;
    strata combined_band; /* Stratify by combined age and calendar time bands */
    weight Weight; /* Apply fine stratification weights */
    hazardratio Generation / cl=wald; /* Include confidence intervals */
run;

/*Get incidence rates, person times and CIs*/

/*Now, calculate your incidence rate*/
proc means data = finally;					/*Calls on the proc means procedure for basic descriptive statistics*/
var event time_years;						/*Specifies the variables we are intetested in to clculate the incidence rate*/
class Generation;							/*Tells SAS to do the cuaculation for each treatment group*/
weight Weight; /* Apply weights */
output out = Inc_rate sum = event py;		/*Creates a new dataset (Irate) to store the results we are about to generate*/
											/*The sum part provides names for the sum of the variables listed above (events and follow-up)*/
quit;

data Inc_rate2;													/*Create a new dataset named 'IR'*/
set Inc_rate;													/*Tell SAS to use data from IRprep*/
Inc_rate2 = event/(py/1000);									/*Calculate the incidence rate*/
LCI=quantile('chisq',0.025, event*2)/((py/100000)*2);			/*Calculate lower confidence interval*/
UCI=quantile('chisq',0.975, (event+1)*2)/((py/100000)*2);run;	/*calculate the upper confidence interval*/
run;

proc print data=Inc_rate2;
run;


/************************************Age 70+****************************************/

/*Extract the year variable*/

data merged_data_new;
    set merged_data_age2;
    t0_year = year(t0); /* Extract the year from the t0 date */
run;

/* Print the first 10 observations to verify */

proc print data=merged_data_new(obs=10); 
    var t0 t0_year;
run;

/*We will create 5-year calendar time bands*/

data merged_data_bands;
    set merged_data_new;
    
    /* Create calendar time bands */
    if t0_year >= 1995 and t0_year < 2000 then time_band = "1995-1999";
    else if t0_year >= 2000 and t0_year < 2005 then time_band = "2000-2004";
    else if t0_year >= 2005 and t0_year < 2010 then time_band = "2005-2009";
    else time_band = "2010-2015";

    /* Combine age and calendar time bands */
    combined_band = time_band;
run;

/*Print a few rows to confirm*/

proc print data=merged_data_bands(obs=10); 
    var time_band t0_year combined_band;
run;

proc sort data=merged_data_bands; by combined_band; run;

/*Estimate the PS*/

proc logistic data=merged_data_bands;
	
	/*Binary and categorical variables*/
    class obesity(ref="0") smoking(ref="0") alcohol(ref="0") male(ref="0") cancer(ref="0")
          dementia(ref="0") nsaids(ref="0") aspirin(ref="0") pvd(ref="0") cad(ref="0") revasc(ref="0") 
          nitrates(ref="0") chf(ref="0") cerebrovascular(ref="0") thyroid(ref="0") statins(ref="0") 
          renal(ref="0") digoxin(ref="0") fibrates(ref="0") clopidogrel(ref="0") warfarin(ref="0") 
          opioids(ref="0") af(ref="0") copd(ref="0") asthma(ref="0") dyslipidemia(ref="0") anaemia(ref="0") 
          angina(ref="0") htn(ref="0") hypotension(ref="0") diabetes(ref="0") mi(ref="0") stroke(ref="0") / param=ref;

    /* Include quadratic and cubic terms for continuous variables */
    model Generation(event="1") = visit visit*visit visit*visit*visit
                                  obesity smoking alcohol male dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
    by combined_band;
    output out=ps_scores p=pscore;
run;

/*Print 30 rows to see propensity scores*/

proc print data=ps_scores(obs=30); run;

/*View overlap*/

proc sgplot data=ps_scores;
    histogram pscore / group=Generation transparency=0.5;
run;


/*Now, we will conduct fine stratification weighting*/

/*Create 50 strata*/

proc rank data=ps_scores out=ps_scores_strata groups=50;
    var pscore; /* The propensity score variable */
    ranks stratum; /* Creates the stratum variable (0 to 49) */
run;

/* Adjust stratum to 1-based index for easier interpretation*/

data ps_scores_strata;
    set ps_scores_strata;
    stratum = stratum + 1;
run;

/*Calculate population counts for each stratum*/

proc sql;
    create table stratum_counts as
    select stratum,
           count(*) as Ntotal,
           sum(Generation = 1) as Nexposed, /* Assuming Generation=1 is the exposed group */
           sum(Generation = 2) as Nunexposed, /* Assuming Generation=2 is the unexposed group */
           (select count(*) from ps_scores_strata where Generation = 1) as Ntotal_exposed,
           (select count(*) from ps_scores_strata where Generation = 2) as Ntotal_unexposed
    from ps_scores_strata
    group by stratum;
quit;

/*Calculate weights*/

/*First, calculate the total number of obs (sum(Ntotal)) across all strata and save it as a macro variable*/

proc sql noprint;
    select sum(Ntotal) into :Total_N from stratum_counts;
quit;

/* Check the macro variable value */

%put &=Total_N;

/*Use precomputed &Total_N to calculate the weights (but sort first)*/

proc sort data=ps_scores_strata; by stratum; run;

proc sort data=stratum_counts; by stratum; run;

/*Calculate PS Weights*/

data ps_scores_weighted;
    merge ps_scores_strata (in=a) stratum_counts;
    by stratum;
    if a then do;
        
/* Calculate weights for treated and untreated */
    if Generation = 1 then 
         Weight = (Ntotal / &Total_N) / (Nexposed / Ntotal_exposed);
    else if Generation = 2 then 
         Weight = (Ntotal / &Total_N) / (Nunexposed / Ntotal_unexposed);
    
    /* Truncate weights above 10 */
        if Weight > 10 then Weight = 10;
    end;
run;

/*Check out the first 20 weights*/

proc print data=ps_scores_weighted(obs=20);
run;

/*Check summary statistics of the weight*/

proc means data=ps_scores_weighted n mean median min max stddev;
    var Weight;
    class Generation; /* Separate summary by exposure group */
run;

/*Check histogram of the weights*/

proc sgplot data=ps_scores_weighted;
    histogram Weight / group=Generation transparency=0.5;
    xaxis label="Weights";
    yaxis label="Frequency";
run;

proc sgplot data=ps_scores_weighted;
    density Weight / group=Generation type=kernel transparency=0.5;
    keylegend / position=topright; /* Legend for clarity */
    xaxis label="Weights";
    yaxis label="Density";
    title "Weight Distribution by Group (Area Plot)";
run;

/*Check how many weights are >10*/

proc sql;
    select count(*) as WeightsAbove10
    from ps_scores_weighted
    where Weight > 10;
quit;

/*Sort and Merge the dataset*/

proc sort data=ps_scores_weighted; by id; run;
proc sort data=final_event; by id; run;

/*Merge*/

data finally;
    merge ps_scores_weighted (in=a) final_event;
    by id;
    if a;
run;

/*Finally, run the cox model*/

proc phreg data=finally;
    class Generation(ref="2") / param=ref;
    model time*event(0) = Generation / ties=efron;
    strata combined_band; /* Stratify by combined age and calendar time bands */
    weight Weight; /* Apply fine stratification weights */
    hazardratio Generation / cl=wald; /* Include confidence intervals */
run;

/*Get incidence rates, person times and CIs*/

/*Now, calculate your incidence rate*/
proc means data = finally;					/*Calls on the proc means procedure for basic descriptive statistics*/
var event time_years;						/*Specifies the variables we are intetested in to clculate the incidence rate*/
class Generation;							/*Tells SAS to do the cuaculation for each treatment group*/
weight Weight; /* Apply weights */
output out = Inc_rate sum = event py;		/*Creates a new dataset (Irate) to store the results we are about to generate*/
											/*The sum part provides names for the sum of the variables listed above (events and follow-up)*/
quit;

data Inc_rate2;													/*Create a new dataset named 'IR'*/
set Inc_rate;													/*Tell SAS to use data from IRprep*/
Inc_rate2 = event/(py/1000);									/*Calculate the incidence rate*/
LCI=quantile('chisq',0.025, event*2)/((py/100000)*2);			/*Calculate lower confidence interval*/
UCI=quantile('chisq',0.975, (event+1)*2)/((py/100000)*2);run;	/*calculate the upper confidence interval*/
run;

proc print data=Inc_rate2;
run;



/**********************************Smoking****************************************/

/*Split data into males and female and basically repeat eveything for each group*/

data merged_data_unknown merged_data_never merged_data_past merged_data_current;
    set merged_data_cov;
    if smoking = 0 then output merged_data_unknown;
    else if smoking = 1 then output merged_data_never;
    else if smoking = 2 then output merged_data_past;
    else if smoking = 3 then output merged_data_current;
run;

/*PART 3: Estimate PS and deal with confounding*/

/*Never*/

/*Extract the year variable*/

data merged_data_new;
    set merged_data_never;
    t0_year = year(t0); /* Extract the year from the t0 date */
run;

/* Print the first 10 observations to verify */

proc print data=merged_data_new(obs=10); 
    var t0 t0_year;
run;

/*We will create 5-year calendar time bands*/

data merged_data_bands;
    set merged_data_new;
    
    /* Create calendar time bands */
    if t0_year >= 1995 and t0_year < 2000 then time_band = "1995-1999";
    else if t0_year >= 2000 and t0_year < 2005 then time_band = "2000-2004";
    else if t0_year >= 2005 and t0_year < 2010 then time_band = "2005-2009";
    else time_band = "2010-2015";

    /* Combine age and calendar time bands */
    combined_band = time_band;
run;

/*Print a few rows to confirm*/

proc print data=merged_data_bands(obs=10); 
    var time_band t0_year combined_band;
run;

proc sort data=merged_data_bands; by combined_band; run;

/*Estimate the PS*/

proc logistic data=merged_data_bands;
	
	/*Binary and categorical variables*/
    class obesity(ref="0") cancer(ref="0") alcohol(ref="0") male(ref="0") 
          dementia(ref="0") nsaids(ref="0") aspirin(ref="0") pvd(ref="0") cad(ref="0") revasc(ref="0") 
          nitrates(ref="0") chf(ref="0") cerebrovascular(ref="0") thyroid(ref="0") statins(ref="0") 
          renal(ref="0") digoxin(ref="0") fibrates(ref="0") clopidogrel(ref="0") warfarin(ref="0") 
          opioids(ref="0") af(ref="0") copd(ref="0") asthma(ref="0") dyslipidemia(ref="0") anaemia(ref="0") 
          angina(ref="0") htn(ref="0") hypotension(ref="0") diabetes(ref="0") mi(ref="0") stroke(ref="0") / param=ref;

    /* Include quadratic and cubic terms for continuous variables */
    model Generation(event="1") = age age*age age*age*age visit visit*visit visit*visit*visit
                                  obesity smoking alcohol male dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
    by combined_band;
    output out=ps_scores p=pscore;
run;

/*Print 30 rows to see propensity scores*/

proc print data=ps_scores(obs=30); run;

/*View overlap*/

proc sgplot data=ps_scores;
    histogram pscore / group=Generation transparency=0.5;
run;


/*Now, we will conduct fine stratification weighting*/

/*Create 50 strata*/

proc rank data=ps_scores out=ps_scores_strata groups=50;
    var pscore; /* The propensity score variable */
    ranks stratum; /* Creates the stratum variable (0 to 49) */
run;

/* Adjust stratum to 1-based index for easier interpretation*/

data ps_scores_strata;
    set ps_scores_strata;
    stratum = stratum + 1;
run;

/*Calculate population counts for each stratum*/

proc sql;
    create table stratum_counts as
    select stratum,
           count(*) as Ntotal,
           sum(Generation = 1) as Nexposed, /* Assuming Generation=1 is the exposed group */
           sum(Generation = 2) as Nunexposed, /* Assuming Generation=2 is the unexposed group */
           (select count(*) from ps_scores_strata where Generation = 1) as Ntotal_exposed,
           (select count(*) from ps_scores_strata where Generation = 2) as Ntotal_unexposed
    from ps_scores_strata
    group by stratum;
quit;

/*Calculate weights*/

/*First, calculate the total number of obs (sum(Ntotal)) across all strata and save it as a macro variable*/

proc sql noprint;
    select sum(Ntotal) into :Total_N from stratum_counts;
quit;

/* Check the macro variable value */

%put &=Total_N;

/*Use precomputed &Total_N to calculate the weights (but sort first)*/

proc sort data=ps_scores_strata; by stratum; run;

proc sort data=stratum_counts; by stratum; run;

/*Calculate PS Weights*/

data ps_scores_weighted;
    merge ps_scores_strata (in=a) stratum_counts;
    by stratum;
    if a then do;
        
/* Calculate weights for treated and untreated */
    if Generation = 1 then 
         Weight = (Ntotal / &Total_N) / (Nexposed / Ntotal_exposed);
    else if Generation = 2 then 
         Weight = (Ntotal / &Total_N) / (Nunexposed / Ntotal_unexposed);
    
    /* Truncate weights above 10 */
        if Weight > 10 then Weight = 10;
    end;
run;

/*Check out the first 20 weights*/

proc print data=ps_scores_weighted(obs=20);
run;

/*Check summary statistics of the weight*/

proc means data=ps_scores_weighted n mean median min max stddev;
    var Weight;
    class Generation; /* Separate summary by exposure group */
run;

/*Check histogram of the weights*/

proc sgplot data=ps_scores_weighted;
    histogram Weight / group=Generation transparency=0.5;
    xaxis label="Weights";
    yaxis label="Frequency";
run;

proc sgplot data=ps_scores_weighted;
    density Weight / group=Generation type=kernel transparency=0.5;
    keylegend / position=topright; /* Legend for clarity */
    xaxis label="Weights";
    yaxis label="Density";
    title "Weight Distribution by Group (Area Plot)";
run;

/*Check how many weights are >10*/

proc sql;
    select count(*) as WeightsAbove10
    from ps_scores_weighted
    where Weight > 10;
quit;

/*Sort and Merge the dataset*/

proc sort data=ps_scores_weighted; by id; run;
proc sort data=final_event; by id; run;

/*Merge*/

data finally;
    merge ps_scores_weighted (in=a) final_event;
    by id;
    if a;
run;

/*Finally, run the cox model*/

proc phreg data=finally;
    class Generation(ref="2") / param=ref;
    model time*event(0) = Generation / ties=efron;
    strata combined_band; /* Stratify by combined age and calendar time bands */
    weight Weight; /* Apply fine stratification weights */
    hazardratio Generation / cl=wald; /* Include confidence intervals */
run;

/*Get incidence rates, person times and CIs*/

/*Now, calculate your incidence rate*/
proc means data = finally;					/*Calls on the proc means procedure for basic descriptive statistics*/
var event time_years;						/*Specifies the variables we are intetested in to clculate the incidence rate*/
class Generation;							/*Tells SAS to do the cuaculation for each treatment group*/
weight Weight; /* Apply weights */
output out = Inc_rate sum = event py;		/*Creates a new dataset (Irate) to store the results we are about to generate*/
											/*The sum part provides names for the sum of the variables listed above (events and follow-up)*/
quit;

data Inc_rate2;													/*Create a new dataset named 'IR'*/
set Inc_rate;													/*Tell SAS to use data from IRprep*/
Inc_rate2 = event/(py/1000);									/*Calculate the incidence rate*/
LCI=quantile('chisq',0.025, event*2)/((py/100000)*2);			/*Calculate lower confidence interval*/
UCI=quantile('chisq',0.975, (event+1)*2)/((py/100000)*2);run;	/*calculate the upper confidence interval*/
run;

proc print data=Inc_rate2;
run;

/*...................................Past....................................................*/

/*PART 3: Estimate PS and deal with confounding*/

/*Past*/

/*Extract the year variable*/

data merged_data_new;
    set merged_data_past;
    t0_year = year(t0); /* Extract the year from the t0 date */
run;

/* Print the first 10 observations to verify */

proc print data=merged_data_new(obs=10); 
    var t0 t0_year;
run;

/*We will create 5-year calendar time bands*/

data merged_data_bands;
    set merged_data_new;
    
    /* Create calendar time bands */
    if t0_year >= 1995 and t0_year < 2000 then time_band = "1995-1999";
    else if t0_year >= 2000 and t0_year < 2005 then time_band = "2000-2004";
    else if t0_year >= 2005 and t0_year < 2010 then time_band = "2005-2009";
    else time_band = "2010-2015";

    /* Combine age and calendar time bands */
    combined_band = time_band;
run;

/*Print a few rows to confirm*/

proc print data=merged_data_bands(obs=10); 
    var time_band t0_year combined_band;
run;

proc sort data=merged_data_bands; by combined_band; run;

/*Estimate the PS*/

proc logistic data=merged_data_bands;
	
	/*Binary and categorical variables*/
    class obesity(ref="0") cancer(ref="0") alcohol(ref="0") male(ref="0") 
          dementia(ref="0") nsaids(ref="0") aspirin(ref="0") pvd(ref="0") cad(ref="0") revasc(ref="0") 
          nitrates(ref="0") chf(ref="0") cerebrovascular(ref="0") thyroid(ref="0") statins(ref="0") 
          renal(ref="0") digoxin(ref="0") fibrates(ref="0") clopidogrel(ref="0") warfarin(ref="0") 
          opioids(ref="0") af(ref="0") copd(ref="0") asthma(ref="0") dyslipidemia(ref="0") anaemia(ref="0") 
          angina(ref="0") htn(ref="0") hypotension(ref="0") diabetes(ref="0") mi(ref="0") stroke(ref="0") / param=ref;

    /* Include quadratic and cubic terms for continuous variables */
    model Generation(event="1") = age age*age age*age*age visit visit*visit visit*visit*visit
                                  obesity smoking alcohol male dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
    by combined_band;
    output out=ps_scores p=pscore;
run;

/*Print 30 rows to see propensity scores*/

proc print data=ps_scores(obs=30); run;

/*View overlap*/

proc sgplot data=ps_scores;
    histogram pscore / group=Generation transparency=0.5;
run;


/*Now, we will conduct fine stratification weighting*/

/*Create 50 strata*/

proc rank data=ps_scores out=ps_scores_strata groups=50;
    var pscore; /* The propensity score variable */
    ranks stratum; /* Creates the stratum variable (0 to 49) */
run;

/* Adjust stratum to 1-based index for easier interpretation*/

data ps_scores_strata;
    set ps_scores_strata;
    stratum = stratum + 1;
run;

/*Calculate population counts for each stratum*/

proc sql;
    create table stratum_counts as
    select stratum,
           count(*) as Ntotal,
           sum(Generation = 1) as Nexposed, /* Assuming Generation=1 is the exposed group */
           sum(Generation = 2) as Nunexposed, /* Assuming Generation=2 is the unexposed group */
           (select count(*) from ps_scores_strata where Generation = 1) as Ntotal_exposed,
           (select count(*) from ps_scores_strata where Generation = 2) as Ntotal_unexposed
    from ps_scores_strata
    group by stratum;
quit;

/*Calculate weights*/

/*First, calculate the total number of obs (sum(Ntotal)) across all strata and save it as a macro variable*/

proc sql noprint;
    select sum(Ntotal) into :Total_N from stratum_counts;
quit;

/* Check the macro variable value */

%put &=Total_N;

/*Use precomputed &Total_N to calculate the weights (but sort first)*/

proc sort data=ps_scores_strata; by stratum; run;

proc sort data=stratum_counts; by stratum; run;

/*Calculate PS Weights*/

data ps_scores_weighted;
    merge ps_scores_strata (in=a) stratum_counts;
    by stratum;
    if a then do;
        
/* Calculate weights for treated and untreated */
    if Generation = 1 then 
         Weight = (Ntotal / &Total_N) / (Nexposed / Ntotal_exposed);
    else if Generation = 2 then 
         Weight = (Ntotal / &Total_N) / (Nunexposed / Ntotal_unexposed);
    
    /* Truncate weights above 10 */
        if Weight > 10 then Weight = 10;
    end;
run;

/*Check out the first 20 weights*/

proc print data=ps_scores_weighted(obs=20);
run;

/*Check summary statistics of the weight*/

proc means data=ps_scores_weighted n mean median min max stddev;
    var Weight;
    class Generation; /* Separate summary by exposure group */
run;

/*Check histogram of the weights*/

proc sgplot data=ps_scores_weighted;
    histogram Weight / group=Generation transparency=0.5;
    xaxis label="Weights";
    yaxis label="Frequency";
run;

proc sgplot data=ps_scores_weighted;
    density Weight / group=Generation type=kernel transparency=0.5;
    keylegend / position=topright; /* Legend for clarity */
    xaxis label="Weights";
    yaxis label="Density";
    title "Weight Distribution by Group (Area Plot)";
run;

/*Check how many weights are >10*/

proc sql;
    select count(*) as WeightsAbove10
    from ps_scores_weighted
    where Weight > 10;
quit;

/*Sort and Merge the dataset*/

proc sort data=ps_scores_weighted; by id; run;
proc sort data=final_event; by id; run;

/*Merge*/

data finally;
    merge ps_scores_weighted (in=a) final_event;
    by id;
    if a;
run;

/*Finally, run the cox model*/

proc phreg data=finally;
    class Generation(ref="2") / param=ref;
    model time*event(0) = Generation / ties=efron;
    strata combined_band; /* Stratify by combined age and calendar time bands */
    weight Weight; /* Apply fine stratification weights */
    hazardratio Generation / cl=wald; /* Include confidence intervals */
run;

/*Get incidence rates, person times and CIs*/

/*Now, calculate your incidence rate*/
proc means data = finally;					/*Calls on the proc means procedure for basic descriptive statistics*/
var event time_years;						/*Specifies the variables we are intetested in to clculate the incidence rate*/
class Generation;							/*Tells SAS to do the cuaculation for each treatment group*/
weight Weight; /* Apply weights */
output out = Inc_rate sum = event py;		/*Creates a new dataset (Irate) to store the results we are about to generate*/
											/*The sum part provides names for the sum of the variables listed above (events and follow-up)*/
quit;

data Inc_rate2;													/*Create a new dataset named 'IR'*/
set Inc_rate;													/*Tell SAS to use data from IRprep*/
Inc_rate2 = event/(py/1000);									/*Calculate the incidence rate*/
LCI=quantile('chisq',0.025, event*2)/((py/100000)*2);			/*Calculate lower confidence interval*/
UCI=quantile('chisq',0.975, (event+1)*2)/((py/100000)*2);run;	/*calculate the upper confidence interval*/
run;

proc print data=Inc_rate2;
run;

/*...................................Current....................................................*/

/*PART 3: Estimate PS and deal with confounding*/

/*Past*/

/*Extract the year variable*/

data merged_data_new;
    set merged_data_current;
    t0_year = year(t0); /* Extract the year from the t0 date */
run;

/* Print the first 10 observations to verify */

proc print data=merged_data_new(obs=10); 
    var t0 t0_year;
run;

/*We will create 5-year calendar time bands*/

data merged_data_bands;
    set merged_data_new;
    
    /* Create calendar time bands */
    if t0_year >= 1995 and t0_year < 2000 then time_band = "1995-1999";
    else if t0_year >= 2000 and t0_year < 2005 then time_band = "2000-2004";
    else if t0_year >= 2005 and t0_year < 2010 then time_band = "2005-2009";
    else time_band = "2010-2015";

    /* Combine age and calendar time bands */
    combined_band = time_band;
run;

/*Print a few rows to confirm*/

proc print data=merged_data_bands(obs=10); 
    var time_band t0_year combined_band;
run;

proc sort data=merged_data_bands; by combined_band; run;

/*Estimate the PS*/

proc logistic data=merged_data_bands;
	
	/*Binary and categorical variables*/
    class obesity(ref="0") cancer(ref="0") alcohol(ref="0") male(ref="0") 
          dementia(ref="0") nsaids(ref="0") aspirin(ref="0") pvd(ref="0") cad(ref="0") revasc(ref="0") 
          nitrates(ref="0") chf(ref="0") cerebrovascular(ref="0") thyroid(ref="0") statins(ref="0") 
          renal(ref="0") digoxin(ref="0") fibrates(ref="0") clopidogrel(ref="0") warfarin(ref="0") 
          opioids(ref="0") af(ref="0") copd(ref="0") asthma(ref="0") dyslipidemia(ref="0") anaemia(ref="0") 
          angina(ref="0") htn(ref="0") hypotension(ref="0") diabetes(ref="0") mi(ref="0") stroke(ref="0") / param=ref;

    /* Include quadratic and cubic terms for continuous variables */
    model Generation(event="1") = age age*age age*age*age visit visit*visit visit*visit*visit
                                  obesity smoking alcohol male dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
    by combined_band;
    output out=ps_scores p=pscore;
run;

/*Print 30 rows to see propensity scores*/

proc print data=ps_scores(obs=30); run;

/*View overlap*/

proc sgplot data=ps_scores;
    histogram pscore / group=Generation transparency=0.5;
run;


/*Now, we will conduct fine stratification weighting*/

/*Create 50 strata*/

proc rank data=ps_scores out=ps_scores_strata groups=50;
    var pscore; /* The propensity score variable */
    ranks stratum; /* Creates the stratum variable (0 to 49) */
run;

/* Adjust stratum to 1-based index for easier interpretation*/

data ps_scores_strata;
    set ps_scores_strata;
    stratum = stratum + 1;
run;

/*Calculate population counts for each stratum*/

proc sql;
    create table stratum_counts as
    select stratum,
           count(*) as Ntotal,
           sum(Generation = 1) as Nexposed, /* Assuming Generation=1 is the exposed group */
           sum(Generation = 2) as Nunexposed, /* Assuming Generation=2 is the unexposed group */
           (select count(*) from ps_scores_strata where Generation = 1) as Ntotal_exposed,
           (select count(*) from ps_scores_strata where Generation = 2) as Ntotal_unexposed
    from ps_scores_strata
    group by stratum;
quit;

/*Calculate weights*/

/*First, calculate the total number of obs (sum(Ntotal)) across all strata and save it as a macro variable*/

proc sql noprint;
    select sum(Ntotal) into :Total_N from stratum_counts;
quit;

/* Check the macro variable value */

%put &=Total_N;

/*Use precomputed &Total_N to calculate the weights (but sort first)*/

proc sort data=ps_scores_strata; by stratum; run;

proc sort data=stratum_counts; by stratum; run;

/*Calculate PS Weights*/

data ps_scores_weighted;
    merge ps_scores_strata (in=a) stratum_counts;
    by stratum;
    if a then do;
        
/* Calculate weights for treated and untreated */
    if Generation = 1 then 
         Weight = (Ntotal / &Total_N) / (Nexposed / Ntotal_exposed);
    else if Generation = 2 then 
         Weight = (Ntotal / &Total_N) / (Nunexposed / Ntotal_unexposed);
    
    /* Truncate weights above 10 */
        if Weight > 10 then Weight = 10;
    end;
run;

/*Check out the first 20 weights*/

proc print data=ps_scores_weighted(obs=20);
run;

/*Check summary statistics of the weight*/

proc means data=ps_scores_weighted n mean median min max stddev;
    var Weight;
    class Generation; /* Separate summary by exposure group */
run;

/*Check histogram of the weights*/

proc sgplot data=ps_scores_weighted;
    histogram Weight / group=Generation transparency=0.5;
    xaxis label="Weights";
    yaxis label="Frequency";
run;

proc sgplot data=ps_scores_weighted;
    density Weight / group=Generation type=kernel transparency=0.5;
    keylegend / position=topright; /* Legend for clarity */
    xaxis label="Weights";
    yaxis label="Density";
    title "Weight Distribution by Group (Area Plot)";
run;

/*Check how many weights are >10*/

proc sql;
    select count(*) as WeightsAbove10
    from ps_scores_weighted
    where Weight > 10;
quit;

/*Sort and Merge the dataset*/

proc sort data=ps_scores_weighted; by id; run;
proc sort data=final_event; by id; run;

/*Merge*/

data finally;
    merge ps_scores_weighted (in=a) final_event;
    by id;
    if a;
run;

/*Finally, run the cox model*/

proc phreg data=finally;
    class Generation(ref="2") / param=ref;
    model time*event(0) = Generation / ties=efron;
    strata combined_band; /* Stratify by combined age and calendar time bands */
    weight Weight; /* Apply fine stratification weights */
    hazardratio Generation / cl=wald; /* Include confidence intervals */
run;

/*Get incidence rates, person times and CIs*/

/*Now, calculate your incidence rate*/
proc means data = finally;					/*Calls on the proc means procedure for basic descriptive statistics*/
var event time_years;						/*Specifies the variables we are intetested in to clculate the incidence rate*/
class Generation;							/*Tells SAS to do the cuaculation for each treatment group*/
weight Weight; /* Apply weights */
output out = Inc_rate sum = event py;		/*Creates a new dataset (Irate) to store the results we are about to generate*/
											/*The sum part provides names for the sum of the variables listed above (events and follow-up)*/
quit;

data Inc_rate2;													/*Create a new dataset named 'IR'*/
set Inc_rate;													/*Tell SAS to use data from IRprep*/
Inc_rate2 = event/(py/1000);									/*Calculate the incidence rate*/
LCI=quantile('chisq',0.025, event*2)/((py/100000)*2);			/*Calculate lower confidence interval*/
UCI=quantile('chisq',0.975, (event+1)*2)/((py/100000)*2);run;	/*calculate the upper confidence interval*/
run;

proc print data=Inc_rate2;
run;


/*.............................Duration of Use........................................*/

/*.............................Duration of Use........................................*/


/*Split data into males and female and basically repeat eveything for each group*/

data merged_data_duration1 merged_data_duration2;
    set merged_data_cov;
    if time_years > 1 then output merged_data_duration1; 
    if time_years >= 3 then output merged_data_duration2; 
run;

/*PART 3: Estimate PS and deal with confounding*/

/*Never*/

/*Extract the year variable*/

data merged_data_new;
    set merged_data_duration1;
    t0_year = year(t0); /* Extract the year from the t0 date */
run;

/* Print the first 10 observations to verify */

proc print data=merged_data_new(obs=10); 
    var t0 t0_year;
run;

/*We will create 5-year calendar time bands*/

data merged_data_bands;
    set merged_data_new;
    
    /* Create calendar time bands */
    if t0_year >= 1995 and t0_year < 2000 then time_band = "1995-1999";
    else if t0_year >= 2000 and t0_year < 2005 then time_band = "2000-2004";
    else if t0_year >= 2005 and t0_year < 2010 then time_band = "2005-2009";
    else time_band = "2010-2015";

    /* Combine age and calendar time bands */
    combined_band = time_band;
run;

/*Print a few rows to confirm*/

proc print data=merged_data_bands(obs=10); 
    var time_band t0_year combined_band;
run;

proc sort data=merged_data_bands; by combined_band; run;

/*Estimate the PS*/

proc logistic data=merged_data_bands;
	
	/*Binary and categorical variables*/
    class obesity(ref="0") cancer(ref="0") alcohol(ref="0") male(ref="0") smoking(ref="0")
          dementia(ref="0") nsaids(ref="0") aspirin(ref="0") pvd(ref="0") cad(ref="0") revasc(ref="0") 
          nitrates(ref="0") chf(ref="0") cerebrovascular(ref="0") thyroid(ref="0") statins(ref="0") 
          renal(ref="0") digoxin(ref="0") fibrates(ref="0") clopidogrel(ref="0") warfarin(ref="0") 
          opioids(ref="0") af(ref="0") copd(ref="0") asthma(ref="0") dyslipidemia(ref="0") anaemia(ref="0") 
          angina(ref="0") htn(ref="0") hypotension(ref="0") diabetes(ref="0") mi(ref="0") stroke(ref="0") / param=ref;

    /* Include quadratic and cubic terms for continuous variables */
    model Generation(event="1") = age age*age age*age*age visit visit*visit visit*visit*visit
                                  obesity smoking alcohol male dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
    by combined_band;
    output out=ps_scores p=pscore;
run;

/*Print 30 rows to see propensity scores*/

proc print data=ps_scores(obs=30); run;

/*View overlap*/

proc sgplot data=ps_scores;
    histogram pscore / group=Generation transparency=0.5;
run;


/*Now, we will conduct fine stratification weighting*/

/*Create 50 strata*/

proc rank data=ps_scores out=ps_scores_strata groups=50;
    var pscore; /* The propensity score variable */
    ranks stratum; /* Creates the stratum variable (0 to 49) */
run;

/* Adjust stratum to 1-based index for easier interpretation*/

data ps_scores_strata;
    set ps_scores_strata;
    stratum = stratum + 1;
run;

/*Calculate population counts for each stratum*/

proc sql;
    create table stratum_counts as
    select stratum,
           count(*) as Ntotal,
           sum(Generation = 1) as Nexposed, /* Assuming Generation=1 is the exposed group */
           sum(Generation = 2) as Nunexposed, /* Assuming Generation=2 is the unexposed group */
           (select count(*) from ps_scores_strata where Generation = 1) as Ntotal_exposed,
           (select count(*) from ps_scores_strata where Generation = 2) as Ntotal_unexposed
    from ps_scores_strata
    group by stratum;
quit;

/*Calculate weights*/

/*First, calculate the total number of obs (sum(Ntotal)) across all strata and save it as a macro variable*/

proc sql noprint;
    select sum(Ntotal) into :Total_N from stratum_counts;
quit;

/* Check the macro variable value */

%put &=Total_N;

/*Use precomputed &Total_N to calculate the weights (but sort first)*/

proc sort data=ps_scores_strata; by stratum; run;

proc sort data=stratum_counts; by stratum; run;

/*Calculate PS Weights*/

data ps_scores_weighted;
    merge ps_scores_strata (in=a) stratum_counts;
    by stratum;
    if a then do;
        
/* Calculate weights for treated and untreated */
    if Generation = 1 then 
         Weight = (Ntotal / &Total_N) / (Nexposed / Ntotal_exposed);
    else if Generation = 2 then 
         Weight = (Ntotal / &Total_N) / (Nunexposed / Ntotal_unexposed);
    
    /* Truncate weights above 10 */
        if Weight > 10 then Weight = 10;
    end;
run;

/*Check out the first 20 weights*/

proc print data=ps_scores_weighted(obs=20);
run;

/*Check summary statistics of the weight*/

proc means data=ps_scores_weighted n mean median min max stddev;
    var Weight;
    class Generation; /* Separate summary by exposure group */
run;

/*Check histogram of the weights*/

proc sgplot data=ps_scores_weighted;
    histogram Weight / group=Generation transparency=0.5;
    xaxis label="Weights";
    yaxis label="Frequency";
run;

proc sgplot data=ps_scores_weighted;
    density Weight / group=Generation type=kernel transparency=0.5;
    keylegend / position=topright; /* Legend for clarity */
    xaxis label="Weights";
    yaxis label="Density";
    title "Weight Distribution by Group (Area Plot)";
run;

/*Check how many weights are >10*/

proc sql;
    select count(*) as WeightsAbove10
    from ps_scores_weighted
    where Weight > 10;
quit;

/*Sort and Merge the dataset*/

proc sort data=ps_scores_weighted; by id; run;
proc sort data=final_event; by id; run;

/*Merge*/

data finally;
    merge ps_scores_weighted (in=a) final_event;
    by id;
    if a;
run;

/*Finally, run the cox model*/

proc phreg data=finally;
    class Generation(ref="2") / param=ref;
    model time*event(0) = Generation / ties=efron;
    strata combined_band; /* Stratify by combined age and calendar time bands */
    weight Weight; /* Apply fine stratification weights */
    hazardratio Generation / cl=wald; /* Include confidence intervals */
run;

/*Get incidence rates, person times and CIs*/

/*Now, calculate your incidence rate*/
proc means data = finally;					/*Calls on the proc means procedure for basic descriptive statistics*/
var event time_years;						/*Specifies the variables we are intetested in to clculate the incidence rate*/
class Generation;							/*Tells SAS to do the cuaculation for each treatment group*/
weight Weight; /* Apply weights */
output out = Inc_rate sum = event py;		/*Creates a new dataset (Irate) to store the results we are about to generate*/
											/*The sum part provides names for the sum of the variables listed above (events and follow-up)*/
quit;

data Inc_rate2;													/*Create a new dataset named 'IR'*/
set Inc_rate;													/*Tell SAS to use data from IRprep*/
Inc_rate2 = event/(py/1000);									/*Calculate the incidence rate*/
LCI=quantile('chisq',0.025, event*2)/((py/100000)*2);			/*Calculate lower confidence interval*/
UCI=quantile('chisq',0.975, (event+1)*2)/((py/100000)*2);run;	/*calculate the upper confidence interval*/
run;

proc print data=Inc_rate2;
run;

/*.................................Duration 2.....................................*/

/*PART 3: Estimate PS and deal with confounding*/

/*Never*/

/*Extract the year variable*/

data merged_data_new;
    set merged_data_duration2;
    t0_year = year(t0); /* Extract the year from the t0 date */
run;

/* Print the first 10 observations to verify */

proc print data=merged_data_new(obs=10); 
    var t0 t0_year;
run;

/*We will create 5-year calendar time bands*/

data merged_data_bands;
    set merged_data_new;
    
    /* Create calendar time bands */
    if t0_year >= 1995 and t0_year < 2000 then time_band = "1995-1999";
    else if t0_year >= 2000 and t0_year < 2005 then time_band = "2000-2004";
    else if t0_year >= 2005 and t0_year < 2010 then time_band = "2005-2009";
    else time_band = "2010-2015";

    /* Combine age and calendar time bands */
    combined_band = time_band;
run;

/*Print a few rows to confirm*/

proc print data=merged_data_bands(obs=10); 
    var time_band t0_year combined_band;
run;

proc sort data=merged_data_bands; by combined_band; run;

/*Estimate the PS*/

proc logistic data=merged_data_bands;
	
	/*Binary and categorical variables*/
    class obesity(ref="0") cancer(ref="0") alcohol(ref="0") male(ref="0") smoking(ref="0")
          dementia(ref="0") nsaids(ref="0") aspirin(ref="0") pvd(ref="0") cad(ref="0") revasc(ref="0") 
          nitrates(ref="0") chf(ref="0") cerebrovascular(ref="0") thyroid(ref="0") statins(ref="0") 
          renal(ref="0") digoxin(ref="0") fibrates(ref="0") clopidogrel(ref="0") warfarin(ref="0") 
          opioids(ref="0") af(ref="0") copd(ref="0") asthma(ref="0") dyslipidemia(ref="0") anaemia(ref="0") 
          angina(ref="0") htn(ref="0") hypotension(ref="0") diabetes(ref="0") mi(ref="0") stroke(ref="0") / param=ref;

    /* Include quadratic and cubic terms for continuous variables */
    model Generation(event="1") = age age*age age*age*age visit visit*visit visit*visit*visit
                                  obesity smoking alcohol male dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
    by combined_band;
    output out=ps_scores p=pscore;
run;

/*Print 30 rows to see propensity scores*/

proc print data=ps_scores(obs=30); run;

/*View overlap*/

proc sgplot data=ps_scores;
    histogram pscore / group=Generation transparency=0.5;
run;


/*Now, we will conduct fine stratification weighting*/

/*Create 50 strata*/

proc rank data=ps_scores out=ps_scores_strata groups=50;
    var pscore; /* The propensity score variable */
    ranks stratum; /* Creates the stratum variable (0 to 49) */
run;

/* Adjust stratum to 1-based index for easier interpretation*/

data ps_scores_strata;
    set ps_scores_strata;
    stratum = stratum + 1;
run;

/*Calculate population counts for each stratum*/

proc sql;
    create table stratum_counts as
    select stratum,
           count(*) as Ntotal,
           sum(Generation = 1) as Nexposed, /* Assuming Generation=1 is the exposed group */
           sum(Generation = 2) as Nunexposed, /* Assuming Generation=2 is the unexposed group */
           (select count(*) from ps_scores_strata where Generation = 1) as Ntotal_exposed,
           (select count(*) from ps_scores_strata where Generation = 2) as Ntotal_unexposed
    from ps_scores_strata
    group by stratum;
quit;

/*Calculate weights*/

/*First, calculate the total number of obs (sum(Ntotal)) across all strata and save it as a macro variable*/

proc sql noprint;
    select sum(Ntotal) into :Total_N from stratum_counts;
quit;

/* Check the macro variable value */

%put &=Total_N;

/*Use precomputed &Total_N to calculate the weights (but sort first)*/

proc sort data=ps_scores_strata; by stratum; run;

proc sort data=stratum_counts; by stratum; run;

/*Calculate PS Weights*/

data ps_scores_weighted;
    merge ps_scores_strata (in=a) stratum_counts;
    by stratum;
    if a then do;
        
/* Calculate weights for treated and untreated */
    if Generation = 1 then 
         Weight = (Ntotal / &Total_N) / (Nexposed / Ntotal_exposed);
    else if Generation = 2 then 
         Weight = (Ntotal / &Total_N) / (Nunexposed / Ntotal_unexposed);
    
    /* Truncate weights above 10 */
        if Weight > 10 then Weight = 10;
    end;
run;

/*Check out the first 20 weights*/

proc print data=ps_scores_weighted(obs=20);
run;

/*Check summary statistics of the weight*/

proc means data=ps_scores_weighted n mean median min max stddev;
    var Weight;
    class Generation; /* Separate summary by exposure group */
run;

/*Check histogram of the weights*/

proc sgplot data=ps_scores_weighted;
    histogram Weight / group=Generation transparency=0.5;
    xaxis label="Weights";
    yaxis label="Frequency";
run;

proc sgplot data=ps_scores_weighted;
    density Weight / group=Generation type=kernel transparency=0.5;
    keylegend / position=topright; /* Legend for clarity */
    xaxis label="Weights";
    yaxis label="Density";
    title "Weight Distribution by Group (Area Plot)";
run;

/*Check how many weights are >10*/

proc sql;
    select count(*) as WeightsAbove10
    from ps_scores_weighted
    where Weight > 10;
quit;

/*Sort and Merge the dataset*/

proc sort data=ps_scores_weighted; by id; run;
proc sort data=final_event; by id; run;

/*Merge*/

data finally;
    merge ps_scores_weighted (in=a) final_event;
    by id;
    if a;
run;

/*Finally, run the cox model*/

proc phreg data=finally;
    class Generation(ref="2") / param=ref;
    model time*event(0) = Generation / ties=efron;
    strata combined_band; /* Stratify by combined age and calendar time bands */
    weight Weight; /* Apply fine stratification weights */
    hazardratio Generation / cl=wald; /* Include confidence intervals */
run;

/*Get incidence rates, person times and CIs*/

/*Now, calculate your incidence rate*/
proc means data = finally;					/*Calls on the proc means procedure for basic descriptive statistics*/
var event time_years;						/*Specifies the variables we are intetested in to clculate the incidence rate*/
class Generation;							/*Tells SAS to do the cuaculation for each treatment group*/
weight Weight; /* Apply weights */
output out = Inc_rate sum = event py;		/*Creates a new dataset (Irate) to store the results we are about to generate*/
											/*The sum part provides names for the sum of the variables listed above (events and follow-up)*/
quit;

data Inc_rate2;													/*Create a new dataset named 'IR'*/
set Inc_rate;													/*Tell SAS to use data from IRprep*/
Inc_rate2 = event/(py/1000);									/*Calculate the incidence rate*/
LCI=quantile('chisq',0.025, event*2)/((py/100000)*2);			/*Calculate lower confidence interval*/
UCI=quantile('chisq',0.975, (event+1)*2)/((py/100000)*2);run;	/*calculate the upper confidence interval*/
run;

proc print data=Inc_rate2;
run;


/*..............................Sensitivity Analysis........................................*/

/*.............................Sensitivity Analyses...........................................*/



/*........................................Sensitivity Analyses...................................*/

/*........................................Sensitivity Analyses...................................*/

/*Varying the lag to 2 years*/

/*Start from PART 2 of the main analysis: Define the Exposure, Lag, Grace periods, Censoring dates etc*/

data cohort1;
 set final_rx;
 by id date;
 if first.id;		/*We select only the first id and assign the dates time 0*/
 rename date=t0;
run;

/*Identify switchers and determine the switch date*/

data switch;
	merge cohort1(in=a) final_rx(where=(date <'31dec2015'd) rename=(Generation=gen));
	by id;
	if a;
	if Generation ne gen;		/*Only the people switch make it into this dataset*/
run;

/*Create a dataset with only the first patient ID and the date of censoring*/

data switch_data;
	set switch;
	by id date;
	if first.id;			/*Keep only first id. Remember you merged with longform final_rx above*/
	keep id date;			/*Keep these 2 variables and drop everything else*/
	rename date=censor1;	
run;

/*Merge 3 datasets. So we now have t0 from 'cohort1' and censoring date (censor1) from 'switch_data' 
added to the long form (final_rx). T0 is fixed for each person while date varies (date of every refill rx)*/

data continuous_data;
	merge cohort1(in=a) final_rx switch_data;
	by id;
	if a;
	if date < censor1 or censor1=.;
run;

/*Apply rxdate + 30 + grace period (365 days)*/

data c1;
	set continuous_data;
	by id date;
	rxdur=date+30+365;				/*rxdur id prescription duration plus grace period*/
	format rxdur lagrxdur date9.;
	lagrxdur=lag(rxdur);			/*lagrxdur is the end date for the previous rx*/
	if first.id then lagrxdur=date;
	retain discon;
	if first.id then discon=0;
	if date > lagrxdur then discon=1;
	if discon=1 and id=lag(id) and discon=lag(discon) then delete; 
run;

/*Calculate the final censoring date (censor2) for each patient. 'Censor1' only identified 
switchers by capturing their first switching date. Censor two accounts for grace periods in both switchers and non-switchers.*/

data c2;
 	set c1;
 	by id date;
 	format censor2 date9.;
 	if discon = 1 then censor2 = lagrxdur + 365; /*For switchers: Added the 365 for the grace period post rx cessation*/
 	else if last.id then censor2 = rxdur + 365;	 /*For those who stop completely: Added the 365 for the grace period post rx cessation*/
 	if last.id;
run;

/*Now, define and merge the events*/

data merged_data;
    merge final_event (in=a) c2(in=b) final_patients(keep=id end); /* Add patients to get 'end' date*/
    by id;
    if b; 
    end_followup = min(eventdate,end,'31dec2015'd,censor2);
    format end date9.;
    if a and t0 < eventdate <= end_followup then event=1; 
	else event = 0;
	time = end_followup - t0 + 1;
	time_years = time/365.25;
	if time < 730 then event=0; /*Apply 1 year lag where we don't count events*/
run;


/*PART 3: Estimate PS and deal with confounding*/

/*First, merge the covariates dataset with the current dataset (but sort first)*/

proc sort data=merged_data; by id; run;
proc sort data=final_covariates; by id; run;

/* Step 2: Merge datasets */

data merged_data2;
    merge merged_data (in=a) final_covariates (in=b);
    by id;
    if a; 	/* Keeps only rows that exist in the main dataset*/
run;

/*Extract the year variable*/

data merged_data3;
    set merged_data2;
    t0_year = year(t0); /* Extract the year from the t0 date */
run;

/* Print the first 10 observations to verify */

proc print data=merged_data3(obs=10); 
    var t0 t0_year;
run;

/*We will create 5-year calendar time bands*/

data merged_data_bands;
    set merged_data3;
    
    /* Create calendar time bands */
    if t0_year >= 1995 and t0_year < 2000 then time_band = "1995-1999";
    else if t0_year >= 2000 and t0_year < 2005 then time_band = "2000-2004";
    else if t0_year >= 2005 and t0_year < 2010 then time_band = "2005-2009";
    else time_band = "2010-2015";

    /* Combine age and calendar time bands */
    combined_band = time_band;
run;

/*Print a few rows to confirm*/

proc print data=merged_data_bands(obs=10); 
    var time_band t0_year combined_band;
run;

proc sort data=merged_data_bands; by combined_band; run;

/*Estimate the PS*/

proc logistic data=merged_data_bands;
	
	/*Binary and categorical variables*/
    class obesity(ref="0") smoking(ref="0") male(ref="0") alcohol(ref="0") cancer(ref="0") 
          dementia(ref="0") nsaids(ref="0") aspirin(ref="0") pvd(ref="0") cad(ref="0") revasc(ref="0") 
          nitrates(ref="0") chf(ref="0") cerebrovascular(ref="0") thyroid(ref="0") statins(ref="0") 
          renal(ref="0") digoxin(ref="0") fibrates(ref="0") clopidogrel(ref="0") warfarin(ref="0") 
          opioids(ref="0") af(ref="0") copd(ref="0") asthma(ref="0") dyslipidemia(ref="0") anaemia(ref="0") 
          angina(ref="0") htn(ref="0") hypotension(ref="0") diabetes(ref="0") mi(ref="0") stroke(ref="0") / param=ref;

    /* Include quadratic and cubic terms for continuous variables */
    model Generation(event="1") = age age*age age*age*age visit visit*visit visit*visit*visit
                                  obesity smoking male alcohol cancer dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
    by combined_band;
    output out=ps_scores p=pscore;
run;

/*Print 30 rows to see propensity scores*/

proc print data=ps_scores(obs=30); run;

/*View overlap*/

proc sgplot data=ps_scores;
    histogram pscore / group=Generation transparency=0.5;
run;


/*Now, we will conduct fine stratification weighting*/

/*Create 50 strata*/

proc rank data=ps_scores out=ps_scores_strata groups=50;
    var pscore; /* The propensity score variable */
    ranks stratum; /* Creates the stratum variable (0 to 49) */
run;

/* Adjust stratum to 1-based index for easier interpretation*/

data ps_scores_strata;
    set ps_scores_strata;
    stratum = stratum + 1;
run;

/*Calculate population counts for each stratum*/

proc sql;
    create table stratum_counts as
    select stratum,
           count(*) as Ntotal,
           sum(Generation = 1) as Nexposed, /* Assuming Generation=1 is the exposed group */
           sum(Generation = 2) as Nunexposed, /* Assuming Generation=2 is the unexposed group */
           (select count(*) from ps_scores_strata where Generation = 1) as Ntotal_exposed,
           (select count(*) from ps_scores_strata where Generation = 2) as Ntotal_unexposed
    from ps_scores_strata
    group by stratum;
quit;

/*Calculate weights*/

/*First, calculate the total number of obs (sum(Ntotal)) across all strata and save it as a macro variable*/

proc sql noprint;
    select sum(Ntotal) into :Total_N from stratum_counts;
quit;

/* Check the macro variable value */

%put &=Total_N;

/*Use precomputed &Total_N to calculate the weights (but sort first)*/

proc sort data=ps_scores_strata; by stratum; run;

proc sort data=stratum_counts; by stratum; run;

/*Calculate PS Weights*/

data ps_scores_weighted;
    merge ps_scores_strata (in=a) stratum_counts;
    by stratum;
    if a then do;
        
/* Calculate weights for treated and untreated */
    if Generation = 1 then 
         Weight = (Ntotal / &Total_N) / (Nexposed / Ntotal_exposed);
    else if Generation = 2 then 
         Weight = (Ntotal / &Total_N) / (Nunexposed / Ntotal_unexposed);
    
    /* Truncate weights above 10 */
        if Weight > 10 then Weight = 10;
    end;
run;

/*Check out the first 20 weights*/

proc print data=ps_scores_weighted(obs=20);
run;

/*Check summary statistics of the weight*/

proc means data=ps_scores_weighted n mean median min max stddev;
    var Weight;
    class Generation; /* Separate summary by exposure group */
run;

/*Check histogram of the weights*/

proc sgplot data=ps_scores_weighted;
    histogram Weight / group=Generation transparency=0.5;
    xaxis label="Weights";
    yaxis label="Frequency";
run;

proc sgplot data=ps_scores_weighted;
    density Weight / group=Generation type=kernel transparency=0.5;
    keylegend / position=topright; /* Legend for clarity */
    xaxis label="Weights";
    yaxis label="Density";
    title "Weight Distribution by Group (Area Plot)";
run;

/*Check how many weights are >10*/

proc sql;
    select count(*) as WeightsAbove10
    from ps_scores_weighted
    where Weight > 10;
quit;

/*Sort and Merge the dataset*/

proc sort data=ps_scores_weighted; by id; run;
proc sort data=final_event; by id; run;

/*Merge*/

data finally;
    merge ps_scores_weighted (in=a) final_event;
    by id;
    if a;
run;

/*Finally, run the cox model*/

proc phreg data=finally;
    class Generation(ref="2") / param=ref;
    model time*event(0) = Generation / ties=efron;
    strata combined_band; /* Stratify by combined age and calendar time bands */
    weight Weight; /* Apply fine stratification weights */
    hazardratio Generation / cl=wald; /* Include confidence intervals */
run;

/*Get incidence rates, person times and CIs*/

/*Now, calculate your incidence rate*/
proc means data = finally;					/*Calls on the proc means procedure for basic descriptive statistics*/
var event time_years;						/*Specifies the variables we are intetested in to clculate the incidence rate*/
class Generation;							/*Tells SAS to do the cuaculation for each treatment group*/
weight Weight; 								/* Apply weights */
output out = Inc_rate sum = event py;		/*Creates a new dataset (Irate) to store the results we are about to generate*/
											/*The sum part provides names for the sum of the variables listed above (events and follow-up)*/
quit;

data Inc_rate2;													/*Create a new dataset named 'IR'*/
set Inc_rate;													/*Tell SAS to use data from IRprep*/
Inc_rate2 = event/(py/1000);									/*Calculate the incidence rate*/
LCI=quantile('chisq',0.025, event*2)/((py/100000)*2);			/*Calculate lower confidence interval*/
UCI=quantile('chisq',0.975, (event+1)*2)/((py/100000)*2);run;	/*calculate the upper confidence interval*/
run;

proc print data=Inc_rate2;
run;


/*........................................Sensitivity Analyses...................................*/

/*........................................Sensitivity Analyses...................................*/

/*Varying the lag to 3 years*/

/*PART 2: Define the Exposure, Lag, Grace periods, Censoring dates etc*/

data cohort1;
 set final_rx;
 by id date;
 if first.id;		/*We select only the first id and assign the dates time 0*/
 rename date=t0;
run;

/*Identify switchers and determine the switch date*/

data switch;
	merge cohort1(in=a) final_rx(where=(date <'31dec2015'd) rename=(Generation=gen));
	by id;
	if a;
	if Generation ne gen;		/*Only the people switch make it into this dataset*/
run;

/*Create a dataset with only the first patient ID and the date of censoring*/

data switch_data;
	set switch;
	by id date;
	if first.id;			/*Keep only first id. Remember you merged with longform final_rx above*/
	keep id date;			/*Keep these 2 variables and drop everything else*/
	rename date=censor1;	
run;

/*Merge 3 datasets. So we now have t0 from 'cohort1' and censoring date (censor1) from 'switch_data' 
added to the long form (final_rx). T0 is fixed for each person while date varies (date of every refill rx)*/

data continuous_data;
	merge cohort1(in=a) final_rx switch_data;
	by id;
	if a;
	if date < censor1 or censor1=.;
run;

/*Apply rxdate + 30 + grace period (365 days)*/

data c1;
	set continuous_data;
	by id date;
	rxdur=date+30+365;				/*rxdur id prescription duration plus grace period*/
	format rxdur lagrxdur date9.;
	lagrxdur=lag(rxdur);			/*lagrxdur is the end date for the previous rx*/
	if first.id then lagrxdur=date;
	retain discon;
	if first.id then discon=0;
	if date > lagrxdur then discon=1;
	if discon=1 and id=lag(id) and discon=lag(discon) then delete; 
run;

/*Calculate the final censoring date (censor2) for each patient. 'Censor1' only identified 
switchers by capturing their first switching date. Censor two accounts for grace periods in both switchers and non-switchers.*/

data c2;
 	set c1;
 	by id date;
 	format censor2 date9.;
 	if discon = 1 then censor2 = lagrxdur + 365; /*For switchers: Added the 365 for the grace period post rx cessation*/
 	else if last.id then censor2 = rxdur + 365;	 /*For those who stop completely: Added the 365 for the grace period post rx cessation*/
 	if last.id;
run;

/*Now, define and merge the events*/

data merged_data;
    merge final_event (in=a) c2(in=b) final_patients(keep=id end); /* Add patients to get 'end' date*/
    by id;
    if b; 
    end_followup = min(eventdate,end,'31dec2015'd,censor2);
    format end date9.;
    if a and t0 < eventdate <= end_followup then event=1; 
	else event = 0;
	time = end_followup - t0 + 1;
	time_years = time/365.25;
	if time < 1095 then event=0; /*Apply 3 year lag where we don't count events*/
run;


/*PART 3: Estimate PS and deal with confounding*/

/*First, merge the covariates dataset with the current dataset (but sort first)*/

proc sort data=merged_data; by id; run;
proc sort data=final_covariates; by id; run;

/* Step 2: Merge datasets */

data merged_data2;
    merge merged_data (in=a) final_covariates (in=b);
    by id;
    if a; 	/* Keeps only rows that exist in the main dataset*/
run;

/*Extract the year variable*/

data merged_data3;
    set merged_data2;
    t0_year = year(t0); /* Extract the year from the t0 date */
run;

/* Print the first 10 observations to verify */

proc print data=merged_data3(obs=10); 
    var t0 t0_year;
run;

/*We will create 5-year calendar time bands*/

data merged_data_bands;
    set merged_data3;
    
    /* Create calendar time bands */
    if t0_year >= 1995 and t0_year < 2000 then time_band = "1995-1999";
    else if t0_year >= 2000 and t0_year < 2005 then time_band = "2000-2004";
    else if t0_year >= 2005 and t0_year < 2010 then time_band = "2005-2009";
    else time_band = "2010-2015";

    /* Combine age and calendar time bands */
    combined_band = time_band;
run;

/*Print a few rows to confirm*/

proc print data=merged_data_bands(obs=10); 
    var time_band t0_year combined_band;
run;

proc sort data=merged_data_bands; by combined_band; run;

/*Estimate the PS*/

proc logistic data=merged_data_bands;
	
	/*Binary and categorical variables*/
    class obesity(ref="0") smoking(ref="0") male(ref="0") alcohol(ref="0") cancer(ref="0") 
          dementia(ref="0") nsaids(ref="0") aspirin(ref="0") pvd(ref="0") cad(ref="0") revasc(ref="0") 
          nitrates(ref="0") chf(ref="0") cerebrovascular(ref="0") thyroid(ref="0") statins(ref="0") 
          renal(ref="0") digoxin(ref="0") fibrates(ref="0") clopidogrel(ref="0") warfarin(ref="0") 
          opioids(ref="0") af(ref="0") copd(ref="0") asthma(ref="0") dyslipidemia(ref="0") anaemia(ref="0") 
          angina(ref="0") htn(ref="0") hypotension(ref="0") diabetes(ref="0") mi(ref="0") stroke(ref="0") / param=ref;

    /* Include quadratic and cubic terms for continuous variables */
    model Generation(event="1") = age age*age age*age*age visit visit*visit visit*visit*visit
                                  obesity smoking male alcohol cancer dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
    by combined_band;
    output out=ps_scores p=pscore;
run;

/*Print 30 rows to see propensity scores*/

proc print data=ps_scores(obs=30); run;

/*View overlap*/

proc sgplot data=ps_scores;
    histogram pscore / group=Generation transparency=0.5;
run;


/*Now, we will conduct fine stratification weighting*/

/*Create 50 strata*/

proc rank data=ps_scores out=ps_scores_strata groups=50;
    var pscore; /* The propensity score variable */
    ranks stratum; /* Creates the stratum variable (0 to 49) */
run;

/* Adjust stratum to 1-based index for easier interpretation*/

data ps_scores_strata;
    set ps_scores_strata;
    stratum = stratum + 1;
run;

/*Calculate population counts for each stratum*/

proc sql;
    create table stratum_counts as
    select stratum,
           count(*) as Ntotal,
           sum(Generation = 1) as Nexposed, /* Assuming Generation=1 is the exposed group */
           sum(Generation = 2) as Nunexposed, /* Assuming Generation=2 is the unexposed group */
           (select count(*) from ps_scores_strata where Generation = 1) as Ntotal_exposed,
           (select count(*) from ps_scores_strata where Generation = 2) as Ntotal_unexposed
    from ps_scores_strata
    group by stratum;
quit;

/*Calculate weights*/

/*First, calculate the total number of obs (sum(Ntotal)) across all strata and save it as a macro variable*/

proc sql noprint;
    select sum(Ntotal) into :Total_N from stratum_counts;
quit;

/* Check the macro variable value */

%put &=Total_N;

/*Use precomputed &Total_N to calculate the weights (but sort first)*/

proc sort data=ps_scores_strata; by stratum; run;

proc sort data=stratum_counts; by stratum; run;

/*Calculate PS Weights*/

data ps_scores_weighted;
    merge ps_scores_strata (in=a) stratum_counts;
    by stratum;
    if a then do;
        
/* Calculate weights for treated and untreated */
    if Generation = 1 then 
         Weight = (Ntotal / &Total_N) / (Nexposed / Ntotal_exposed);
    else if Generation = 2 then 
         Weight = (Ntotal / &Total_N) / (Nunexposed / Ntotal_unexposed);
    
    /* Truncate weights above 10 */
        if Weight > 10 then Weight = 10;
    end;
run;

/*Check out the first 20 weights*/

proc print data=ps_scores_weighted(obs=20);
run;

/*Check summary statistics of the weight*/

proc means data=ps_scores_weighted n mean median min max stddev;
    var Weight;
    class Generation; /* Separate summary by exposure group */
run;

/*Check histogram of the weights*/

proc sgplot data=ps_scores_weighted;
    histogram Weight / group=Generation transparency=0.5;
    xaxis label="Weights";
    yaxis label="Frequency";
run;

proc sgplot data=ps_scores_weighted;
    density Weight / group=Generation type=kernel transparency=0.5;
    keylegend / position=topright; /* Legend for clarity */
    xaxis label="Weights";
    yaxis label="Density";
    title "Weight Distribution by Group (Area Plot)";
run;

/*Check how many weights are >10*/

proc sql;
    select count(*) as WeightsAbove10
    from ps_scores_weighted
    where Weight > 10;
quit;

/*Sort and Merge the dataset*/

proc sort data=ps_scores_weighted; by id; run;
proc sort data=final_event; by id; run;

/*Merge*/

data finally;
    merge ps_scores_weighted (in=a) final_event;
    by id;
    if a;
run;

/*Finally, run the cox model*/

proc phreg data=finally;
    class Generation(ref="2") / param=ref;
    model time*event(0) = Generation / ties=efron;
    strata combined_band; /* Stratify by combined age and calendar time bands */
    weight Weight; /* Apply fine stratification weights */
    hazardratio Generation / cl=wald; /* Include confidence intervals */
run;

/*Get incidence rates, person times and CIs*/

/*Now, calculate your weighted incidence rate*/
proc means data = finally;					/*Calls on the proc means procedure for basic descriptive statistics*/
var event time_years;						/*Specifies the variables we are intetested in to clculate the incidence rate*/
class Generation;							/*Tells SAS to do the cuaculation for each treatment group*/
weight Weight; /* Apply weights */
output out = Inc_rate sum = event py;		/*Creates a new dataset (Irate) to store the results we are about to generate*/
											/*The sum part provides names for the sum of the variables listed above (events and follow-up)*/
quit;

data Inc_rate2;													/*Create a new dataset named 'IR'*/
set Inc_rate;													/*Tell SAS to use data from IRprep*/
Inc_rate2 = event/(py/1000);									/*Calculate the incidence rate*/
LCI=quantile('chisq',0.025, event*2)/((py/100000)*2);			/*Calculate lower confidence interval*/
UCI=quantile('chisq',0.975, (event+1)*2)/((py/100000)*2);run;	/*calculate the upper confidence interval*/
run;

proc print data=Inc_rate2;
run;




/*.......................Sensitivity Analyses....................................*/

/*.......................Sensitivity Analyses....................................*/

/*Varying the grace period to 6 months*/

/*PART 2: Define the Exposure, Lag, Grace periods, Censoring dates etc*/

data cohort1;
 set final_rx;
 by id date;
 if first.id;		/*We select only the first id and assign the dates time 0*/
 rename date=t0;
run;

/*Identify switchers and determine the switch date*/

data switch;
	merge cohort1(in=a) final_rx(where=(date <'31dec2015'd) rename=(Generation=gen));
	by id;
	if a;
	if Generation ne gen;		/*Only the people switch make it into this dataset*/
run;

/*Create a dataset with only the first patient ID and the date of censoring*/

data switch_data;
	set switch;
	by id date;
	if first.id;			/*Keep only first id. Remember you merged with longform final_rx above*/
	keep id date;			/*Keep these 2 variables and drop everything else*/
	rename date=censor1;	
run;

/*Merge 3 datasets. So we now have t0 from 'cohort1' and censoring date (censor1) from 'switch_data' 
added to the long form (final_rx). T0 is fixed for each person while date varies (date of every refill rx)*/

data continuous_data;
	merge cohort1(in=a) final_rx switch_data;
	by id;
	if a;
	if date < censor1 or censor1=.;
run;

/*Apply rxdate + 30 + grace period (365 days)*/

data c1;
	set continuous_data;
	by id date;
	rxdur=date+30+183;				/*rxdur id prescription duration plus grace period*/
	format rxdur lagrxdur date9.;
	lagrxdur=lag(rxdur);			/*lagrxdur is the end date for the previous rx*/
	if first.id then lagrxdur=date;
	retain discon;
	if first.id then discon=0;
	if date > lagrxdur then discon=1;
	if discon=1 and id=lag(id) and discon=lag(discon) then delete; 
run;

/*Calculate the final censoring date (censor2) for each patient. 'Censor1' only identified 
switchers by capturing their first switching date. Censor two accounts for grace periods in both switchers and non-switchers.*/

data c2;
 	set c1;
 	by id date;
 	format censor2 date9.;
 	if discon = 1 then censor2 = lagrxdur + 183; /*For switchers: Added the 365 for the grace period post rx cessation*/
 	else if last.id then censor2 = rxdur + 183;	 /*For those who stop completely: Added the 365 for the grace period post rx cessation*/
 	if last.id;
run;

/*Now, define and merge the events*/

data merged_data;
    merge final_event (in=a) c2(in=b) final_patients(keep=id end); /* Add patients to get 'end' date*/
    by id;
    if b; 
    end_followup = min(eventdate,end,'31dec2015'd,censor2);
    format end date9.;
    if a and t0 < eventdate <= end_followup then event=1; 
	else event = 0;
	time = end_followup - t0 + 1;
	time_years = time/365.25;
	if time < 365.25 then event=0; /*Apply 3 year lag where we don't count events*/
run;


/*PART 3: Estimate PS and deal with confounding*/

/*First, merge the covariates dataset with the current dataset (but sort first)*/

proc sort data=merged_data; by id; run;
proc sort data=final_covariates; by id; run;

/* Step 2: Merge datasets */

data merged_data2;
    merge merged_data (in=a) final_covariates (in=b);
    by id;
    if a; 	/* Keeps only rows that exist in the main dataset*/
run;

/*Extract the year variable*/

data merged_data3;
    set merged_data2;
    t0_year = year(t0); /* Extract the year from the t0 date */
run;

/* Print the first 10 observations to verify */

proc print data=merged_data3(obs=10); 
    var t0 t0_year;
run;

/*We will create 5-year calendar time bands*/

data merged_data_bands;
    set merged_data3;
    
    /* Create calendar time bands */
    if t0_year >= 1995 and t0_year < 2000 then time_band = "1995-1999";
    else if t0_year >= 2000 and t0_year < 2005 then time_band = "2000-2004";
    else if t0_year >= 2005 and t0_year < 2010 then time_band = "2005-2009";
    else time_band = "2010-2015";

    /* Combine age and calendar time bands */
    combined_band = time_band;
run;

/*Print a few rows to confirm*/

proc print data=merged_data_bands(obs=10); 
    var time_band t0_year combined_band;
run;

proc sort data=merged_data_bands; by combined_band; run;

/*Estimate the PS*/

proc logistic data=merged_data_bands;
	
	/*Binary and categorical variables*/
    class obesity(ref="0") smoking(ref="0") male(ref="0") alcohol(ref="0") cancer(ref="0") 
          dementia(ref="0") nsaids(ref="0") aspirin(ref="0") pvd(ref="0") cad(ref="0") revasc(ref="0") 
          nitrates(ref="0") chf(ref="0") cerebrovascular(ref="0") thyroid(ref="0") statins(ref="0") 
          renal(ref="0") digoxin(ref="0") fibrates(ref="0") clopidogrel(ref="0") warfarin(ref="0") 
          opioids(ref="0") af(ref="0") copd(ref="0") asthma(ref="0") dyslipidemia(ref="0") anaemia(ref="0") 
          angina(ref="0") htn(ref="0") hypotension(ref="0") diabetes(ref="0") mi(ref="0") stroke(ref="0") / param=ref;

    /* Include quadratic and cubic terms for continuous variables */
    model Generation(event="1") = age age*age age*age*age visit visit*visit visit*visit*visit
                                  obesity smoking male alcohol cancer dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
    by combined_band;
    output out=ps_scores p=pscore;
run;

/*Print 30 rows to see propensity scores*/

proc print data=ps_scores(obs=30); run;

/*View overlap*/

proc sgplot data=ps_scores;
    histogram pscore / group=Generation transparency=0.5;
run;


/*Now, we will conduct fine stratification weighting*/

/*Create 50 strata*/

proc rank data=ps_scores out=ps_scores_strata groups=50;
    var pscore; /* The propensity score variable */
    ranks stratum; /* Creates the stratum variable (0 to 49) */
run;

/* Adjust stratum to 1-based index for easier interpretation*/

data ps_scores_strata;
    set ps_scores_strata;
    stratum = stratum + 1;
run;

/*Calculate population counts for each stratum*/

proc sql;
    create table stratum_counts as
    select stratum,
           count(*) as Ntotal,
           sum(Generation = 1) as Nexposed, /* Assuming Generation=1 is the exposed group */
           sum(Generation = 2) as Nunexposed, /* Assuming Generation=2 is the unexposed group */
           (select count(*) from ps_scores_strata where Generation = 1) as Ntotal_exposed,
           (select count(*) from ps_scores_strata where Generation = 2) as Ntotal_unexposed
    from ps_scores_strata
    group by stratum;
quit;

/*Calculate weights*/

/*First, calculate the total number of obs (sum(Ntotal)) across all strata and save it as a macro variable*/

proc sql noprint;
    select sum(Ntotal) into :Total_N from stratum_counts;
quit;

/* Check the macro variable value */

%put &=Total_N;

/*Use precomputed &Total_N to calculate the weights (but sort first)*/

proc sort data=ps_scores_strata; by stratum; run;

proc sort data=stratum_counts; by stratum; run;

/*Calculate PS Weights*/

data ps_scores_weighted;
    merge ps_scores_strata (in=a) stratum_counts;
    by stratum;
    if a then do;
        
/* Calculate weights for treated and untreated */
    if Generation = 1 then 
         Weight = (Ntotal / &Total_N) / (Nexposed / Ntotal_exposed);
    else if Generation = 2 then 
         Weight = (Ntotal / &Total_N) / (Nunexposed / Ntotal_unexposed);
    
    /* Truncate weights above 10 */
        if Weight > 10 then Weight = 10;
    end;
run;

/*Check out the first 20 weights*/

proc print data=ps_scores_weighted(obs=20);
run;

/*Check summary statistics of the weight*/

proc means data=ps_scores_weighted n mean median min max stddev;
    var Weight;
    class Generation; /* Separate summary by exposure group */
run;

/*Check histogram of the weights*/

proc sgplot data=ps_scores_weighted;
    histogram Weight / group=Generation transparency=0.5;
    xaxis label="Weights";
    yaxis label="Frequency";
run;

proc sgplot data=ps_scores_weighted;
    density Weight / group=Generation type=kernel transparency=0.5;
    keylegend / position=topright; /* Legend for clarity */
    xaxis label="Weights";
    yaxis label="Density";
    title "Weight Distribution by Group (Area Plot)";
run;

/*Check how many weights are >10*/

proc sql;
    select count(*) as WeightsAbove10
    from ps_scores_weighted
    where Weight > 10;
quit;

/*Sort and Merge the dataset*/

proc sort data=ps_scores_weighted; by id; run;
proc sort data=final_event; by id; run;

/*Merge*/

data finally;
    merge ps_scores_weighted (in=a) final_event;
    by id;
    if a;
run;

/*Finally, run the cox model*/

proc phreg data=finally;
    class Generation(ref="2") / param=ref;
    model time*event(0) = Generation / ties=efron;
    strata combined_band; /* Stratify by combined age and calendar time bands */
    weight Weight; /* Apply fine stratification weights */
    hazardratio Generation / cl=wald; /* Include confidence intervals */
run;

/*Get incidence rates, person times and CIs*/

/*Now, calculate your incidence rate*/
proc means data = finally;					/*Calls on the proc means procedure for basic descriptive statistics*/
var event time_years;						/*Specifies the variables we are intetested in to clculate the incidence rate*/
class Generation;							/*Tells SAS to do the cuaculation for each treatment group*/
weight Weight; 								/* Apply weights */
output out = Inc_rate sum = event py;		/*Creates a new dataset (Irate) to store the results we are about to generate*/
											/*The sum part provides names for the sum of the variables listed above (events and follow-up)*/
quit;

data Inc_rate2;													/*Create a new dataset named 'IR'*/
set Inc_rate;													/*Tell SAS to use data from IRprep*/
Inc_rate2 = event/(py/1000);									/*Calculate the incidence rate*/
LCI=quantile('chisq',0.025, event*2)/((py/100000)*2);			/*Calculate lower confidence interval*/
UCI=quantile('chisq',0.975, (event+1)*2)/((py/100000)*2);run;	/*calculate the upper confidence interval*/
run;

proc print data=Inc_rate2;
run;


/*.......................Sensitivity Analyses....................................*/

/*.......................Sensitivity Analyses....................................*/

/*Varying the grace period to 2 years*/

/*PART 2: Define the Exposure, Lag, Grace periods, Censoring dates etc*/

data cohort1;
 set final_rx;
 by id date;
 if first.id;		/*We select only the first id and assign the dates time 0*/
 rename date=t0;
run;

/*Identify switchers and determine the switch date*/

data switch;
	merge cohort1(in=a) final_rx(where=(date <'31dec2015'd) rename=(Generation=gen));
	by id;
	if a;
	if Generation ne gen;		/*Only the people switch make it into this dataset*/
run;

/*Create a dataset with only the first patient ID and the date of censoring*/

data switch_data;
	set switch;
	by id date;
	if first.id;			/*Keep only first id. Remember you merged with longform final_rx above*/
	keep id date;			/*Keep these 2 variables and drop everything else*/
	rename date=censor1;	
run;

/*Merge 3 datasets. So we now have t0 from 'cohort1' and censoring date (censor1) from 'switch_data' 
added to the long form (final_rx). T0 is fixed for each person while date varies (date of every refill rx)*/

data continuous_data;
	merge cohort1(in=a) final_rx switch_data;
	by id;
	if a;
	if date < censor1 or censor1=.;
run;

/*Apply rxdate + 30 + grace period (365 days)*/

data c1;
	set continuous_data;
	by id date;
	rxdur=date+30+730;				/*rxdur id prescription duration plus grace period*/
	format rxdur lagrxdur date9.;
	lagrxdur=lag(rxdur);			/*lagrxdur is the end date for the previous rx*/
	if first.id then lagrxdur=date;
	retain discon;
	if first.id then discon=0;
	if date > lagrxdur then discon=1;
	if discon=1 and id=lag(id) and discon=lag(discon) then delete; 
run;

/*Calculate the final censoring date (censor2) for each patient. 'Censor1' only identified 
switchers by capturing their first switching date. Censor two accounts for grace periods in both switchers and non-switchers.*/

data c2;
 	set c1;
 	by id date;
 	format censor2 date9.;
 	if discon = 1 then censor2 = lagrxdur + 730; /*For switchers: Added the 365 for the grace period post rx cessation*/
 	else if last.id then censor2 = rxdur + 730;	 /*For those who stop completely: Added the 365 for the grace period post rx cessation*/
 	if last.id;
run;

/*Now, define and merge the events*/

data merged_data;
    merge final_event (in=a) c2(in=b) final_patients(keep=id end); /* Add patients to get 'end' date*/
    by id;
    if b; 
    end_followup = min(eventdate,end,'31dec2015'd,censor2);
    format end date9.;
    if a and t0 < eventdate <= end_followup then event=1; 
	else event = 0;
	time = end_followup - t0 + 1;
	time_years = time/365.25;
	if time < 365.25 then event=0; /*Apply 3 year lag where we don't count events*/
run;


/*PART 3: Estimate PS and deal with confounding*/

/*First, merge the covariates dataset with the current dataset (but sort first)*/

proc sort data=merged_data; by id; run;
proc sort data=final_covariates; by id; run;

/* Step 2: Merge datasets */

data merged_data2;
    merge merged_data (in=a) final_covariates (in=b);
    by id;
    if a; 	/* Keeps only rows that exist in the main dataset*/
run;

/*Extract the year variable*/

data merged_data3;
    set merged_data2;
    t0_year = year(t0); /* Extract the year from the t0 date */
run;

/* Print the first 10 observations to verify */

proc print data=merged_data3(obs=10); 
    var t0 t0_year;
run;

/*We will create 5-year calendar time bands*/

data merged_data_bands;
    set merged_data3;
    
    /* Create calendar time bands */
    if t0_year >= 1995 and t0_year < 2000 then time_band = "1995-1999";
    else if t0_year >= 2000 and t0_year < 2005 then time_band = "2000-2004";
    else if t0_year >= 2005 and t0_year < 2010 then time_band = "2005-2009";
    else time_band = "2010-2015";

    /* Combine age and calendar time bands */
    combined_band = time_band;
run;

/*Print a few rows to confirm*/

proc print data=merged_data_bands(obs=10); 
    var time_band t0_year combined_band;
run;

proc sort data=merged_data_bands; by combined_band; run;

/*Estimate the PS*/

proc logistic data=merged_data_bands;
	
	/*Binary and categorical variables*/
    class obesity(ref="0") smoking(ref="0") male(ref="0") alcohol(ref="0") cancer(ref="0") 
          dementia(ref="0") nsaids(ref="0") aspirin(ref="0") pvd(ref="0") cad(ref="0") revasc(ref="0") 
          nitrates(ref="0") chf(ref="0") cerebrovascular(ref="0") thyroid(ref="0") statins(ref="0") 
          renal(ref="0") digoxin(ref="0") fibrates(ref="0") clopidogrel(ref="0") warfarin(ref="0") 
          opioids(ref="0") af(ref="0") copd(ref="0") asthma(ref="0") dyslipidemia(ref="0") anaemia(ref="0") 
          angina(ref="0") htn(ref="0") hypotension(ref="0") diabetes(ref="0") mi(ref="0") stroke(ref="0") / param=ref;

    /* Include quadratic and cubic terms for continuous variables */
    model Generation(event="1") = age age*age age*age*age visit visit*visit visit*visit*visit
                                  obesity smoking male alcohol cancer dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
    by combined_band;
    output out=ps_scores p=pscore;
run;

/*Print 30 rows to see propensity scores*/

proc print data=ps_scores(obs=30); run;

/*View overlap*/

proc sgplot data=ps_scores;
    histogram pscore / group=Generation transparency=0.5;
run;


/*Now, we will conduct fine stratification weighting*/

/*Create 50 strata*/

proc rank data=ps_scores out=ps_scores_strata groups=50;
    var pscore; /* The propensity score variable */
    ranks stratum; /* Creates the stratum variable (0 to 49) */
run;

/* Adjust stratum to 1-based index for easier interpretation*/

data ps_scores_strata;
    set ps_scores_strata;
    stratum = stratum + 1;
run;

/*Calculate population counts for each stratum*/

proc sql;
    create table stratum_counts as
    select stratum,
           count(*) as Ntotal,
           sum(Generation = 1) as Nexposed, /* Assuming Generation=1 is the exposed group */
           sum(Generation = 2) as Nunexposed, /* Assuming Generation=2 is the unexposed group */
           (select count(*) from ps_scores_strata where Generation = 1) as Ntotal_exposed,
           (select count(*) from ps_scores_strata where Generation = 2) as Ntotal_unexposed
    from ps_scores_strata
    group by stratum;
quit;

/*Calculate weights*/

/*First, calculate the total number of obs (sum(Ntotal)) across all strata and save it as a macro variable*/

proc sql noprint;
    select sum(Ntotal) into :Total_N from stratum_counts;
quit;

/* Check the macro variable value */

%put &=Total_N;

/*Use precomputed &Total_N to calculate the weights (but sort first)*/

proc sort data=ps_scores_strata; by stratum; run;

proc sort data=stratum_counts; by stratum; run;

/*Calculate PS Weights*/

data ps_scores_weighted;
    merge ps_scores_strata (in=a) stratum_counts;
    by stratum;
    if a then do;
        
/* Calculate weights for treated and untreated */
    if Generation = 1 then 
         Weight = (Ntotal / &Total_N) / (Nexposed / Ntotal_exposed);
    else if Generation = 2 then 
         Weight = (Ntotal / &Total_N) / (Nunexposed / Ntotal_unexposed);
    
    /* Truncate weights above 10 */
        if Weight > 10 then Weight = 10;
    end;
run;

/*Check out the first 20 weights*/

proc print data=ps_scores_weighted(obs=20);
run;

/*Check summary statistics of the weight*/

proc means data=ps_scores_weighted n mean median min max stddev;
    var Weight;
    class Generation; /* Separate summary by exposure group */
run;

/*Check histogram of the weights*/

proc sgplot data=ps_scores_weighted;
    histogram Weight / group=Generation transparency=0.5;
    xaxis label="Weights";
    yaxis label="Frequency";
run;

proc sgplot data=ps_scores_weighted;
    density Weight / group=Generation type=kernel transparency=0.5;
    keylegend / position=topright; /* Legend for clarity */
    xaxis label="Weights";
    yaxis label="Density";
    title "Weight Distribution by Group (Area Plot)";
run;

/*Check how many weights are >10*/

proc sql;
    select count(*) as WeightsAbove10
    from ps_scores_weighted
    where Weight > 10;
quit;

/*Sort and Merge the dataset*/

proc sort data=ps_scores_weighted; by id; run;
proc sort data=final_event; by id; run;

/*Merge*/

data finally;
    merge ps_scores_weighted (in=a) final_event;
    by id;
    if a;
run;

/*Finally, run the cox model*/

proc phreg data=finally;
    class Generation(ref="2") / param=ref;
    model time*event(0) = Generation / ties=efron;
    strata combined_band; /* Stratify by combined age and calendar time bands */
    weight Weight; /* Apply fine stratification weights */
    hazardratio Generation / cl=wald; /* Include confidence intervals */
run;

/*Get incidence rates, person times and CIs*/

/*Now, calculate your incidence rate*/
proc means data = finally;					/*Calls on the proc means procedure for basic descriptive statistics*/
var event time_years;						/*Specifies the variables we are intetested in to clculate the incidence rate*/
class Generation;							/*Tells SAS to do the cuaculation for each treatment group*/
weight Weight; 								/* Apply weights */
output out = Inc_rate sum = event py;		/*Creates a new dataset (Irate) to store the results we are about to generate*/
											/*The sum part provides names for the sum of the variables listed above (events and follow-up)*/
quit;

data Inc_rate2;													/*Create a new dataset named 'IR'*/
set Inc_rate;													/*Tell SAS to use data from IRprep*/
Inc_rate2 = event/(py/1000);									/*Calculate the incidence rate*/
LCI=quantile('chisq',0.025, event*2)/((py/100000)*2);			/*Calculate lower confidence interval*/
UCI=quantile('chisq',0.975, (event+1)*2)/((py/100000)*2);run;	/*calculate the upper confidence interval*/
run;

proc print data=Inc_rate2;
run;


/*..............................Sensitivity ITT Analysis......................................*/

/*..............................Sensitivity ITT Analysis......................................*/

/*PART 2: Define the Exposure, Lag, Grace periods, Censoring dates etc*/

data cohort1;
 set final_rx;
 by id date;
 if first.id;		/*We select only the first id and assign the dates time 0*/
 rename date=t0;
run;

data mergedcohort;										/*We are merging both cohorts to create a new cohort named 'mergedcohort'*/
merge final_event(in = a) cohort1 (in = b); 			/*Create temporary variables 'a' and 'b' that tell us whether an onservation from each dataset is included in the merge*/
by ID; 													/*Tell SAS we are sorting by ID*/
if b;													/*We are keeping everyone in the prescription dataset (b) even if they did not have an event because we are interested in those who did and did not have the outcome*/
end_followup = min('31dec2015'd, end, eventdate);		/*Create a new variable end whose value will be the minimum of these 3 dates: one year from time 0, end of our follow up or event date - whichever occurs first*/
format end date9. end_followup date9.;					/*Change the date format or style*/
if a and t0 < eventdate <= end_followup then event=1;	/*Create a condition if event occured after trestment commenced and before study ended, event = 1*/
else event = 0;											/*Or else, event = 0*/
time = end_followup - t0+1; 							/*We calculate follow up time by subtracting the start date from the end date minus 1. Note our calculation is all in days.*/
follow = time/365.25;
if time < 365 then event = 0; 							/*We are converting time to years. (.25 to account for leap yeears)*/
run;

/*PART 3: Estimate PS and deal with confounding*/

/*First, merge the covariates dataset with the current dataset (but sort first)*/

proc sort data=merged_data; by id; run;
proc sort data=final_covariates; by id; run;

/* Step 2: Merge datasets */

data merged_data2;
    merge merged_data (in=a) final_covariates (in=b);
    by id;
    if a; 	/* Keeps only rows that exist in the main dataset*/
run;

/*Extract the year variable*/

data merged_data3;
    set merged_data2;
    t0_year = year(t0); /* Extract the year from the t0 date */
run;

/* Print the first 10 observations to verify */

proc print data=merged_data3(obs=10); 
    var t0 t0_year;
run;

/*We will create 5-year calendar time bands*/

data merged_data_bands;
    set merged_data3;
    
    /* Create calendar time bands */
    if t0_year >= 1995 and t0_year < 2000 then time_band = "1995-1999";
    else if t0_year >= 2000 and t0_year < 2005 then time_band = "2000-2004";
    else if t0_year >= 2005 and t0_year < 2010 then time_band = "2005-2009";
    else time_band = "2010-2015";

    /* Combine age and calendar time bands */
    combined_band = time_band;
run;

/*Print a few rows to confirm*/

proc print data=merged_data_bands(obs=10); 
    var time_band t0_year combined_band;
run;

proc sort data=merged_data_bands; by combined_band; run;

/*Estimate the PS*/

proc logistic data=merged_data_bands;
	
	/*Binary and categorical variables*/
    class obesity(ref="0") smoking(ref="0") male(ref="0") alcohol(ref="0") cancer(ref="0") 
          dementia(ref="0") nsaids(ref="0") aspirin(ref="0") pvd(ref="0") cad(ref="0") revasc(ref="0") 
          nitrates(ref="0") chf(ref="0") cerebrovascular(ref="0") thyroid(ref="0") statins(ref="0") 
          renal(ref="0") digoxin(ref="0") fibrates(ref="0") clopidogrel(ref="0") warfarin(ref="0") 
          opioids(ref="0") af(ref="0") copd(ref="0") asthma(ref="0") dyslipidemia(ref="0") anaemia(ref="0") 
          angina(ref="0") htn(ref="0") hypotension(ref="0") diabetes(ref="0") mi(ref="0") stroke(ref="0") / param=ref;

    /* Include quadratic and cubic terms for continuous variables */
    model Generation(event="1") = age age*age age*age*age visit visit*visit visit*visit*visit
                                  obesity smoking male alcohol cancer dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
    by combined_band;
    output out=ps_scores p=pscore;
run;

/*Print 30 rows to see propensity scores*/

proc print data=ps_scores(obs=30); run;

/*View overlap*/

proc sgplot data=ps_scores;
    histogram pscore / group=Generation transparency=0.5;
run;


/*Now, we will conduct fine stratification weighting*/

/*Create 50 strata*/

proc rank data=ps_scores out=ps_scores_strata groups=50;
    var pscore; /* The propensity score variable */
    ranks stratum; /* Creates the stratum variable (0 to 49) */
run;

/* Adjust stratum to 1-based index for easier interpretation*/

data ps_scores_strata;
    set ps_scores_strata;
    stratum = stratum + 1;
run;

/*Calculate population counts for each stratum*/

proc sql;
    create table stratum_counts as
    select stratum,
           count(*) as Ntotal,
           sum(Generation = 1) as Nexposed, /* Assuming Generation=1 is the exposed group */
           sum(Generation = 2) as Nunexposed, /* Assuming Generation=2 is the unexposed group */
           (select count(*) from ps_scores_strata where Generation = 1) as Ntotal_exposed,
           (select count(*) from ps_scores_strata where Generation = 2) as Ntotal_unexposed
    from ps_scores_strata
    group by stratum;
quit;

/*Calculate weights*/

/*First, calculate the total number of obs (sum(Ntotal)) across all strata and save it as a macro variable*/

proc sql noprint;
    select sum(Ntotal) into :Total_N from stratum_counts;
quit;

/* Check the macro variable value */

%put &=Total_N;

/*Use precomputed &Total_N to calculate the weights (but sort first)*/

proc sort data=ps_scores_strata; by stratum; run;

proc sort data=stratum_counts; by stratum; run;

/*Calculate PS Weights*/

data ps_scores_weighted;
    merge ps_scores_strata (in=a) stratum_counts;
    by stratum;
    if a then do;
        
/* Calculate weights for treated and untreated */
    if Generation = 1 then 
         Weight = (Ntotal / &Total_N) / (Nexposed / Ntotal_exposed);
    else if Generation = 2 then 
         Weight = (Ntotal / &Total_N) / (Nunexposed / Ntotal_unexposed);
    
    /* Truncate weights above 10 */
        if Weight > 10 then Weight = 10;
    end;
run;

/*Check out the first 20 weights*/

proc print data=ps_scores_weighted(obs=20);
run;

/*Check summary statistics of the weight*/

proc means data=ps_scores_weighted n mean median min max stddev;
    var Weight;
    class Generation; /* Separate summary by exposure group */
run;

/*Check histogram of the weights*/

proc sgplot data=ps_scores_weighted;
    histogram Weight / group=Generation transparency=0.5;
    xaxis label="Weights";
    yaxis label="Frequency";
run;

proc sgplot data=ps_scores_weighted;
    density Weight / group=Generation type=kernel transparency=0.5;
    keylegend / position=topright; /* Legend for clarity */
    xaxis label="Weights";
    yaxis label="Density";
    title "Weight Distribution by Group (Area Plot)";
run;

/*Check how many weights are >10*/

proc sql;
    select count(*) as WeightsAbove10
    from ps_scores_weighted
    where Weight > 10;
quit;

/*Sort and Merge the dataset*/

proc sort data=ps_scores_weighted; by id; run;
proc sort data=final_event; by id; run;

/*Merge*/

data finally;
    merge ps_scores_weighted (in=a) final_event;
    by id;
    if a;
run;

/*Finally, run the cox model*/

proc phreg data=finally;
    class Generation(ref="2") / param=ref;
    model time*event(0) = Generation / ties=efron;
    strata combined_band; /* Stratify by combined age and calendar time bands */
    weight Weight; /* Apply fine stratification weights */
    hazardratio Generation / cl=wald; /* Include confidence intervals */
run;

/*Get incidence rates, person times and CIs*/

/*Now, calculate your incidence rate*/
proc means data = finally;					/*Calls on the proc means procedure for basic descriptive statistics*/
var event time_years;						/*Specifies the variables we are intetested in to clculate the incidence rate*/
class Generation;							/*Tells SAS to do the cuaculation for each treatment group*/
weight Weight; 								/* Apply weights */
output out = Inc_rate sum = event py;		/*Creates a new dataset (Irate) to store the results we are about to generate*/
											/*The sum part provides names for the sum of the variables listed above (events and follow-up)*/
quit;

data Inc_rate2;													/*Create a new dataset named 'IR'*/
set Inc_rate;													/*Tell SAS to use data from IRprep*/
Inc_rate2 = event/(py/1000);									/*Calculate the incidence rate*/
LCI=quantile('chisq',0.025, event*2)/((py/100000)*2);			/*Calculate lower confidence interval*/
UCI=quantile('chisq',0.975, (event+1)*2)/((py/100000)*2);run;	/*calculate the upper confidence interval*/
run;

proc print data=Inc_rate2;
run;


/*..............................Sensitivity Modified ITT Analysis......................................*/

/*..............................Sensitivity Modified ITT Analysis......................................*/

/*PART 2: Define the Exposure, Lag, Grace periods, Censoring dates etc*/

data cohort1;
 set final_rx;
 by id date;
 if first.id;		/*We select only the first id and assign the dates time 0*/
 rename date=t0;
run;

/*Identify switchers and determine the switch date*/

data switch;
	merge cohort1(in=a) final_rx(where=(date <'31dec2015'd) rename=(Generation=gen));
	by id;
	if a;
	if Generation ne gen;		/*Only the people switch make it into this dataset*/
run;

/*Create a dataset with only the first patient ID and the date of censoring*/

data switch_data;
	set switch;
	by id date;
	if first.id;			/*Keep only first id. Remember you merged with longform final_rx above*/
	keep id date;			/*Keep these 2 variables and drop everything else*/
	rename date=censor1;	
run;

/*Merge 3 datasets. So we now have t0 from 'cohort1' and censoring date (censor1) from 'switch_data' 
added to the long form (final_rx). T0 is fixed for each person while date varies (date of every refill rx)*/

data continuous_data;
	merge cohort1(in=a) final_rx switch_data;
	by id;
	if a;
	if date < censor1 or censor1=.;
run;

/*Apply rxdate + 30 + grace period (365 days)*/

data c1;
	set continuous_data;
	by id date;
	rxdur=date+30+365;				/*rxdur id prescription duration plus grace period*/
	format rxdur lagrxdur date9.;
	lagrxdur=lag(rxdur);			/*lagrxdur is the end date for the previous rx*/
	if first.id then lagrxdur=date;
	retain discon;
	if first.id then discon=0;
	if date > lagrxdur then discon=1;
	if discon=1 and id=lag(id) and discon=lag(discon) then delete; 
run;

/*Calculate the final censoring date (censor2) for each patient. 'Censor1' only identified 
switchers by capturing their first switching date. Censor two accounts for grace periods in both switchers and non-switchers.*/

data c2;
 	set c1;
 	by id date;
 	format censor2 date9.;
 	if discon = 1 then censor2 = lagrxdur + 365; /*For switchers: Added the 365 for the grace period post rx cessation*/
 	else if last.id then censor2 =.;	 /*For those who stop completely: Added the 365 for the grace period post rx cessation*/
 	if last.id;
run;

/*Now, define and merge the events*/

data merged_data;
    merge final_event (in=a) c2(in=b) final_patients(keep=id end); /* Add patients to get 'end' date*/
    by id;
    if b; 
    if missing(censor2) then end_followup = min(eventdate, end, '31dec2015'd); /* Stoppers */
    else end_followup = min(eventdate, end, '31dec2015'd, censor2); /* Switchers and others */
    format end date9.;
    if a and t0 < eventdate <= end_followup then event=1; 
	else event = 0;
	time = end_followup - t0 + 1;
	time_years = time/365.25;
	if time < 365.25 then event=0; /*Apply 1 year lag where we don't count events*/
run;


/*PART 3: Estimate PS and deal with confounding*/

/*First, merge the covariates dataset with the current dataset (but sort first)*/

proc sort data=merged_data; by id; run;
proc sort data=final_covariates; by id; run;

/* Step 2: Merge datasets */

data merged_data2;
    merge merged_data (in=a) final_covariates (in=b);
    by id;
    if a; 	/* Keeps only rows that exist in the main dataset*/
run;

/*Extract the year variable*/

data merged_data3;
    set merged_data2;
    t0_year = year(t0); /* Extract the year from the t0 date */
run;

/* Print the first 10 observations to verify */

proc print data=merged_data3(obs=10); 
    var t0 t0_year;
run;

/*We will create 5-year calendar time bands*/

data merged_data_bands;
    set merged_data3;
    
    /* Create calendar time bands */
    if t0_year >= 1995 and t0_year < 2000 then time_band = "1995-1999";
    else if t0_year >= 2000 and t0_year < 2005 then time_band = "2000-2004";
    else if t0_year >= 2005 and t0_year < 2010 then time_band = "2005-2009";
    else time_band = "2010-2015";

    /* Combine age and calendar time bands */
    combined_band = time_band;
run;

/*Print a few rows to confirm*/

proc print data=merged_data_bands(obs=10); 
    var time_band t0_year combined_band;
run;

proc sort data=merged_data_bands; by combined_band; run;

/*Estimate the PS*/

proc logistic data=merged_data_bands;
	
	/*Binary and categorical variables*/
    class obesity(ref="0") smoking(ref="0") male(ref="0") alcohol(ref="0") cancer(ref="0") 
          dementia(ref="0") nsaids(ref="0") aspirin(ref="0") pvd(ref="0") cad(ref="0") revasc(ref="0") 
          nitrates(ref="0") chf(ref="0") cerebrovascular(ref="0") thyroid(ref="0") statins(ref="0") 
          renal(ref="0") digoxin(ref="0") fibrates(ref="0") clopidogrel(ref="0") warfarin(ref="0") 
          opioids(ref="0") af(ref="0") copd(ref="0") asthma(ref="0") dyslipidemia(ref="0") anaemia(ref="0") 
          angina(ref="0") htn(ref="0") hypotension(ref="0") diabetes(ref="0") mi(ref="0") stroke(ref="0") / param=ref;

    /* Include quadratic and cubic terms for continuous variables */
    model Generation(event="1") = age age*age age*age*age visit visit*visit visit*visit*visit
                                  obesity smoking male alcohol cancer dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
    by combined_band;
    output out=ps_scores p=pscore;
run;

/*Print 30 rows to see propensity scores*/

proc print data=ps_scores(obs=30); run;

/*View overlap*/

proc sgplot data=ps_scores;
    histogram pscore / group=Generation transparency=0.5;
run;


/*Now, we will conduct fine stratification weighting*/

/*Create 50 strata*/

proc rank data=ps_scores out=ps_scores_strata groups=50;
    var pscore; /* The propensity score variable */
    ranks stratum; /* Creates the stratum variable (0 to 49) */
run;

/* Adjust stratum to 1-based index for easier interpretation*/

data ps_scores_strata;
    set ps_scores_strata;
    stratum = stratum + 1;
run;

/*Calculate population counts for each stratum*/

proc sql;
    create table stratum_counts as
    select stratum,
           count(*) as Ntotal,
           sum(Generation = 1) as Nexposed, /* Assuming Generation=1 is the exposed group */
           sum(Generation = 2) as Nunexposed, /* Assuming Generation=2 is the unexposed group */
           (select count(*) from ps_scores_strata where Generation = 1) as Ntotal_exposed,
           (select count(*) from ps_scores_strata where Generation = 2) as Ntotal_unexposed
    from ps_scores_strata
    group by stratum;
quit;

/*Calculate weights*/

/*First, calculate the total number of obs (sum(Ntotal)) across all strata and save it as a macro variable*/

proc sql noprint;
    select sum(Ntotal) into :Total_N from stratum_counts;
quit;

/* Check the macro variable value */

%put &=Total_N;

/*Use precomputed &Total_N to calculate the weights (but sort first)*/

proc sort data=ps_scores_strata; by stratum; run;

proc sort data=stratum_counts; by stratum; run;

/*Calculate PS Weights*/

data ps_scores_weighted;
    merge ps_scores_strata (in=a) stratum_counts;
    by stratum;
    if a then do;
        
/* Calculate weights for treated and untreated */
    if Generation = 1 then 
         Weight = (Ntotal / &Total_N) / (Nexposed / Ntotal_exposed);
    else if Generation = 2 then 
         Weight = (Ntotal / &Total_N) / (Nunexposed / Ntotal_unexposed);
    
    /* Truncate weights above 10 */
        if Weight > 10 then Weight = 10;
    end;
run;

/*Check out the first 20 weights*/

proc print data=ps_scores_weighted(obs=20);
run;

/*Check summary statistics of the weight*/

proc means data=ps_scores_weighted n mean median min max stddev;
    var Weight;
    class Generation; /* Separate summary by exposure group */
run;

/*Check histogram of the weights*/

proc sgplot data=ps_scores_weighted;
    histogram Weight / group=Generation transparency=0.5;
    xaxis label="Weights";
    yaxis label="Frequency";
run;

proc sgplot data=ps_scores_weighted;
    density Weight / group=Generation type=kernel transparency=0.5;
    keylegend / position=topright; /* Legend for clarity */
    xaxis label="Weights";
    yaxis label="Density";
    title "Weight Distribution by Group (Area Plot)";
run;

/*Check how many weights are >10*/

proc sql;
    select count(*) as WeightsAbove10
    from ps_scores_weighted
    where Weight > 10;
quit;

/*Sort and Merge the dataset*/

proc sort data=ps_scores_weighted; by id; run;
proc sort data=final_event; by id; run;

/*Merge*/

data finally;
    merge ps_scores_weighted (in=a) final_event;
    by id;
    if a;
run;

/*Finally, run the cox model*/

proc phreg data=finally;
    class Generation(ref="2") / param=ref;
    model time*event(0) = Generation / ties=efron;
    strata combined_band; /* Stratify by combined age and calendar time bands */
    weight Weight; /* Apply fine stratification weights */
    hazardratio Generation / cl=wald; /* Include confidence intervals */
run;

/*Get incidence rates, person times and CIs (these are unweighted and we didnt use them)*/

/*............ Weighted incidence rate (what we used)..........*/

proc means data=finally sum; 
var event time_years;
class Generation;
weight Weight; /* Apply weights */
output out=Inc_rate sum=event py;
run;

data Inc_rate2;													/*Create a new dataset named 'IR'*/
set Inc_rate;													/*Tell SAS to use data from IRprep*/
Inc_rate2 = event/(py/1000);									/*Calculate the incidence rate*/
LCI=quantile('chisq',0.025, event*2)/((py/100000)*2);			/*Calculate lower confidence interval*/
UCI=quantile('chisq',0.975, (event+1)*2)/((py/100000)*2);run;	/*calculate the upper confidence interval*/
run;

proc print data=Inc_rate2;
run;


/*....................After thought sensitivity_Truncation at 99th percentile........................*/

/*..................................After thought sensitivity.......................................*/

/*PART 3: Estimate PS and deal with confounding*/

/*First, merge the covariates dataset with the current dataset (but sort first)*/

proc sort data=merged_data; by id; run;
proc sort data=final_covariates; by id; run;

/* Step 2: Merge datasets */

data merged_data2;
    merge merged_data (in=a) final_covariates (in=b);
    by id;
    if a; 	/* Keeps only rows that exist in the main dataset*/
run;

/*Extract the year variable*/

data merged_data3;
    set merged_data2;
    t0_year = year(t0); /* Extract the year from the t0 date */
run;

/* Print the first 10 observations to verify */

proc print data=merged_data3(obs=10); 
    var t0 t0_year;
run;

/*We will create 5-year calendar time bands*/

data merged_data_bands;
    set merged_data3;
    
    /* Create calendar time bands */
    if t0_year >= 1995 and t0_year < 2000 then time_band = "1995-1999";
    else if t0_year >= 2000 and t0_year < 2005 then time_band = "2000-2004";
    else if t0_year >= 2005 and t0_year < 2010 then time_band = "2005-2009";
    else time_band = "2010-2015";

    /* Combine age and calendar time bands */
    combined_band = time_band;
run;

/*Print a few rows to confirm*/

proc print data=merged_data_bands(obs=10); 
    var time_band t0_year combined_band;
run;

proc sort data=merged_data_bands; by combined_band; run;

/*Estimate the PS*/

proc logistic data=merged_data_bands;
	
	/*Binary and categorical variables*/
    class obesity(ref="0") smoking(ref="0") male(ref="0") alcohol(ref="0") cancer(ref="0") 
          dementia(ref="0") nsaids(ref="0") aspirin(ref="0") pvd(ref="0") cad(ref="0") revasc(ref="0") 
          nitrates(ref="0") chf(ref="0") cerebrovascular(ref="0") thyroid(ref="0") statins(ref="0") 
          renal(ref="0") digoxin(ref="0") fibrates(ref="0") clopidogrel(ref="0") warfarin(ref="0") 
          opioids(ref="0") af(ref="0") copd(ref="0") asthma(ref="0") dyslipidemia(ref="0") anaemia(ref="0") 
          angina(ref="0") htn(ref="0") hypotension(ref="0") diabetes(ref="0") mi(ref="0") stroke(ref="0") / param=ref;

    /* Include quadratic and cubic terms for continuous variables */
    model Generation(event="1") = age age*age age*age*age visit visit*visit visit*visit*visit
                                  obesity smoking male alcohol cancer dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
    by combined_band;
    output out=ps_scores p=pscore;
run;

/*Print 30 rows to see propensity scores*/

proc print data=ps_scores(obs=30); run;

/*View overlap*/

proc sgplot data=ps_scores;
    histogram pscore / group=Generation transparency=0.5;
run;


/*Now, we will conduct fine stratification weighting*/

/*Create 50 strata*/

proc rank data=ps_scores out=ps_scores_strata groups=50;
    var pscore; /* The propensity score variable */
    ranks stratum; /* Creates the stratum variable (0 to 49) */
run;

/* Adjust stratum to 1-based index for easier interpretation*/

data ps_scores_strata;
    set ps_scores_strata;
    stratum = stratum + 1;
run;

/*Calculate population counts for each stratum*/

proc sql;
    create table stratum_counts as
    select stratum,
           count(*) as Ntotal,
           sum(Generation = 1) as Nexposed, /* Assuming Generation=1 is the exposed group */
           sum(Generation = 2) as Nunexposed, /* Assuming Generation=2 is the unexposed group */
           (select count(*) from ps_scores_strata where Generation = 1) as Ntotal_exposed,
           (select count(*) from ps_scores_strata where Generation = 2) as Ntotal_unexposed
    from ps_scores_strata
    group by stratum;
quit;

/*Calculate weights*/

/*First, calculate the total number of obs (sum(Ntotal)) across all strata and save it as a macro variable*/

proc sql noprint;
    select sum(Ntotal) into :Total_N from stratum_counts;
quit;

/* Check the macro variable value */

%put &=Total_N;

/*Use precomputed &Total_N to calculate the weights (but sort first)*/

proc sort data=ps_scores_strata; by stratum; run;

proc sort data=stratum_counts; by stratum; run;

/*Now, calculate the weghts*/

/* Step 1: Calculate the 99th Percentile of Weights */
proc univariate data=ps_scores_weighted noprint;
    var Weight;
    output out=pctl_results pctlpts=99 pctlpre=Weight_;
run;

/* Step 2: Calculate and Truncate Weights at the 99th Percentile */
data ps_scores_weighted_truncated;
    merge ps_scores_strata (in=a) stratum_counts;
    by stratum;

    if a then do;

        /* Calculate weights for treated and untreated */
        if Generation = 1 then 
             Weight = (Ntotal / &Total_N) / (Nexposed / Ntotal_exposed);
        else if Generation = 2 then 
             Weight = (Ntotal / &Total_N) / (Nunexposed / Ntotal_unexposed);

        /* Merge the 99th percentile value from pctl_results */
        if _N_ = 1 then set pctl_results;

        /* Truncate weights above the 99th percentile */
        if Weight > Weight_99 then Weight = Weight_99;
    end;
run;


/*Check out the first 20 weights*/

proc print data=ps_scores_weighted_truncated(obs=20);
run;

/*Check summary statistics of the weight*/

proc means data=ps_scores_weighted_truncated n mean median min max stddev;
    var Weight;
    class Generation; /* Separate summary by exposure group */
run;

/*Check histogram of the weights*/

proc sgplot data=ps_scores_weighted_truncated;
    histogram Weight / group=Generation transparency=0.5;
    xaxis label="Weights";
    yaxis label="Frequency";
run;

proc sgplot data=ps_scores_weighted_truncated;
    density Weight / group=Generation type=kernel transparency=0.5;
    keylegend / position=topright; /* Legend for clarity */
    xaxis label="Weights";
    yaxis label="Density";
    title "Weight Distribution by Group (Area Plot)";
run;

/*Check how many weights are >10*/

proc sql;
    select count(*) as WeightsAbove10
    from ps_scores_weighted_truncated
    where Weight > 10;
quit;

/*Sort and Merge the dataset*/

proc sort data=ps_scores_weighted_truncated; by id; run;
proc sort data=final_event; by id; run;

/*Merge*/

data finally;
    merge ps_scores_weighted_truncated (in=a) final_event;
    by id;
    if a;
run;

/*Finally, run the cox model*/

proc phreg data=finally;
    class Generation(ref="2") / param=ref;
    model time*event(0) = Generation / ties=efron;
    strata combined_band; /* Stratify by combined age and calendar time bands */
    weight Weight; /* Apply fine stratification weights */
    hazardratio Generation / cl=wald; /* Include confidence intervals */
run;

/*Get incidence rates, person times and CIs (these are unweighted and we didnt use them)*/

/*............ Weighted incidence rate (what we used)..........*/

proc means data=finally sum; 
var event time_years;
class Generation;
weight Weight; /* Apply weights */
output out=Inc_rate sum=event py;
run;

data Inc_rate2;													/*Create a new dataset named 'IR'*/
set Inc_rate;													/*Tell SAS to use data from IRprep*/
Inc_rate2 = event/(py/1000);									/*Calculate the incidence rate*/
LCI=quantile('chisq',0.025, event*2)/((py/100000)*2);			/*Calculate lower confidence interval*/
UCI=quantile('chisq',0.975, (event+1)*2)/((py/100000)*2);run;	/*calculate the upper confidence interval*/
run;

proc print data=Inc_rate2;
run;

/*Some diagnostics since event numbers are low*/

/*Total number of events that occured after censoring = 117 events */

proc sql;
    create table events_after_censor as
    select * from final_event as fe
    left join finally as f
    on fe.id = f.id
    where fe.eventdate > f.end_followup;
quit;

proc sql;
    select count(*) as events_after_censor from events_after_censor;
quit;

/*Total number of events that occured within the one year lag = 2884 events */

proc sql;
    create table first_year_events as
    select * from final_event as fe
    left join finally as f
    on fe.id = f.id
    where f.time < 365.25;
quit;

proc sql;
    select count(*) as excluded_first_year_events from first_year_events;
quit;


/*....................................Bootstrapped Condience Intervals............................*/


%let n_boot = 1000; /* Number of bootstrap iterations */

/* Step 1: Generate Bootstrap Samples */
proc surveyselect data=finally out=boot_samples method=urs samprate=1 outhits reps=&n_boot seed=12345;
    id id Generation combined_band Weight time event; /* Include all required variables */
run;

/* Step 2: Fit the Cox Model for Each Bootstrap Sample */
proc phreg data=boot_samples noprint outest=boot_estimates;
    by Replicate; /* Fit model separately for each bootstrap sample */
    class Generation(ref="2") / param=ref; /* Define categorical variable */
    strata combined_band; /* Stratify by combined bands if applicable */
    weight Weight; /* Apply fine stratification weights */
    model time*event(0) = Generation / ties=efron; /* Cox model */
run;

/* Step 3: Exponentiate to Get Hazard Ratios */
data boot_estimates;
    set boot_estimates;
    HR = exp(Generation1); /* Convert log HR to HR */
run;

/* Step 4: Compute Bootstrapped Confidence Intervals */
proc univariate data=boot_estimates noprint;
    var HR; /* Use the bootstrap HR values */
    output out=boot_ci pctlpts=2.5 97.5 pctlpre=ci_; /* 2.5th and 97.5th percentiles */
run;

/* Step 5: Display Bootstrapped Confidence Intervals */
proc print data=boot_ci;
    title "Bootstrapped Confidence Intervals for Hazard Ratios";
run;


/*...............................Table 1...............................................*/


/*...........Before weighting.....................*/

proc sort data=finally out = c;by generation;run;
ods output basicmeasures=m1;

proc univariate data=c;
var age visit;
by generation;
*weight weight;
quit;

data m2; set m1;
where locmeasure='Mean';
variable=varname;
res=strip(put(locvalue,5.1)) || " (" || strip(put(varvalue,5.1)) || ")";
m=locvalue;
s=varvalue;
keep variable generation res m s;
run;

ods output onewayfreqs=r1;
proc freq data=c;
tables obesity smoking male alcohol cancer dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
by generation;
*weight weight;
run;

 

data r2;
set r1(drop=f_:);
res=strip(put(frequency,comma8.)) || " ( " || strip(put(percent,5.1)) || ")";
variable=scan(table,2);
value=min(obesity, smoking, male, alcohol, cancer, dementia, NSAIDs, aspirin, PVD, CAD, revasc, 
nitrates, CHF, cerebrovascular, thyroid, statins, renal, digoxin, fibrates, clopidogrel, warfarin, 
opioids, af, COPD, asthma, dyslipidemia, anaemia, angina, htn, hypotension, diabetes, mi, stroke);
if variable not in ('smoking', 'obesity') and value = 0 then delete;
run;
 
data r3;
set r2 m2;
run;

proc sort data=r3;
by variable value;
run;

/*standardized diff*/

data r4;
merge r3(where=(generation=1) rename=(percent=p1 frequency=f1 m=m1 s=s1 res=exposed))
        r3(where=(generation=2) rename=(percent=p2 frequency=f2 m=m2 s=s2 res=reference));
by variable value;
if f1=. then f1=0;
if f2=. then f2=0;
if variable in('age' 'visit') then stddiff=abs((m1-m2)/sqrt((s1**2+s2**2)/2));
else if f1>0 and f2>0 and p1<100 and p2<100 then stddiff=abs((p1-p2)/sqrt((p1*(100-p1)+p2*(100-p2))/2));
else stddiff=.;
keep variable exposed reference value stddiff;
run;


/*........................................After weighting..................................*/



proc sort data=finally out = c;by generation;run;
ods output basicmeasures=m1;

proc univariate data=c;
var age visit;
by generation;
weight weight;
quit;

data m2; set m1;
where locmeasure='Mean';
variable=varname;
res=strip(put(locvalue,5.1)) || " (" || strip(put(varvalue,5.1)) || ")";
m=locvalue;
s=varvalue;
keep variable generation res m s;
run;

ods output onewayfreqs=r1;
proc freq data=c;
tables obesity smoking male alcohol cancer dementia nsaids aspirin 
                                  pvd cad revasc nitrates chf cerebrovascular thyroid statins 
                                  renal digoxin fibrates clopidogrel warfarin opioids af copd 
                                  asthma dyslipidemia anaemia angina htn hypotension diabetes mi stroke;
by generation;
weight weight;
run;

 

data r2;
set r1(drop=f_:);
res=strip(put(frequency,comma8.)) || " ( " || strip(put(percent,5.1)) || ")";
variable=scan(table,2);
value=min(obesity, smoking, male, alcohol, cancer, dementia, NSAIDs, aspirin, PVD, CAD, revasc, 
nitrates, CHF, cerebrovascular, thyroid, statins, renal, digoxin, fibrates, clopidogrel, warfarin, 
opioids, af, COPD, asthma, dyslipidemia, anaemia, angina, htn, hypotension, diabetes, mi, stroke);
if variable not in ('smoking', 'obesity') and value = 0 then delete;
run;
 
data r3;
set r2 m2;
run;

proc sort data=r3;
by variable value;
run;

/*standardized diff*/

data r4;
merge r3(where=(generation=1) rename=(percent=p1 frequency=f1 m=m1 s=s1 res=exposed))
        r3(where=(generation=2) rename=(percent=p2 frequency=f2 m=m2 s=s2 res=reference));
by variable value;
if f1=. then f1=0;
if f2=. then f2=0;
if variable in('age' 'visit') then stddiff=abs((m1-m2)/sqrt((s1**2+s2**2)/2));
else if f1>0 and f2>0 and p1<100 and p2<100 then stddiff=abs((p1-p2)/sqrt((p1*(100-p1)+p2*(100-p2))/2));
else stddiff=.;
keep variable exposed reference value stddiff;
run;


/*.................................Cumulative Incidence curves..................................*/


/*use proclifetest, Include strata statement (generation), use 1- survival probability, */

/*Afterbthe proclife code, oitput the dataset and recreate one (using 1-prbability)
then plot.*/ 

/* Step 1: Run PROC LIFETEST to Output Survival Estimates */
proc lifetest data=finally outsurv=test method=km plots=none;
    time time*event(0); /* Define time-to-event and censoring */
    weight weight; /* Apply weights if applicable */
    strata generation / test=logrank; /* Stratify by generation */
run;

/* Step 2: Define a Custom Format for Treatment Labels */
proc format;
    value Treatment_fmt
        1 = "PEAP"
        2 = "PSAP";
run;

/* Step 3: Prepare Data for Plotting */
data test_cif;
    set test_cif; /* Assuming 'test_cif' is the dataset with calculated CIF values */
    Treatment = generation; /* Rename 'generation' to 'treatment' for clarity */
    cumulative_incidence = 1 - survival; /* Compute cumulative incidence */
run;

/* Step 4: Plot Cumulative Incidence Curves */
proc sgplot data=test_cif;
    series x=time_years y=cumulative_incidence / group=Treatment;
    format Treatment Treatment_fmt.; /* Apply the custom format */
    xaxis label="Time since treatment initiation, years" values=(0 to 10 by 1); /* Custom tick marks for each year */
    yaxis label="Cumulative Incidence of Breast Cancer";
    title "Weighted Cumulative Incidence Curves of Breast Cancer for PEAP versus PSAP Use";
run;


