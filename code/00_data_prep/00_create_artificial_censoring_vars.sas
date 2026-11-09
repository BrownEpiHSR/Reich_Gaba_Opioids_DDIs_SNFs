/**************************************************************************
Purpose: Generate the date of artificial censoring for each
clone in the cohort.

The treatment strategies for this study are as follows:

	1. Opioids-only strategy (assign = 0)
	   Administered any dose of opioids (without gabapentinoids) within 3 days of SNF admission.

	2. Gabapentinoids + opioids strategy (assign = 1)
	   Administered any dose of opioids and gabapentinoids within 3 days of SNF admission.

	- Clones are artificially censored as soon as they deviate from their assigned treatment
      strategy.
		- For the opioids-only strategy (assign = 0), the clone is censored as soon as they
		  receive a gabapentinoid (in the grace period), or the clone is censored on the 
		  last day of the grace period if they are not administered an opioid. 
		- For the gabapentinoids + opioids strategy (assign = 1), the clone is censored on 
          the last day of the grace period if they are not administered opioids and gabapentinoids
		  by that day.

Input dataset:
	Cloned dataset (one record per clone) with the following variables:
		- id_var: the unique ID variable for each eligible patient.
		- assign: (0/1) the assigned treatment strategy for the clone.
		- start_date_d: date of day 0 (i.e., date of SNF admission in DATE9. format).
		- time_to_gaba: the number of days (0, 1, etc.) between day 0 and the first
                        gabapentinoid administration.
			- Can be missing, 0 (if administration is on same day as day 0), or a positive integer.
			- If missing, patient was never administered a gabapentinoid on or after day 0 in an LTCDC-participating SNF.
		- first_gaba: date of first gabapentinoid administration following SNF admission; in DATE9. format.
			- Can be missing if no gabapentinoid administration during SNF stay.
		- time_to_opioid: the number of days (0, 1, etc.) between day 0 and the first opioid administration.
			- Can be missing, 0 (if administration occurred on the same day as day 0), or a positive integer.
			- If missing, the patient was never administered an opioid on or after day 0 in an LTCDC-participating SNF.

Output:
	Cloned dataset (one record per clone) with an artificial censoring date variable (art_cens_d; in DATE9. format)
		- Note that for a given id_var, art_cens_d may differ for each of its clones, since
		  the date of artificial censoring is based on when the clone deviated from that particular treatment
          strategy the clone was assigned.
			- If art_cens_d is missing for a given clone, they were not artificially censored during follow-up.
		- Clones with assign = 1 and art_cens_d = . should have a value for trt_init_d, which is the date by which both
		  medications were initiated during the grace period.

**************************************************************************/

/*** Insert libraries or %include statements here ***/

/*** Options ***/
options mprint mlogic symbolgen;

/*** Macro to create artificial censoring date variable ***/
%macro create_artificial_censoring_vars(
	in_ds=,
	grace_day_end=,
	out_ds=,
	rtf_out=
);

	/* Optionally write to RTF file */
	%if %length(&rtf_out) > 0 %then %do;
		ods rtf file = "&rtf_out.";
	%end;

	/********************************************************************/
	/*** 1. Create art_cens_d variable according to each assign value ***/
	/********************************************************************/

	/* Note: art_cens_d is the date of artificial censoring due to deviation from the assigned treatment strategy */

	data stay2;
		set &in_ds.;

		/* opioid-only treatment strategy */
		if assign = 0 then do;
			/* if patient has gaba admin during grace period, set artificial censoring date to day of first gaba admin */
			if time_to_gaba ne . and time_to_gaba <= &grace_day_end. then do;
				art_flag_1 = 1;
				art_cens_d_1 = first_gaba;
			end;
			/* if patient does not have opioid admin during grace period, set artificial censoring date to last day of grace period */
			if time_to_opioid = . or time_to_opioid > &grace_day_end. then do;
				art_flag_2 = 1;
				art_cens_d_2 = start_date_d + &grace_day_end.;
			end;
			
			/* set artificial censoring date to the minimum date established */
			if art_flag_1 = 1 or art_flag_2 = 1 then do;
				art_cens_d = min(art_cens_d_1, art_cens_d_2); /* takes the lowest non-missing value */
			end;

			/* if no artificial censoring, set to missing */
			if art_flag_1 ne 1 and art_flag_2 ne 1 then do;
				art_flag_3 = 1;
				art_cens_d = .;
			end;

		end;

		/* gabapentinoids + opioids treatment strategy */
		if assign = 1 then do;
			/* if patient does not have an opioid admin during the grace period, set artificial censoring date to last day of grace period */
			if time_to_opioid = . or time_to_opioid > &grace_day_end. then do;
				art_flag_1 = 1;
				art_cens_d_1 = start_date_d + &grace_day_end.;
			end;
			/* if patient does not have gaba admin during the grace period, set artificial censoring date to last day of grace period */
			if time_to_gaba = . or time_to_gaba > &grace_day_end. then do;
				art_flag_2 = 1;
				art_cens_d_2 = start_date_d + &grace_day_end.;
			end;

			/* set artificial censoring date to the minimum date established */
			if art_flag_1 = 1 or art_flag_2 = 1 then do;
				art_cens_d = min(art_cens_d_1, art_cens_d_2);
			end;
			
			/* if no artificial censoring, set to missing AND create date for treatment strategy initiation */
			if art_flag_1 ne 1 and art_flag_2 ne 1 then do;
				art_flag_3 = 1;
				art_cens_d = .;

				/* Create date for treatment strategy initiation (i.e., date second drug was taken) */
				opioid_d = time_to_opioid + start_date_d;
				gaba_d = time_to_gaba + start_date_d;
				trt_init_d = max(opioid_d, gaba_d);

			end;


		end;
		
		format art_cens_d art_cens_d_1 art_cens_d_2 opioid_d gaba_d trt_init_d date9.;

		label
			art_cens_d = "Date of artificial censoring (according to assign value)"
			trt_init_d = "Date of treatment strategy initiation for assign=1 where art_cens_d=."
		;

	run;

	/***********************/
	/*** 2. Housekeeping ***/
	/***********************/

	/* Save dataset and drop temporary variables */
	data &out_ds.;
		set stay2 (drop = art_cens_d_1 art_cens_d_2 opioid_d gaba_d art_flag_1 art_flag_2 art_flag_3);
	run;

	/* Close RTF file */
	%if %length(&rtf_out) > 0 %then %do;
	  ods rtf close;
	%end;

	%put NOTE: create_artificial_censoring_vars completed on %sysfunc(today(), worddate.).;

%mend;

/*** Example macro call ***/
%create_artificial_censoring_vars(
	in_ds=lib.gaba_opioids_clones_ex,                                 /* Input dataset */
	grace_day_end=2,                                                  /* Last day of the grace period (day 0, 1, *2*) */
	out_ds=lib.gaba_ops_artcens_ex,                                   /* Output dataset */
	rtf_out=/your/path/create_artificial_censoring_vars_&RUNDATE..rtf /* RTF output */
);

/* End of create_artificial_censoring_vars macro */
