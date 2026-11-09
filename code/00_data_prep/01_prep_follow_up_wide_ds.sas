/**************************************************************************
Purpose: Prep variables related to follow-up for outcome-specific wide-format datasets.
		 You can think of this dataset as a stepping stone to the person-period dataset
         with discretized intervals of time.
		 The dataset output from this program can also be used to fit a naive, 
		 unadjusted Kaplan-Meier survival estimator (i.e., using the event, time_final, and assign variables).

Input dataset:
	Cloned cohort dataset with the following variables:
		- id_var: the unique ID variable for each eligible patient.
		- assign: (0/1) indicates treatment assignment for the clone.
		- start_date_d: date of day 0 (i.e., date of SNF admission in DATE9. format).
		- art_cens_d: date of artificial censoring in DATE9. format. 
					  If missing, clone was not artificially censored.
			- 00_create_artificial_censoring_vars.sas shows how art_cens_d was generated.
		- outcome variables: variables with the date a given outcome occurred in DATE9. format.
			- the date can be anytime after day 0, or it can be set to missing.
			- these outcome variable names are as follows:
				- outcome_death: all-cause mortality
				- outcome_hosp: all-cause rehospitalization
				- outcome_fri: fall-related injury rehospitalization
				- outcome_snf_nohosp: SNF discharge without hospital transfer
		- natural censoring variables: variables with the date natural censoring events occurred.
		  (e.g., end of follow-up, study period end, medicare disenrollment, etc.)
			- the date can be anytime on or after day 0, or it can be set to missing.
			- these natural censoring variables are as follows:
				- censor_follow_end: Last date of follow-up (i.e., 99 days)
				- censor_study_end: December 31, 2023
				- censor_disenroll: Last date of Medicare enrollment (or first day of enrollment in Medicare Advantage)
				- outcome_death: For non-mortality outcomes, follow-up ends at date of death.
				- outcome_snf_hosp: For the SNF discharge outcome, follow-up ends on the date of hospital transfer.

	This program uses the following hierarchy for events:
		- outcome > artificial censoring > natural censoring
			- This means that, in instances where multiple events occur on the same date, we use
			  the above hierarchy to decide which event occurred (i.e., resolve ties).
				- For example, if the clone experienced the outcome and was artificially censored on the same day, 
				  we consider the outcome to have occurred (and don't count it as artificial censoring).

Output dataset:
	Cloned cohort dataset for each of the outcomes with the following variables added:
		- time_final_d: Last date of follow-up for the clone; in DATE9. format.
		- time_final: Last day of follow-up relative to day 0 (e.g., 2, 5, 30, etc.).
		- event: (0/1) Indicator of whether the outcome occurred during follow-up.
		- art_cens_on_final: (0/1) Indicator of whether artificial censoring occurred during follow-up.
		- nat_cens_on_final: (0/1) Indicator of whether natural censoring occurred during follow-up.

		Note that only one of the event, art_cens_on_final, and nat_cens_on_final variables will be 1 for a given clone.
		Ties have been resolved.

**************************************************************************/

/*** Insert libraries or %include statements here ***/


/*** Options ***/
options mprint mlogic symbolgen;

/*** Macro to prep variables related to follow-up in outcome-specific datasets ***/
%macro prep_follow_up_wide_ds(
	in_ds=,
	out_prefix=,
	rtf_out=
);

	/* Optionally write to RTF file */
	%if %length(&rtf_out) > 0 %then %do;
		ods rtf file = "&rtf_out.";
	%end;

	/*******************************************************************************************************************/
	/*** 2. Output separate datasets for each outcome                                                                ***/
	/***    (all-cause mortality, all-cause rehospitalization, fall-related injury rehospitalization, SNF discharge) ***/
	/*******************************************************************************************************************/

	data d_outcome_death d_outcome_hosp d_outcome_fri d_outcome_snf_nohosp;
		set &in_ds.;
	run;

	/********************************************/
	/*** 3. Identify minimum event by outcome ***/
	/********************************************/
	
	%macro min_event_outcome(outcome);

		/*** All-cause mortality ***/
		%if &outcome. = outcome_death %then %do;

			/* Death */
			data &out_prefix._&outcome.;
				set d_&outcome.;

				/* Determine last day of follow-up */
				time_final_d = min(of art_cens_d outcome_death censor_follow_end censor_disenroll censor_study_end);

				/* Determine minimum natural censoring event (before dealing with ties) */
				min_nat_cens = min(of censor_follow_end censor_disenroll censor_study_end);

				format time_final_d min_nat_cens date9.;

				/* Create variable with last day of follow-up as a day relative to day 0 */
				time_final = time_final_d - start_date_d; /* start_date_d is the date of day 0 */

				/* Event occurs ONLY if:
					- outcome exists
					- outcome is exactly time_final_d
				*/
					/* Outcome takes precedence over artificial censoring or natural censoring */
				if outcome_death ne . and
				   outcome_death = time_final_d
				then event = 1;
				else event = 0;

				/* Flag artificial censoring on final day */
				if event = 0 and time_final_d = art_cens_d then 
					art_cens_on_final = 1;
				else art_cens_on_final = 0;

				/* Flag natural censoring on final day */
				if event = 0 and art_cens_on_final = 0 then
					nat_cens_on_final = 1;
				else nat_cens_on_final = 0;
				
				label
					time_final_d = "Last date of follow-up based on artificial censoring, natural censoring, and outcome"
					min_nat_cens = "Minimum date of natural censoring before dealing with ties"
					time_final = "time_final_d - start_date_d"
					event = "1=event occurred on time_final"
					art_cens_on_final = "(0/1) 1=artificial censoring occurred during grace period (ties resolved)"
					nat_cens_on_final = "(0/1) 1=natural censoring occurred during follow-up (ties resolved)"
				;

			run;

		%end;

		/*** SNF discharge without hospital transfer ***/
		%else %if &outcome. = outcome_snf_nohosp %then %do;

			data &out_prefix._&outcome.;
				set d_&outcome.;

				/* Determine last day of follow-up */
				time_final_d = min(of art_cens_d &outcome. outcome_death outcome_snf_hosp censor_follow_end censor_disenroll censor_study_end);

				/* Determine minimum natural censoring event (before dealing with ties) */
				min_nat_cens = min(of outcome_death outcome_snf_hosp censor_follow_end censor_disenroll censor_study_end);

				format time_final_d min_nat_cens date9.;

				/* Create variable with last day of follow-up as a day relative to day 0 */
				time_final = time_final_d - start_date_d;

				/* Event occurs ONLY if:
					- outcome exists
					- outcome is exactly time_final_d
				*/
					/* Outcome takes precedence over artificial censoring or natural censoring */
				if &outcome. ne . and
				   &outcome. = time_final_d
				then event = 1;
				else event = 0;

				/* Flag artificial censoring on final day */
				if event = 0 and time_final_d = art_cens_d then 
					art_cens_on_final = 1;
				else art_cens_on_final = 0;

				/* Flag natural censoring on final day */
				if event = 0 and art_cens_on_final = 0 then
					nat_cens_on_final = 1;
				else nat_cens_on_final = 0;

				label
					time_final_d = "Last date of follow-up based on artificial censoring, natural censoring, and outcome"
					min_nat_cens = "Minimum date of natural censoring before dealing with ties"
					time_final = "time_final_d - start_date_d"
					event = "1=event occurred on time_final"
					art_cens_on_final = "(0/1) 1=artificial censoring occurred during grace period (ties resolved)"
					nat_cens_on_final = "(0/1) 1=natural censoring occurred during follow-up (ties resolved)"
				;

			run;

		%end;

		/*** All-cause rehospitalization and fall-related injury rehospitalization ***/
		%else %do;

			data &out_prefix._&outcome.;
				set d_&outcome.;

				/* Determine last date of follow-up */
				time_final_d = min(of art_cens_d &outcome. outcome_death censor_follow_end censor_disenroll censor_study_end);

				/* Determine minimum natural censoring event (before dealing with ties) */
				min_nat_cens = min(of outcome_death censor_follow_end censor_disenroll censor_study_end);

				format time_final_d min_nat_cens date9.;

				/* Create variable with last day of follow-up as a day relative to day 0 */
				time_final = time_final_d - start_date_d;

				/* Event occurs ONLY if:
						- outcome exists 
						- outcome is exactly time_final_d
				*/
					/* Outcome takes precedence over artificial censoring or natural censoring */
				if &outcome. ne . and 
				   &outcome. = time_final_d
				then event = 1;
				else event = 0;

				/* Flag artificial censoring on final day */
				if event = 0 and time_final_d = art_cens_d then 
					art_cens_on_final = 1;
				else art_cens_on_final = 0;

				/* Flag natural censoring on final day */
				if event = 0 and art_cens_on_final = 0 then
					nat_cens_on_final = 1;
				else nat_cens_on_final = 0;

				label
					time_final_d = "Last date of follow-up based on artificial censoring, natural censoring, and outcome"
					min_nat_cens = "Minimum date of natural censoring before dealing with ties"
					time_final = "time_final_d - start_date_d"
					event = "1=event occurred on time_final"
					art_cens_on_final = "(0/1) 1=artificial censoring occurred during grace period (ties resolved)"
					nat_cens_on_final = "(0/1) 1=natural censoring occurred during follow-up (ties resolved)"
				;

			run;

		%end;

		title1 "Number of outcome events, natural censoring events, and artificial censoring events";
		title2 "Dataset: &out_prefix._&outcome.";
		proc freq data = &out_prefix._&outcome.;
			tables event art_cens_on_final nat_cens_on_final event*art_cens_on_final*nat_cens_on_final / list missing;
		run;

	%mend;

	/* Run for each outcome */
	%min_event_outcome(outcome_death);      /* All-cause mortality */
	%min_event_outcome(outcome_hosp);       /* All-cause rehospitalization */
	%min_event_outcome(outcome_fri);        /* fall-related injury rehospitalization */
	%min_event_outcome(outcome_snf_nohosp); /* SNF discharge without hospital transfer */

	/***********************/
	/*** 4. Housekeeping ***/
	/***********************/

	/* Close RTF file */
	%if %length(&rtf_out) > 0 %then %do;
	  ods rtf close;
	%end;

	%put NOTE: prep_follow_up_wide_ds completed on %sysfunc(today(), worddate.).;

%mend;

/*** Example macro call ***/
%prep_follow_up_wide_ds(
	in_ds=lib.gaba_ops_artcens_ex,                           /* Input dataset */                         
	out_prefix=lib.wide_follow,                              /* Prefix of outcome-specific output dataset (e.g., lib.wide_follow_outcome_death) */
	rtf_out=/your/path/prep_follow_up_wide_ds_&RUNDATE..rtf  /* RTF output; can be left empty if you don't want to produce a .rtf file */
);

/* End of prep_follow_up_wide_ds macro */
