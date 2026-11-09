/**************************************************************************
Purpose: Create outcome-specific person-period datasets for pooled logistic regression models.

Input dataset:

	Cloned, outcome-specific wide-format datasets with the following variables:
		- &id_var.: the unique ID variable for each eligible patient.
		- assign: (0/1) indicates treatment assignment for the clone.
		- start_date_d: date of day 0 (i.e., date of SNF admission in DATE9. format).
		- art_cens_d: date of artificial censoring in DATE9. format. 
					  If missing, clone was not artificially censored.
			- 00_create_artificial_censoring_vars.sas shows how art_cens_d was generated.
		- trt_init_d: for the gabapentinoids + opioids strategy (assign = 1), trt_init_d
					  is the date that both gabapentinoids and opioids were administered by
					  during the grace period; in DATE9. format.
		- outcome variables: variables with the date a given outcome occurred in DATE9. format.
			- the date can be anytime after day 0, or it can be set to missing.
			- these outcome variable names are as follows:
				- outcome_death: all-cause mortality
				- outcome_hosp: all-cause rehospitalization
				- outcome_fri: fall-related injury rehospitalization
				- outcome_snf_nohosp: SNF discharge without hospital transfer
		- min_nat_cens: the first date of natural censoring during follow-up; in DATE9. format.
		- time_final_d: last date of follow-up for the clone; in DATE9. format.

This program uses the following hierarchy for events:
	- outcome > artificial censoring > natural censoring
		- This means that, in instances where multiple events occur on the same date, we use the above hierarchy
		  to decide which event occurred (i.e., resolve ties).
			- For example, if the clone experienced the outcome and was artificially censored on the same date, we consider
			  the outcome to have taken place (and not artificial censoring).
		- If more than one event occurred within a discrete-time interval, but one clearly occurred on an earlier day in said interval,
		  that event takes precedence.	

Output dataset:

	Cloned, outcome-specific person-period datasets.
		- Each row represents an interval of time the clone contributed follow-up. 
		- The following variables are included:
			- &id_var: the unique ID variable for each eligible patient.
			- assign: (0/1) indicates treatment assignment for the clone.
			- artificial_cens: (0/1) Indicator of whether artificial censoring occurred in that interval.
			- trt_init: (0/1) For assign = 1, whether the strategy was initiated in that interval.
			- trt_init_prior: (0/1) For assign = 1, flags strategy initiation before the final grace-period day, from that initiation interval onward;
			                   intended for restriction at the final grace-period interval.
			- natural_cens: (0/1) Indicator of whether natural censoring occurred in that interval.
			- event: (0/1) Indicator of whether the outcome occurred in that interval.
			- t_intrv: The # for the discrete-time interval (e.g., 0, 1, 2, etc.).
			- intrv_start_day: The first day in the given discrete-time interval (e.g., 0).
			- intrv_end_day: The last day in the given discrete-time interval (e.g., 2).

**************************************************************************/

/*** Insert libraries or %include statements here ***/

/*** Options ***/
options mprint mlogic symbolgen;

/*** Macro to create person-period datasets ***/
%macro prep_person_period_ds(
	in_ds_prefix=,
	out_ds_prefix=,
	id_var=,
	outcomes=,
	grace_day_end=,
	pre_grace_length=,
	post_grace_length=,
	study_end_day=,
	xlsx_out=
);

	ods excel file="&xlsx_out."
			  options(sheet_interval = "none"
                      embedded_titles = "yes"
			  );

	/**************************************************************/
	/*** 1. Create separate row for each discrete-time interval ***/
	/**************************************************************/

	%macro by_outcome(outcome);
		
		data pp_&outcome.;
			set &in_ds_prefix._&outcome.;
			format intrv_start_d intrv_end_d date9.;
			
			/* Generate discrete-time intervals */
			/* Initialize variables */
			t_intrv = 0;
			intrv_start_day = 0;

			do while (intrv_start_day <= &study_end_day.);

				/* Pre-grace period discrete-time intervals */
				if intrv_start_day <= &grace_day_end. then do;
					intrv_end_day = intrv_start_day + &pre_grace_length. - 1;

					/* Ensure pre-grace interval doesn't spill into post-grace period */
					if intrv_end_day > &grace_day_end. then intrv_end_day = &grace_day_end.;

				end;

				/* Post-grace period discrete-time intervals */
				else do;
					intrv_end_day = intrv_start_day + &post_grace_length. - 1;
				end;

				/* Cap the final interval at the end of follow-up */
				if intrv_end_day > &study_end_day. then intrv_end_day = &study_end_day.;

				/* Calendar date calculations */
				intrv_start_d = start_date_d + intrv_start_day;
				intrv_end_d   = start_date_d + intrv_end_day;

				output;

				/* Increment for the next interval */
				intrv_start_day = intrv_end_day + 1;
				t_intrv = t_intrv + 1;

			end;

			label
				intrv_start_day="Index day for start of interval"
				intrv_end_day="Index day for end of interval"
				intrv_start_d="Date for start of interval"
				intrv_end_d="Date for end of interval"
			;

		run;

		/*** Check that the person-period dataset was generated correctly ***/
		proc sql;
			create table pp_check as
			select distinct t_intrv,
			       intrv_start_day,
				   intrv_end_day
			from pp_&outcome.;
		quit;

		ods excel options(sheet_interval = "now" sheet_name = "&Outcome. Intervals");

		title1 "Discrete-time intervals in the person-period dataset";
		title2 "Dataset: pp_check; Outcome: &outcome.";
		proc print data = pp_check;
		run;
		title;

		/**********************************************************************************************************************************************/
		/*** 2. Remove rows where the start of the discrete-time interval is after time_final_d (date of outcome or censoring event for that clone) ***/
		/**********************************************************************************************************************************************/

		data pp_&outcome.2;
			set pp_&outcome.;

			if time_final_d < intrv_start_d then delete;

		run;

		ods excel options(sheet_interval = "now" sheet_name = "&Outcome. Select Obs");

		title1 "50 observations (after removing rows after outcome/censoring event)";
		title2 "Dataset: pp_&outcome.2; Outcome: &outcome.";
		proc print data = pp_&outcome.2 (obs = 50);
			var &id_var. assign t_intrv intrv_start_d intrv_end_d time_final_d;
		run;

		/***********************************************************************************************************************************/
		/*** 3. Create indicator (0/1) variables for artificial censoring, natural censoring, and the outcome (before dealing with ties) ***/
		/***********************************************************************************************************************************/
		
		proc sort data = pp_&outcome.2;
			by &id_var. assign;
		run;

		data pp_&outcome.3;
			set pp_&outcome.2;
			by &id_var. assign;

			retain trt_init_prior 0;

			/* Reset trt_init_prior to 0 at the beginning of each clone arm */
			if first.assign then trt_init_prior = 0;

			/* artificial censoring indicator variable */
			if art_cens_d ne . and intrv_start_d <= art_cens_d <= intrv_end_d then artificial_cens = 1;
			else artificial_cens = 0;

			/* natural censoring indicator variable */
			if min_nat_cens ne . and intrv_start_d <= min_nat_cens <= intrv_end_d then natural_cens = 1;
			else natural_cens = 0;

			/* outcome indicator variable */
			if &outcome. ne . and intrv_start_d <= &outcome. <= intrv_end_d then outcome = 1;
			else outcome = 0;

			/* treatment initiation variable (only relevant for assign = 1) */
				/* treatment initiation in current interval */
			if trt_init_d ne . and intrv_start_d <= trt_init_d <= intrv_end_d then trt_init = 1;
			else trt_init = 0;

			/* make trt_init_prior variable at end of grace period for assign = 1 */
			/* Record if initiation occurred BEFORE the final day of the grace period */
			if assign = 1 and trt_init = 1 and intrv_end_day < &grace_day_end. then trt_init_prior = 1;
				
			label
				artificial_cens = "(0/1) Whether artificial censoring happened during interval"
				natural_cens = "(0/1) Whether natural censoring happened during interval"
				outcome = "(0/1) Whether outcome happened during interval"
				trt_init =  "(0/1) For assign = 1; whether treatment was initiated in the given interval"
				trt_init_prior = "(0/1) For assign = 1; whether treatment was initiated before the last day of the grace period"
			;

		run;

		/*** Check how often more than one indicator is 1 in a given interval (i.e., there's a tie within an interval) ***/
		proc sort data = pp_&outcome.3;
			by assign &id_var. t_intrv;
		run;

		data pp_&outcome.4;
			set pp_&outcome.3;

			if artificial_cens = 1 and natural_cens = 1 then cens_tie = 1;
			else cens_tie = 0;

			if artificial_cens = 1 and outcome = 1 then art_out_tie = 1;
			else art_out_tie = 0;

			if natural_cens = 1 and outcome = 1 then nat_out_tie = 1;
			else nat_out_tie = 0;

			if sum(of artificial_cens natural_cens outcome) >= 2 then any_tie = 1;
			else any_tie = 0;

			label
				cens_tie = "Artificial censoring and natural censoring occur in the same interval"
				art_out_tie = "Artificial censoring and outcome occur in the same interval"
				nat_out_tie = "Natural censoring and outcome occur in the same interval"
				any_tie = "Any tie in censoring or outcome events in the same interval"
			;

		run;

		/*** Group by clone to check overall frequency of these ties ***/
		proc sql;
			create table check_ties as
			select &id_var.,
				   assign,
				   max(cens_tie) as cens_tie,
				   max(art_out_tie) as art_out_tie,
				   max(nat_out_tie) as nat_out_tie,
				   max(any_tie) as any_tie
			from pp_&outcome.4
			group by &id_var., assign;
		quit;

		ods excel options(sheet_interval = "now" sheet_name = "&Outcome. Ties");

		title1 "Frequency of censoring/outcome ties per clone";
		title2 "Dataset: check_ties; Outcome: &outcome.";
		proc freq data = check_ties;
			tables cens_tie art_out_tie nat_out_tie any_tie / list missing;
		run;

		/****************************************************************************************************************************************/
		/*** 4. In instances of ties, do the following:                                                                                       ***/
		/***    1. If one censoring/outcome event occurred earlier in the interval, keep that = 1.                                            ***/
		/***    2. For censoring/outcome events that occur on the same day, give precedence in order of outcome > artificial > natural        ***/
		/****************************************************************************************************************************************/

		data pp_&outcome.5;
			set pp_&outcome.4;

			if any_tie = 1 then do;

				/* Find the earliest observed date */
				earliest_date = min(
					art_cens_d,
					min_nat_cens,
					&outcome.
				);
				
				/* Count how many events occur on that earliest date */
				n_earliest = 0;
				if art_cens_d          = earliest_date then n_earliest + 1;
				if min_nat_cens        = earliest_date then n_earliest + 1;
				if &outcome.           = earliest_date then n_earliest + 1;

				/* Initialize indicators */
				resolved_art_cens = 0;
				resolved_nat_cens = 0;
				resolved_outcome  = 0;
				day_level_tie     = 0;

				/* Case 1: true day-level tie */
				if n_earliest > 1 then do;
					day_level_tie = 1;

					/* Precedence: outcome > artificial > natural */
					if &outcome.               = earliest_date then resolved_outcome = 1;
					else if art_cens_d         = earliest_date then resolved_art_cens = 1;
					else if min_nat_cens       = earliest_date then resolved_nat_cens = 1;
				end;

				/* Case 2: one event clearly occurs first */
				else do;
					day_level_tie = 0;

					if art_cens_d              = earliest_date then resolved_art_cens = 1;
					else if min_nat_cens       = earliest_date then resolved_nat_cens = 1;
					else if &outcome.          = earliest_date then resolved_outcome = 1;
				end;
				
			end;

			/* Clean up records without a tie */
			else do;
				day_level_tie = 0;

				resolved_art_cens = artificial_cens;
				resolved_nat_cens = natural_cens;
				resolved_outcome = outcome;
			end;

			format earliest_date date9.;

			label
				day_level_tie = "Censoring or outcome tie on the same exact day"
				resolved_art_cens = "(0/1) Whether artificial censoring happened during interval (after resolving ties)"
				resolved_nat_cens = "(0/1) Whether natural censoring happened during interval (after resolving ties)"
				resolved_outcome = "(0/1) Whether the outcome happened during interval (after resolving ties)"
			;

		run;

		/*******************************************************************************/
		/*** 5. Save outcome-specific person-period datasets for subsequent modeling ***/
		/*******************************************************************************/

		data &out_ds_prefix._&outcome. (keep = &id_var. assign artificial_cens natural_cens trt_init trt_init_prior event t_intrv intrv_start_day intrv_end_day);
			set pp_&outcome.5 (drop = artificial_cens natural_cens outcome event);

			artificial_cens = resolved_art_cens;
			natural_cens = resolved_nat_cens;
			event = resolved_outcome;

		run;

	%mend;

	/* Run for each outcome */
	%macro run_pp(outcomes_list=);

		%local i current_outcome;

		/* Loop through each outcome in the list */
		%do i = 1 %to %sysfunc(countw(&outcomes.));

			/* Extract the current outcome name */
			%let current_outcome = %scan(&outcomes_list, &i);

			/* Execute the macro for this outcome */
			%by_outcome(&current_outcome.);

		%end;

	%mend;

	%run_pp(outcomes_list=&outcomes.);

	ods excel close;

	/***********************/
	/*** 6. Housekeeping ***/
	/***********************/

	%put NOTE: prep_person_period_ds completed on %sysfunc(today(), worddate.).;

%mend;

/* List of outcome variables */
%let outcomes_list = 
	outcome_death       /* All-cause mortality */
	outcome_hosp        /* All-cause rehospitalization */
	outcome_fri         /* fall-related injury rehospitalization */
	outcome_snf_nohosp  /* SNF discharge without hospital transfer */
;

/*** Example macro call ***/
%prep_person_period_ds(
	in_ds_prefix=lib.wide_follow,     /* Prefix of outcome-specific wide-format input datasets (from 01_prep_follow_up_wide_ds.sas) */
	out_ds_prefix=lib.pp,             /* Prefix of outcome-specific person-period output datasets (e.g., lib.pp_outcome_death) */
	id_var=id_var,                    /* Name of ID variable that uniquely identifies each patient */
	outcomes=&outcomes_list.,         /* List of outcome variable names */
	grace_day_end=2,                  /* Last day of the grace period (day 0, 1, *2*) */
	pre_grace_length=1,               /* Length of discrete-time intervals during the grace period (e.g., 1-day intervals) */
	post_grace_length=3,              /* Length of discrete-time intervals after the end of the grace period (e.g., 3-day intervals) */
	study_end_day=98,                 /* Last day of follow-up (day 98) */
	xlsx_out=/your/path/prep_pp_dataset_grace2_&RUNDATE..xlsx /* Excel output */
);

/* End of prep_person_period_ds macro */
