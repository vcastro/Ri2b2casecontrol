#' case_control
#'
#' The main case control function to run a case control analysis from i2b2-generated data files
#'
#' @param dem_tbl A data frame of patients demographics
#' @param data_tbl A data frame of patient data in a tall file
#' @param dict_tbl A data frame of descriptions for each concept_cd in the data_tbl
#' @param exposure_cds A string vector of exposures codes to test
#' @param visit_cd The code used for visits (default is V)
#' @param outcome_cd The code used for the outcome in the data_tbl (default is O)
#' @param riskWindow_daysPreIndex_start The number of days prior to the index date to start the risk window (default is 365)
#' @param riskWindow_daysPreIndex_end The number of days prior to the index date to end the risk window (default is 1)
#'
#' @return A list with the elements
#' \item{result}{A data frame with results including a row for each defined exposure}
#' \item{match_results}{A list returned from the matching function }
#' \item{cohort_file}{The full patient data frame with case control status after matching
#' @export
#'
#' @examples
#' #Not run
case_control <- function (dem_tbl,
                          data_tbl,
                          dict_tbl,
                          exposure_cds,
                          visit_cd = "V",
                          outcome_cd = "O",
                          riskWindow_daysPreIndex_start = 365,
                          riskWindow_daysPreIndex_end = 1
                         )
{

outcome_name <- as.character(dict[which(dict$concept_cd == outcome_cd), "name_char"])

#patient cleanup%>%
#filter(index_date < as.Date(death_date) | is.na(death_date))

visits <- d %>% filter(concept_cd == visit_cd) %>%
  mutate(visit_date = as.Date(concept_date),
         visit_year = year(visit_date)) %>%
  select(patient_num, visit_date, concept_order, visit_year)


cases <- d %>% filter(concept_cd == outcome_cd & concept_order == 1) %>%
  mutate(index_date = as.Date(concept_date),
         riskWindow_start = index_date - riskWindow_daysPreIndex_start,
         riskWindow_end = index_date - riskWindow_daysPreIndex_end) %>%
  select(patient_num, concept_cd, index_date, riskWindow_start, riskWindow_end) %>%
  left_join(p, by="patient_num") %>%
  mutate(age_at_index = as.integer((index_date - as.Date(birth_date))/365.25))
#remove age_at_index outliers in cases


cases_final <- cases %>% left_join(visits, by="patient_num") %>%
  filter(visit_date >= riskWindow_start & visit_date <= riskWindow_end) %>%
  mutate(cohort = as.factor("Case"),
         index_year = year(index_date)) %>%
  group_by(patient_num, cohort, index_date, index_year, age_at_index, riskWindow_start, riskWindow_end, SEX_CD, RACE_CD) %>%
  summarise(visits_in_window = n()) %>%
  ungroup()

#filter cases by having 1+ visit in observation_window
insufficient_observation_cases <- cases %>% anti_join(cases, by="patient_num") %>% select(patient_num)


eligible_controls <- p %>% filter(!patient_num %in% cases$patient_num)


controlpool <- visits %>% inner_join(eligible_controls, by="patient_num") %>%
  mutate(age_at_index = as.integer((visit_date - as.Date(birth_date))/365.25)) %>%
  group_by(patient_num) %>%
  sample_n(1) %>%
  mutate(index_date = visit_date,
         riskWindow_start = index_date - riskWindow_daysPreIndex_start,
         riskWindow_end = index_date - riskWindow_daysPreIndex_end) %>%
  select(patient_num, index_date, visit_year, age_at_index, riskWindow_start, riskWindow_end) %>%
  ungroup() %>%
  left_join(visits, by="patient_num") %>%
  rename(index_year = visit_year.x) %>% select(-visit_year.y) %>%
  filter(visit_date >= riskWindow_start & visit_date <= riskWindow_end) %>%
  group_by(patient_num, index_date, index_year, age_at_index, riskWindow_start, riskWindow_end) %>%
  summarise(visits_in_window = n())


cohort_file <- controlpool %>% ungroup() %>%
  mutate(cohort = as.factor("Control")) %>%
  left_join(p, by="patient_num") %>%
  select(patient_num, cohort, index_date, index_year, age_at_index, riskWindow_start, riskWindow_end, SEX_CD, RACE_CD, visits_in_window) %>%
  rbind(cases_final) %>%
  as.data.frame()


#do the matching
m <- Ri2b2matchcontrols::matchcontrols(cohort_file, controls_to_match = 3, match_variables = c("age_at_index", "SEX_CD", "RACE_CD", "index_year"))

final_cohorts <- cohort_file %>%
  inner_join(m$match_data, by=c("patient_num","SEX_CD","RACE_CD","age_at_index","cohort", "index_year")) %>%
  filter((cohort == "Case" & matched == TRUE) | (cohort == "Control" & k2k_control == TRUE)) %>%
  select(index_year, patient_num, cohort, match_strata, index_date) %>%
  left_join(cohort_file, by=c("patient_num", "index_year", "index_date", "cohort"))


analysis_results <- as_data_frame()

for (exp in exposure_cds) {

  df_exp <- d %>% filter(concept_cd == exp) %>% mutate(drug_date = as.Date(concept_date))

  exp_name <- as.character(dict[which(dict$concept_cd == exp), "name_char"])

  t <- final_cohorts %>%
    left_join(df_exp, by="patient_num") %>%
    filter(drug_date >= riskWindow_start & drug_date <= riskWindow_end) %>%
    group_by(patient_num) %>%
    summarise(min_drug_order = min(concept_order)) %>%
    right_join(final_cohorts, by=c("patient_num")) %>%
    mutate(exposed = factor(ifelse(is.na(min_drug_order), 0, 1), levels=c(0,1), labels=c("Unexposed", "Exposed")))

  r <- epitools::riskratio.boot(table(t$exposed, t$cohort))

  exp.logit <- glm(exposed ~ cohort+index_year+log(visits_in_window), data=t, family="binomial")
  or <- exp(cbind(OR = coef(exp.logit), confint.default(exp.logit)))
  s<-summary(exp.logit)

  exp.clogit <- survival::clogit(as.integer(t$exposed) ~ cohort +index_year+log(visits_in_window)+strata(match_strata), data=t)
  cor <- exp(cbind(OR = coef(exp.clogit), confint.default(exp.clogit)))
  cs<-summary(exp.clogit)


  result <- data_frame(drug_cd=exp,
                       drug_name=exp_name,
                       outcome_name=outcome_name,
                       a=r$data["Unexposed","Control"],
                       b=r$data["Unexposed","Case"],
                       c=r$data["Exposed","Control"],
                       d=r$data["Exposed","Case"],
                       RR=r$measure["Exposed", "estimate"],
                       LL=r$measure["Exposed", "lower"],
                       UL=r$measure["Exposed", "upper"],
                       fisherP=r$p.value["Exposed","fisher.exact"],
                       logit_OR=or["cohortCase","OR"],
                       logit_LL=or["cohortCase","2.5 %"],
                       logit_UL=or["cohortCase","97.5 %"],
                       logit_P=s$coefficients["cohortCase", "Pr(>|z|)"],
                       clogit_OR=cor["cohortCase","OR"],
                       clogit_LL=cor["cohortCase","2.5 %"],
                       clogit_UL=cor["cohortCase","97.5 %"],
                       clogit_P=cs$coefficients["cohortCase", "Pr(>|z|)"]
  )

  analysis_results <- rbind(result, analysis_results)

}

list(results=analysis_results, match_results=m, cohort_file=final_cohorts)

}
