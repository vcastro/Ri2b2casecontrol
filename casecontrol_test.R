library(dplyr)
library(lubridate)
library(Ri2b2matchcontrols)


p <- read.csv("c:/projects/test_data/test_casecontrol_patients.csv", na.strings = "NULL")
d <- read.csv("c:/projects/test_data/test_casecontrol_analysistable.csv", na.strings = "NULL")

names(d)[1] <- "patient_num"
names(p)[1] <- "patient_num"


outcome <- "ICD10:M81"
visit_cd <- "L_!1802"
exposure_cds <- c("RXNORM:8640", "RXNORM:5492","NDFRT:N0000029116","NDFRT:N0000029132","NDFRT:N0000029168","NDFRT:N0000029178")


observationWindow_daysPreIndex_start <- 730
observationWindow_daysPreIndex_end <- 30


#patient cleanup%>% 
#filter(index_date < as.Date(death_date) | is.na(death_date))



visits <- d %>% filter(concept_cd == visit_cd) %>% 
  mutate(visit_date = as.Date(concept_date),
         visit_year = year(visit_date)) %>% 
  select(patient_num, visit_date, concept_order, visit_year)


cases <- d %>% filter(concept_cd == outcome & concept_order == 1) %>% 
  mutate(index_date = as.Date(substring(concept_date, 0, 10)),
         observationWindow_start = index_date - observationWindow_daysPreIndex_start,
         observationWindow_end = index_date - observationWindow_daysPreIndex_end) %>% 
  select(patient_num, concept_cd, index_date, observationWindow_start, observationWindow_end) %>% 
  left_join(p, by="patient_num") %>% 
  mutate(age_at_index = as.integer((index_date - as.Date(birth_date))/365.25)) 
#remove age_at_index outliers in cases


eligible_cases <- cases %>% left_join(visits, by="patient_num") %>% 
  filter(visit_date >= observationWindow_start & visit_date <= observationWindow_end) %>% 
  select(patient_num) %>% 
  distinct(patient_num)


cases_final <- cases %>% inner_join(eligible_cases, by="patient_num") %>% 
  mutate(cohort = as.factor("Case"),
         index_year = year(index_date)) %>% 
  select(patient_num, cohort, index_date, index_year, age_at_index, observationWindow_start, observationWindow_end, SEX_CD, RACE_CD)


eligible_controls <- p %>% filter(!patient_num %in% cases$patient_num)


controlpool <- visits %>% inner_join(eligible_controls, by="patient_num") %>% 
  mutate(age_at_index = as.integer((visit_date - as.Date(birth_date))/365.25)) %>% 
  group_by(patient_num, visit_year) %>% 
  sample_n(1) %>% 
  mutate(index_date = visit_date,
         observationWindow_start = index_date - observationWindow_daysPreIndex_start,
         observationWindow_end = index_date - observationWindow_daysPreIndex_end) %>% 
  select(patient_num, index_date, visit_year, age_at_index, observationWindow_start, observationWindow_end) %>% 
  ungroup() %>% 
  left_join(visits, by="patient_num") %>% 
  rename(index_year = visit_year.x) %>% select(-visit_year.y) %>% 
  filter(visit_date >= observationWindow_start & visit_date <= observationWindow_end) %>% 
  group_by(patient_num, index_date, index_year, age_at_index, observationWindow_start, observationWindow_end) %>% 
  summarise(visits_in_window = n())  


cohort_file <- controlpool %>% ungroup() %>% 
      mutate(cohort = as.factor("Control")) %>% 
      left_join(p, by="patient_num") %>% 
      select(patient_num, cohort, index_date, index_year, age_at_index, observationWindow_start, observationWindow_end, SEX_CD, RACE_CD) %>% 
      rbind(cases_final)


#filter cases by having 1+ visit in observation_window

insufficient_observation_cases <- cases %>% anti_join(eligible_cases, by="patient_num") %>% select(patient_num)



temp.cohort_file <- as.data.frame(cohort_file)
final_cohorts <- as_data_frame()
results <- as_data_frame()


for (year in sample(unique(cases_final$index_year)))
  {

    df_for_matching <- temp.cohort_file %>% filter(index_year == year)

      
    #do the matching
    m <- Ri2b2matchcontrols::matchcontrols(df_for_matching, controls_to_match = 3, match_variables = c("age_at_index", "SEX_CD", "RACE_CD"))
    

    print(paste0("Index year ", year, "|Cases=", nrow(df_for_matching[which(df_for_matching$cohort == "Case"),]), ";Controls=",nrow(df_for_matching[which(df_for_matching$cohort == "Control"),])))
    
    matched_df <- df_for_matching %>% inner_join(m$match_data, by=c("patient_num","SEX_CD","RACE_CD","age_at_index","cohort")) %>% 
      filter((cohort == "Case" & matched == TRUE) | (cohort == "Control" & k2k_control == TRUE)) %>% 
      select(index_year, patient_num, cohort, match_strata, index_date)
    
    temp.cohort_file <- temp.cohort_file %>% filter(!patient_num %in% matched_df$patient_num)
    
    final_cohorts <- rbind(final_cohorts, matched_df)

    #results <- rbind(results, data_frame(index_year = year, match_results = m))
        
  }
    

a <- final_cohorts %>% 
    mutate(fullstrata = as.integer(paste0(as.character(index_year), as.character(match_strata)))) %>% 
    left_join(cohort_file, by=c("patient_num", "index_year", "index_date", "cohort")) %>% 
    group_by(fullstrata)

analysis_results <- as_data_frame()


for (exp in exposure_cds) {
  
  df_exp <- d %>% filter(concept_cd == exp) %>% mutate(drug_date = as.Date(concept_date))

  t <- a %>%  
    left_join(df_exp, by="patient_num") %>% 
    filter(drug_date >= observationWindow_start & drug_date <= observationWindow_end) %>% 
    group_by(patient_num) %>% 
    summarise(min_drug_order = min(concept_order)) %>% 
    right_join(a, by=c("patient_num")) %>% 
    mutate(exposed = factor(ifelse(is.na(min_drug_order), 0, 1), levels=c(0,1), labels=c("Unexposed", "Exposed")))

  r <- riskratio.boot(table(t$exposed, t$cohort))
  

    result <- data_frame(drug_cd=exp,
                         a=r$data["Unexposed","Control"],
                         b=r$data["Unexposed","Case"],
                         c=r$data["Exposed","Control"],
                         d=r$data["Exposed","Case"],
                         RR=r$measure["Exposed", "estimate"],
                         LL=r$measure["Exposed", "lower"],
                         UL=r$measure["Exposed", "upper"],
                         fisherP=r$p.value["Exposed","fisher.exact"])
           
    analysis_results <- rbind(result, analysis_results) 

}


analysis_results


plot(jitter(a$age_at_index,0.3), a$index_year, col=a$cohort)

#tableone::CreateTableOne(data=a, strata = "cohort")
#table(substring(cases_final$index_date,0,4))
#length(unique(final_cohorts$patient_num))
#table(final_cohorts$cohort, final_cohorts$index_year)
#table(cases_final$cohort, substring(cases_final$index_date,0,4))



#histogram of index_year
#hist(as.numeric(substring(cases_final$index_date, 0,4)))
#hist(cases$age_at_index)
#barplot(cases$SEX_CD)
#barplot(prop.table(table(a$age_at_index, a$cohort),2), beside=TRUE)
#barplot(prop.table(table(a$SEX_CD, a$cohort),2), beside=TRUE)
#hist(a$cohort, a$age_at_index)

