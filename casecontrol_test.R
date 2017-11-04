library(Ri2b2casecontrol)

system.time({

  dem <- read.csv("c:/projects/test_data/test_casecontrol_patients.csv", na.strings = "NULL", fileEncoding = "UTF-8-BOM")
  df <- read.csv("c:/projects/test_data/test_casecontrol_analysistable.csv", na.strings = "NULL", fileEncoding = "UTF-8-BOM")
  dict_tbl <- read.csv("c:/projects/test_data/test_casecontrol_dictionary.csv", na.strings= "NULL", fileEncoding = "UTF-8-BOM")

  c <- case_control(p=dem,
                    d=df,
                    dict=dict_tbl,
                    exposure_cds = c('RXNORM:8640', 'RXNORM:5492','NDFRT:N0000029116','NDFRT:N0000029132','NDFRT:N0000029168',
                                     'NDFRT:N0000029178','RXNORM:7646','RXNORM:40790','RXNORM:283742'),
                    riskWindow_daysPreIndex_start=730,
                    riskWindow_daysPreIndex_end=30)

})


write.csv(c$results, "c:/projects/test_data/casecontrol_results.csv")

View(c$results %>% arrange(desc(clogit_OR)))

#plot(jitter(a$age_at_index,0.3), a$index_year, col=a$cohort)

#tableone::CreateTableOne(data=c$cohort_file, strata = "cohort")
#d %>% group_by(concept_cd) %>% summarise(n())
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


