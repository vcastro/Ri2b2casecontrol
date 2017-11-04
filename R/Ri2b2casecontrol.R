#' \code{Ri2b2casecontrol} package
#'
#' i2b2 R package to run case-control analysis
#'
#' See the README on
#' \href{https://gitlab.partners.org/vc070/Ri2b2casecontrol}{GitHub}
#'
#' @docType package
#' @name Ri2b2casecontrol
#' @importFrom dplyr %>%
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(
  c(".",
  "RACE_CD", "SEX_CD", "age_at_index", "birth_date", "cohort", "concept_cd", "concept_date",
  "concept_order", "drug_date", "index_date", "index_year", "k2k_control",
  "match_strata", "matched","min_drug_order","patient_num","riskWindow_end",
  "riskWindow_start","visit_date","visit_year","visit_year.x", "visit_year.y",
  "visits_in_window")
)
