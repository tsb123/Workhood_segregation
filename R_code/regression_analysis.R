##version for publication
library(tidyverse)
library(seg)
library(sf)
library(ggplot2)
library(rmapshaper)
library(tigris)
library(gridExtra)
library(glmnet)
library(plm)
library(grid)
library(gridExtra)
library(sjPlot)

load(file = "datafiles_new/segregration_msa_subset382_joined_RAC_nospatial.RData")
msa_rac<-msa2
load(file = "datafiles_new/segregration_msa_subset382_joined_wac_nospatial.RData")
msa_wac<-msa2
rm(msa2)


### (A) collate dataset ####
rac_colstokeep<-c("CBSAFP", "year", "multigroupd_race_rac", "d_nonhispanichispanic_rac", "multigroupd_education_rac", "multigroupd_profession_rac" ) 
wac_colstokeep<-c("CBSAFP", "year", "multigroupd_race_wac","d_nonhispanichispanic_wac", "multigroupd_education_wac" ,"multigroupd_profession_wac" )
MSA_popprofiletokeep<-c("total","age55_older", "income1250_less",
                        "white", "black", "asian", "hispanic", "female", 
                        "education_college", "education_lesshighschool",
                       "agriculture","mining","utilities","construction","manufacturing","wholesale","retail","transportation",
                       "information","finance","realestate","professional","management","waste","education","healthcare","arts",
                       "accommodation","otherservices","publicadmin" )
msa_ss<-merge(msa_rac[,rac_colstokeep], msa_wac[,c(wac_colstokeep,MSA_popprofiletokeep)], by=c("CBSAFP", "year"),all=T)

##join OTHER VARIABLES 
#  (AI) transit variables #####
source('R_code/road_network_variables.R')
msa_ss[c("patch_bupl","patch_bua","all_bupl","all_bua","distance","mean_local_gridness")]<-NULL
msa_ss<-merge(msa_ss, road_network_vars, by.x=c("CBSAFP", "year"), by.y=c("msaid", "year"), all.x=T)
sum(is.na(subset(msa_ss, year==2016)$distance)) ## 7 missing MSAs
msa_ss$distance_k<-msa_ss$distance/1000

# (AII) ACS data #####
load(file = "datafiles_new/suburban_CBSA_acs.RData")
county_2020_suburban$percent_suburban<-round(county_2020_suburban$suburbantotal/county_2020_suburban$totalsum*100,2)
msa_ss<-merge(msa_ss, county_2020_suburban[c("percent_suburban", "CBSA_Code", "year")], by.x=c("CBSAFP", "year"), by.y=c("CBSA_Code", "year"), all.x=T)

msa_acs<-read.csv("datafiles_new/residential_seg_acs5_2011_2018_msa_race_ethnicityseparate.csv")
msa_edu_acs<-read.csv("datafiles_new/residential_seg_acs5_2011_2018_msa_edu.cvs")

msa_edu_acs<-st_drop_geometry(msa_edu_acs)

msa_acs<-merge(msa_acs[c("CBSAFP","year", "multigroupd_raceseparate_rac", "d_nothispanicethnicityhispanicethnicity_rac")], msa_edu_acs[c("CBSAFP","year","multigroupd_edu_rac")], by=c("CBSAFP", "year"))
rm(msa_edu_acs)

load(file="datafiles_new/census_covariates_msa.RData")
grouped_dfs<-list(all_2011_msa, all_2012_msa, all_2012_msa, all_2013_msa, all_2014_msa, all_2015_msa, all_2016_msa, all_2017_msa, all_2018_msa)
wut<-lapply(seq_along(grouped_dfs), function(i){
  df <- grouped_dfs[[i]] %>% ungroup()
  df$year <- 2010+i
  return(df)
})
ACSvars<-do.call(rbind, wut)
ACSvars<-ACSvars%>%
  rename( white = B03002_003E, total_pop = B02001_001E, total_laborforce=B23025_003E, unemployed_laborforce=B23025_005E, gini=B19083_001E)%>%
  mutate(percent_white=round(white/total_pop*100,2),
         percent_unemployed=round(unemployed_laborforce/total_laborforce*100,2),
         unemployed_laborforce_white=C23002H_008E+C23002H_013E+C23002H_021E+C23002H_026E,
         unemployed_laborforce_nonwhite=unemployed_laborforce- unemployed_laborforce_white,
         total_laborforce_white=unemployed_laborforce_white+C23002H_007E+C23002H_012E+C23002H_020E+C23002H_025E,
         total_laborforce_notwhite=total_laborforce-total_laborforce_white,
         percent_whitelabor_unemployed=round(unemployed_laborforce_white/total_laborforce_white*100,2),
         percent_nonwhitelabor_unemployed=round(unemployed_laborforce_nonwhite/total_laborforce_notwhite*100,2),
         delta_nonwhite_white_unemployed=percent_nonwhitelabor_unemployed-percent_whitelabor_unemployed,
         average_time_work=B08013_001E/B08134_001E,
         percent_publictransport=round(B08134_061E/B08134_001E*100,2),
         pct_carcommutes_over45m=round((B08134_019E+B08134_020E)/B08134_011E*100,2),
         pct_ptcommutes_over45m=ifelse( B08134_061E==0,0,round((B08134_069E+B08134_070E)/B08134_061E*100,2)),
         delta_pct_ptcar_longcommute=pct_ptcommutes_over45m-pct_carcommutes_over45m)
rm(wut)
msa_geom <- st_read('data_files/tl_2020_us_cbsa/tl_2020_us_cbsa.shp')
# Note the relevant shapefile can be downloaded here: https://catalog.data.gov/dataset/tiger-line-shapefile-2020-nation-u-s-core-based-statistical-areas-cbsa/resource/5f1a32d3-e516-41bb-b6fe-3ca19c5958b9
## original shapefile is not included in github repository
msa_geom$areasqkm<-msa_geom$ALAND/1000000  
msa_geom<-st_drop_geometry(msa_geom[c("GEOID", "NAME","areasqkm")])
ACSvars<-merge(ACSvars, msa_geom, by.x="CBSA_Code", by.y="GEOID", all.x=T)                             
rm(grouped_dfs, all_2011_msa, all_2012_msa, all_2013_msa, all_2014_msa, all_2015_msa, all_2016_msa, all_2017_msa, all_2018_msa)
ACSvars$popdensity_sqkm=as.numeric(ACSvars$total_pop/ACSvars$areasqkm)

msa_acs<-merge(msa_acs, ACSvars, by.x=c("CBSAFP", "year"), by.y=c("CBSA_Code", "year"), all.x=T)
unique(msa_rac$CBSAFP[!msa_rac$CBSAFP%in%ACSvars$CBSA_Code]) ## perfect join
msa_acs<-msa_acs%>%rename(multigroupd_race_rac=multigroupd_raceseparate_rac,
                          d_nonhispanichispanic_rac=d_nothispanicethnicityhispanicethnicity_rac,
                          multigroupd_education_rac=multigroupd_edu_rac)
#JOIN BACK TO MAIN 
msa_acs_ss<-st_drop_geometry(msa_acs[c("CBSAFP","year", "multigroupd_race_rac", "d_nonhispanichispanic_rac","multigroupd_education_rac" ,"total_pop", "percent_white", "delta_nonwhite_white_unemployed", 
                                       'gini', 'percent_unemployed',"popdensity_sqkm", 'average_time_work', "percent_publictransport","pct_carcommutes_over45m","pct_ptcommutes_over45m",
                                       "delta_pct_ptcar_longcommute" )])
msa_ss<-merge(msa_ss, msa_acs_ss, by=c("CBSAFP","year"), 
              all.x=T, suffixes=c("","_acs"))

# (AIII) BEA data #####
#Downloaded from https://www.bea.gov/
BEA_GDP<-read.csv('datafiles_new/BEA_GDP_MSA.csv', stringsAsFactors = F)[1:385, c("GeoFips","GeoName",paste0("X", 2011:2018))]
##note:GDP numbers are in 1000, chained 2017 dollars. 
BEA_GDP<-BEA_GDP %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "year",
    names_prefix = "X",
    values_to = "GDP",
    values_drop_na = TRUE
  )

BEA_POP<-read.csv('datafiles_new/BEA_totalpop_MSA.csv', stringsAsFactors = F)[1:385, c("GeoFips","GeoName",paste0("X", 2011:2018))]
BEA_POP<-BEA_POP %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "year",
    names_prefix = "X",
    values_to = "BEA_Pop",
    values_drop_na = TRUE
  )

BEA_all<-merge(BEA_GDP,BEA_POP, by=c("GeoFips","year"))
BEA_all[c("GeoName.x","GeoName.y")]<-NULL
BEA_all$year<-as.numeric(BEA_all$year)
BEA_all$GDP_percapita<-BEA_all$GDP/BEA_all$BEA_Pop

msa_ss<-merge(msa_ss, BEA_all, by.x=c("CBSAFP", "year"), by.y=c("GeoFips","year"), all.x=T)
sum(is.na(msa_ss$GDP_percapita)) ## all perfectly joined, no zero

## (B) Cleaning up variables  ####
MSA_popprofiletokeep<-c("age55_older", "income1250_less",
                         "asian", "hispanic","female",
                        "education_college", "education_lesshighschool") 

AgriMiningUtilCons<-c("agriculture","mining","utilities","construction")
Manuf<-c("manufacturing")
TradeCommerceTransport<-c("wholesale","retail","transportation")
InfoFinProf_services<-c("information","finance","realestate","professional","management","waste")
PublicSocial_services<-c("education","healthcare", "publicadmin")
Entertainment<-c("arts","accommodation")
#Others<-"otherservices"

msa_ss_full<-msa_ss%>% mutate(across(all_of(MSA_popprofiletokeep), ~ . / total * 100, .names = "pct_workers_{.col}"))%>%
                               select(-all_of(MSA_popprofiletokeep))
                             
msa_ss_full<-msa_ss_full%>%
  mutate(
    AgriMiningUtilCons = rowSums(select(., all_of(AgriMiningUtilCons))),
    Manuf = rowSums(select(., all_of(Manuf))),
    TradeCommerceTransport = rowSums(select(., all_of(TradeCommerceTransport))),
    InfoFinProf_services = rowSums(select(., all_of(InfoFinProf_services))),
    PublicSocial_services = rowSums(select(., all_of(PublicSocial_services))),
    Entertainment = rowSums(select(., all_of(Entertainment)))
  )%>% mutate(
    pct_workers_AgriMiningUtilCons = AgriMiningUtilCons / total * 100,
    pct_workers_Manuf = Manuf / total * 100,
    pct_workers_TradeCommerceTransport = TradeCommerceTransport / total * 100,
    pct_workers_InfoFinProf_services = InfoFinProf_services / total * 100,
    pct_workers_PublicSocial_services = PublicSocial_services / total * 100,
    pct_workers_Entertainment = Entertainment / total * 100,
    pct_workers_Jobs_others = otherservices / total * 100
  )
#msa_ss_full$kvmt_per_hh_ami =msa_ss_full$vmt_per_hh_ami/1000
msa_ss_full$popdensity_1000persqkm =msa_ss_full$popdensity_sqkm/1000

msa_ss_full$resiworkrace_delta<-round(msa_ss_full$multigroupd_race_rac-msa_ss_full$multigroupd_race_wac,2) 
msa_ss_full$resiworkrace_delta_outlier<-ifelse(msa_ss_full$resiworkrace_delta<0,1,0)
table(msa_ss_full$resiworkrace_delta_outlier, msa_ss_full$year) #16 outliers in 2018

msa_ss_full$resiworkedu_delta<-round(msa_ss_full$multigroupd_education_rac-msa_ss_full$multigroupd_education_wac,2) 
msa_ss_full$resiworkedu_delta_outlier<-ifelse(msa_ss_full$resiworkedu_delta>0,1,0)
table(msa_ss_full$resiworkedu_delta_outlier, msa_ss_full$year) #5 outliers in 2018

msa_ss_full$resiworkethnicity_delta<-round(msa_ss_full$d_nonhispanichispanic_rac-msa_ss_full$d_nonhispanichispanic_wac,2) 
msa_ss_full$resiworkethnicity_delta_outlier<-ifelse(msa_ss_full$resiworkethnicity_delta<0,1,0)
table(msa_ss_full$resiworkethnicity_delta_outlier, msa_ss_full$year) #41 outliers in 2018

rm(county_2020_suburban, msa_rac,msa_wac, msa_acs)

### (C) DESCRIPTIVE STATS ######## 
source('R_code/label_dict.R')

MSA_poppctvariables<-c(unlist(lapply(c("age55_older", "income1250_less", "education_college", "female"
),function(x) paste0("pct_workers_",x))),"percent_suburban", "GDP_percapita") #"education_college", "education_lesshighschool" " choose either to avoid multicollinearity

transportvariables1<-c("t_ami", "kvmt_per_hh_ami" , "pct_transit_commuters_ami", "compact_ndx", "emp_ovrll_ndx", "Rank_transit_connectivity")
transportvariables2<-c("percent_publictransport","pct_ptcommutes_over45m", "pct_carcommutes_over45m","popdensity_1000persqkm") #"average_time_work",  remove because collinear with car commute factor, #"pct_ptcommutes_over45m",delta_pct_ptcar_longcommute

socialnetworkvariables<-c("percent_white", "delta_nonwhite_white_unemployed" ,"gini")
MSA_jobvariables<-c("percent_unemployed",
                    "pct_workers_AgriMiningUtilCons",
                    "pct_workers_Manuf", 
                    "pct_workers_TradeCommerceTransport", #,
                    "pct_workers_InfoFinProf_services",  
                    "pct_workers_PublicSocial_services",
                    "pct_workers_Entertainment")


allvariables<-c(transportvariables2, MSA_jobvariables, socialnetworkvariables, MSA_poppctvariables )
correlationsmatrix<-cor(subset(msa_ss_full, year==2018 &resiworkethnicity_delta_outlier==0)[allvariables], use = "pairwise.complete.obs")
colnames(correlationsmatrix) <- label_dict[colnames(correlationsmatrix)]
rownames(correlationsmatrix) <- label_dict[rownames(correlationsmatrix)]

# Following code generates Supplementary Annex 5, Figure 5.1, Correlation Plot
corrplot::corrplot(correlationsmatrix, method="number", tl.cex = 0.5, number.cex=0.7, number.digits=1, tl.col="black", bg="grey", type="lower")

# Following code generates Supplementary Annex 5, Table 5.1, Distribution of MSA-Level Characteristics
summaryfunction <- list(
  SRange = ~paste0(round(min(.x, na.rm = TRUE),2), ", ", round(max(.x, na.rm = TRUE),2)),
  SMean_sd= ~paste0(round(mean(.x, na.rm = TRUE),2)," (", round(sd(.x, na.rm = TRUE),2), ")")
)

summarydf<-msa_ss_full%>%group_by(year)%>%
  summarise( across(all_of(c("multigroupd_race_rac","multigroupd_race_wac", "multigroupd_education_rac", "multigroupd_education_wac",
                             transportvariables2, socialnetworkvariables, MSA_jobvariables,MSA_poppctvariables)), summaryfunction)) 

summarydf_t<-t(summarydf)
summarydf_t <-as.data.frame(cbind(rownames(summarydf_t),summarydf_t))
rownames(summarydf_t)<-NULL
colnames(summarydf_t)<-summarydf_t[1,]
summarydf_t<-summarydf_t[2:nrow(summarydf_t),]
summarydf_t[c("variable", "summary_stat")] <- do.call(rbind, strsplit(summarydf_t$year, "_S"))
summarydf_t$year<-NULL
summarydf_t$variable<-label_dict[summarydf_t$variable]
summarydf_t<-summarydf_t[c("variable", "summary_stat", "2011","2018")]
summarydf_t<-cbind(subset(summarydf_t, summary_stat=="Range"), subset(summarydf_t, summary_stat=="Mean_sd"))
summarydf_t<-summarydf_t[c(1,7,3, 8,4)]
colnames(summarydf_t)<-c("variable","2011 Mean(SD)","2011 Range","2018 Mean(SD)","2018 Range")

 ### (CI)  outlier analyses: 
msa_ss_full$resiworkrace_delta_neg<-ifelse(msa_ss_full$resiworkrace_delta<0,1,0)
msa_ss_full$resiworkedu_delta_pos<-ifelse(msa_ss_full$resiworkedu_delta>0,1,0)

msa_ss_full%>%group_by(resiworkrace_delta_neg, year)%>%summarise(n=n())
msa_ss_full%>%group_by(resiworkedu_delta_pos, year)%>%summarise(n=n())

summaryfunction <- list(
  SRange = ~paste0(round(min(.x, na.rm = TRUE),2), ", ", round(max(.x, na.rm = TRUE),2)),
  SMean_sd= ~paste0(round(mean(.x, na.rm = TRUE),2)," (", round(sd(.x, na.rm = TRUE),2), ")")
)

generate_outliersummaries<-function(outcome_variable){
  summarydf<-msa_ss_full%>%group_by(year, !!!syms(outcome_variable))%>%
    summarise(across(all_of(c("multigroupd_race_rac","multigroupd_race_wac", "multigroupd_education_rac", "multigroupd_education_wac",
                               transportvariables2, socialnetworkvariables, MSA_jobvariables,MSA_poppctvariables)), summaryfunction)) 
  summarydf_t<-t(summarydf)
  summarydf_t <-as.data.frame(cbind(rownames(summarydf_t),summarydf_t))
  rownames(summarydf_t)<-NULL
  colnames(summarydf_t)<-paste0(summarydf_t[1,],"_",summarydf_t[2,])
  summarydf_t<-summarydf_t[3:nrow(summarydf_t),]
  summarydf_t[c("variable", "summary_stat")] <- do.call(rbind, strsplit(summarydf_t$year, "_S"))
  summarydf_t$year<-NULL
  summarydf_t$variable<-as.character(label_dict[summarydf_t$variable])
  summarydf_t<-summarydf_t[c("variable", "summary_stat", "2018_0","2018_1")]
  summarydf_t<-subset(summarydf_t, summary_stat=="Mean_sd")[,c(1,3,4)]
  colnames(summarydf_t)<-c("Variable","Mean(SD) of Non-Outlier MSAs", "Mean(SD) of Outlier MSAs")
  summarydf_t
}
generate_outliersummariesallyears<-function(outcome_variable){
  summarydf<-msa_ss_full%>%group_by( !!!syms(outcome_variabl)%>%
    summarise(across(all_of(c("multigroupd_race_rac","multigroupd_race_wac", "multigroupd_education_rac", "multigroupd_education_wac",
                              transportvariables2, socialnetworkvariables, MSA_jobvariables,MSA_poppctvariables)), summaryfunction)) 
  summarydf_t<-t(summarydf)
  summarydf_t <-as.data.frame(summarydf_t)
  summarydf_t$variable<-rownames(summarydf_t)
  rownames(summarydf_t)<-NULL  #colnames(summarydf_t)<-paste0(summarydf_t[1,],"_",summarydf_t[2,])
  summarydf_t<-summarydf_t[2:nrow(summarydf_t),]
  summarydf_t[c("variable", "summary_stat")] <- do.call(rbind, strsplit(summarydf_t$variable, "_S"))
  #summarydf_t$year<-NULL
  summarydf_t$variable<-as.character(label_dict[summarydf_t$variable])
  summarydf_t<-subset(summarydf_t, summary_stat=="Mean_sd")[,c(3,1,2)]
  colnames(summarydf_t)<-c("Variable","Mean(SD) of Non-Outlier MSAs", "Mean(SD) of Outlier MSAs")
  summarydf_t
}
# Following code generates Supplementary Annex 2, Table 2.1, and Table 2.2
race_outlier<-generate_outliersummaries("resiworkrace_delta_neg")
subset(msa_ss_full,year==2018 & resiworkrace_delta_neg==1)$GeoName

class_outlier<-generate_outliersummaries("resiworkedu_delta_pos")
subset(msa_ss_full,year==2018 & resiworkedu_delta_pos==1)$GeoName

## 2FE REGRESSION analysis####
### 1. Single Moderator Variable 2FE regressions ######
## Following code generates Supplementary Annex 3, Figure 3.1 and Figure 3.2: Forest plots 

#set up dataframe
pdata <- pdata.frame(msa_ss_full, index = c("CBSAFP", "year"))
testing<-msa_ss_full[c("CBSAFP","year","total_pop", "popdensity_1000persqkm")]
testing$area<-round(testing$total_pop/testing$popdensity_1000persqkm)
testing<-testing%>%group_by(CBSAFP)%>%mutate(n=n_distinct(area)) ## essentially the two should be collinear

MSA_poppctvariables<-c(unlist(lapply(c("age55_older", "income1250_less", "education_lesshighschool", "female"
),function(x) paste0("pct_workers_",x))),"percent_suburban", "GDP_percapita") #"education_college", "education_lesshighschool" " choose either to avoid multicollinearity

transportvariables2<-c("percent_publictransport","pct_ptcommutes_over45m", "pct_carcommutes_over45m","popdensity_1000persqkm") #"average_time_work",  remove because collinear with car commute factor, #"pct_ptcommutes_over45m",delta_pct_ptcar_longcommute

socialnetworkvariables<-c("percent_white", "delta_nonwhite_white_unemployed" ,"gini")
MSA_jobvariables<-c("percent_unemployed",
                    "pct_workers_AgriMiningUtilCons",
                    "pct_workers_Manuf", 
                    "pct_workers_TradeCommerceTransport", #,
                    "pct_workers_InfoFinProf_services",  
                    "pct_workers_PublicSocial_services",
                    "pct_workers_Entertainment")

sigvars_race<-c()
sigvars_education<-c()

forestplot_race_list<-list()
forestplot_class_list<-list()

for(j in 1:3){
  chosen_variables<-list( transportvariables2, MSA_jobvariables, socialnetworkvariables )[[j]] #transportvariables1,
  suffix<-c( "transportvars2", "jobvars", "networkvars")[j] 
  forestplottitle<-c( "Transportation & Urban Form Moderators", "Labour Market Moderators", "Diversity & Inclusiveness Moderators")[j]
  resultslist_race<-list()
  resultslist_class<-list()
  resultslist_ethnic<-list()
  
  for(i in 1:length(chosen_variables)){
    testvar<-chosen_variables[i]
    formula_input1<-paste0("multigroupd_race_wac~multigroupd_race_rac*", testvar,"+", paste(MSA_poppctvariables,collapse="+"))
    formula_input2<-paste0("multigroupd_education_wac~multigroupd_education_rac*", testvar,"+", paste(MSA_poppctvariables,collapse="+"))
    formula_input3<-paste0("d_nonhispanichispanic_wac~d_nonhispanichispanic_rac*", testvar,"+", paste(MSA_poppctvariables,collapse="+"))
    
    m1<-plm(as.formula(formula_input1), data=pdata, model = "within", effect = "twoways")
    m2<-plm(as.formula(formula_input2), data=pdata, model = "within", effect = "twoways")
    m3<-plm(as.formula(formula_input3), data=pdata, model = "within", effect = "twoways")
    
    resultslist_race[[i]]<-m1
    resultslist_class[[i]]<-m2
    resultslist_ethnic[[i]]<-m3
  }
  preferredorderterms1<-c("multigroupd_race_rac", chosen_variables, unlist(lapply(chosen_variables, function(x) paste0("multigroupd_race_rac:", x))),
                          MSA_poppctvariables, "(Intercept)")
  variablelabels1<-unlist(lapply(preferredorderterms1, function(x) label_dict[x]))
  stargazer::stargazer(resultslist_race, type="html", single.row=T,
                       order= paste0("^", preferredorderterms1 , "$"),
                       covariate.labels = variablelabels1,
                       star.cutoffs=c(0.05, 0.01, 0.001),
                       dep.var.labels="Workplace Racial/Ethnic Segregation Score",
                       out=paste0("Regression_results/FE/workhome_racialseg_",suffix,"_FE.html"))
  
  preferredorderterms2<-c("multigroupd_education_rac", chosen_variables, unlist(lapply(chosen_variables, function(x) paste0("multigroupd_education_rac:", x))),
                          MSA_poppctvariables, "(Intercept)")
  variablelabels2<-unlist(lapply(preferredorderterms2, function(x) label_dict[x]))
  stargazer::stargazer(resultslist_class, type="html", single.row=T,
                       order= paste0("^", preferredorderterms2 , "$"),
                       covariate.labels = variablelabels2,
                       star.cutoffs=c(0.05, 0.01, 0.001),
                       dep.var.labels="Workplace Educational Segregation Score",
                       out=paste0("Regression_results/FE/workhome_educationalseg_",suffix,"_FE.html"))
  
  preferredorderterms3<-c("d_nonhispanichispanic_rac", chosen_variables, unlist(lapply(chosen_variables, function(x) paste0("multigroupd_education_rac:", x))),
                          MSA_poppctvariables, "(Intercept)")
  variablelabels3<-unlist(lapply(preferredorderterms3, function(x) label_dict[x]))
  stargazer::stargazer(resultslist_ethnic, type="html", single.row=T,
                       order= paste0("^", preferredorderterms3 , "$"),
                       covariate.labels = variablelabels3,
                       star.cutoffs=c(0.05, 0.01, 0.001),
                       dep.var.labels="Workplace Educational Segregation Score",
                       out=paste0("Regression_results/FE/workhome_ethniclseg_",suffix,"_FE.html"))

  
  resultslist_race_sigvars<-list()
  resultslist_class_sigvars<-list()

  for(i in 1:length(resultslist_race)){
    coefficients_table <- summary(resultslist_race[[i]])$coefficients
    significant_variables <- rownames(coefficients_table)[coefficients_table[, "Pr(>|t|)"] < 0.05]
    significant_variables<-significant_variables[significant_variables%in%c(chosen_variables, paste0("multigroupd_race_rac:",chosen_variables))]
    resultslist_race_sigvars[[i]]<-significant_variables}
  
  for(i in 1:length(resultslist_class)){
    coefficients_table <- summary(resultslist_class[[i]])$coefficients
    significant_variables <- rownames(coefficients_table)[coefficients_table[, "Pr(>|t|)"] < 0.05]
    significant_variables<-significant_variables[significant_variables%in%c(chosen_variables, paste0("multigroupd_education_rac:",chosen_variables))]
    resultslist_class_sigvars[[i]]<-significant_variables}
  
  sigvars_race<-c(sigvars_race,unlist(resultslist_race_sigvars))
  sigvars_education<-c(sigvars_education, unlist(resultslist_class_sigvars))
  
  name_vars_label<-unlist(lapply(c(chosen_variables, paste0("multigroupd_race_rac:",chosen_variables)), function(x) label_dict[x]))
  forestplot_race_list[[j]]<-plot_models(resultslist_race,  rm.terms=c("multigroupd_race_rac", MSA_poppctvariables), 
                                axis.labels=name_vars_label , p.threshold = c(0.05),colors="bw", show.legend=FALSE,
                                legend.title=NULL,title=forestplottitle,
                                show.values = FALSE, show.p = FALSE, p.shape = TRUE) + ylim(-0.5, 0.6)
  
  
  name_vars_label<-unlist(lapply(c(chosen_variables, paste0("multigroupd_education_rac:",chosen_variables)), function(x) label_dict[x]))
  forestplot_class_list[[j]]<-plot_models(resultslist_class,  rm.terms=c("multigroupd_education_rac", MSA_poppctvariables), 
                                axis.labels=name_vars_label , p.threshold = c(0.05),colors="bw", show.legend=FALSE,
                                legend.title=NULL,title=forestplottitle,
                                show.values = FALSE, show.p = FALSE, p.shape = TRUE)+ ylim(-0.5, 0.5)

}

title <- textGrob("Class Segregation Regression Models \n testing Single Moderator Variable", 
                  gp = gpar(fontsize = 14, fontface = "bold"))
combined_plot <- grid.arrange(
  title,
  arrangeGrob(forestplot_class_list[[1]], forestplot_class_list[[2]],forestplot_class_list[[3]],ncol = 1),
  ncol = 1,
  heights = c(0.1, 1)  # Adjust the heights as needed
)
ggsave("Regression_results/Plots/class_bivariate_moderators.png", plot = combined_plot, width = 6, height = 16)

title <- textGrob("Racial Segregation Regression Models\n testing Single Moderator Variable", 
                  gp = gpar(fontsize = 14, fontface = "bold"))
combined_plot <- grid.arrange(
  title,
  arrangeGrob(forestplot_race_list[[1]], forestplot_race_list[[2]],forestplot_race_list[[3]],ncol = 1),
  ncol = 1,
  heights = c(0.1, 1)  # Adjust the heights as needed
)
ggsave("Regression_results/Plots/race_bivariate_moderators.png", plot = combined_plot, width = 6, height = 16)

plot_models(resultslist_class,  rm.terms=c("multigroupd_education_rac", MSA_poppctvariables), 
            axis.labels=name_vars_label , p.threshold = c(0.05),colors="bw", show.legend=TRUE,
            legend.title=NULL,title=forestplottitle,
            show.values = FALSE, show.p = FALSE, p.shape = TRUE)+ ylim(-0.5, 0.5)

rm(forestplot_class_list,forestplot_race_list,combined_plot)

### run combined FE model FOR RACE ####

### Following code generates Supplementary Annex 3, Table 3.1: Two-way fixed effects regression results, showing the associations between workhood racial segregation, residential racial segregation, potential moderators and other covariates
MSA_poppctvariables<-c(unlist(lapply(c("age55_older", "income1250_less", "education_college", "female"),
                                     function(x) paste0("pct_workers_",x))),"percent_suburban", "GDP_percapita") #"education_college", "education_lesshighschool" " choose either to avoid multicollinearity
chosen_interactvariables_all <- sub(".*:", "", sigvars_race[grepl(":", sigvars_race)])
chosen_mainvariables_all<-unique(c(sigvars_race[!grepl(":", sigvars_race)],chosen_interactvariables_all))

formula_input<-paste0("multigroupd_race_wac~", 
                      paste( paste("multigroupd_race_rac*",chosen_interactvariables_all, sep=""), collapse="+"),"+", paste(chosen_mainvariables_all, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m1<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")

coefficients_table <- summary(m1)$coefficients
# Extract the names of the variables that are significant
significant_variables <- rownames(coefficients_table)[coefficients_table[, "Pr(>|t|)"] < 0.05]
significant_variables<-significant_variables[significant_variables%in%c(chosen_mainvariables_all, paste0("multigroupd_race_rac:",chosen_interactvariables_all))]

chosen_interactvariables_ss <- sub(".*:", "", significant_variables[grepl(":", significant_variables)])
chosen_mainvariables_ss<-unique(c(significant_variables[!grepl(":", significant_variables)],chosen_interactvariables_ss))

formula_input<-paste0("multigroupd_race_wac~", 
                      paste( paste("multigroupd_race_rac*",chosen_interactvariables_ss, sep=""), collapse="+"),"+", paste(chosen_mainvariables_ss, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m2<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")

formula_input<-paste0("multigroupd_race_wac~multigroupd_education_rac+", 
                      paste( paste("multigroupd_race_rac*",chosen_interactvariables_ss, sep=""), collapse="+"),"+", paste(chosen_mainvariables_ss, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m3<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")

formula_input<-paste0("multigroupd_race_wac~multigroupd_education_wac+", 
                      paste( paste("multigroupd_race_rac*",chosen_interactvariables_ss, sep=""), collapse="+"),"+", paste(chosen_mainvariables_ss, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m4<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")

preferredordertermFULL<-c("multigroupd_race_rac", "multigroupd_education_rac", "multigroupd_education_wac",
                          transportvariables2, paste0("multigroupd_race_rac:", transportvariables2),
                          MSA_jobvariables, paste0("multigroupd_race_rac:", MSA_jobvariables),
                          socialnetworkvariables, paste0("multigroupd_race_rac:", socialnetworkvariables),
                          MSA_poppctvariables )

preferredorderterms<-c("multigroupd_race_rac", "multigroupd_education_rac", "multigroupd_education_wac",
                       unique(c(chosen_mainvariables_all,chosen_interactvariables_all)), 
                       unlist(lapply(chosen_interactvariables_all, function(x) paste0("multigroupd_race_rac:", x))),
                       MSA_poppctvariables)

preferredorderterms<-preferredordertermFULL[preferredordertermFULL%in%preferredorderterms]
variablelabels<-unlist(lapply(preferredorderterms, function(x) label_dict[x]))
stargazer::stargazer(m1, m2, m3, m4, type="html", single.row=T,
                     order= paste0("^", preferredorderterms , "$"),
                     covariate.labels = variablelabels,
                     dep.var.labels="Workplace Racial Segregation Score",
                     star.cutoffs=c(0.05, 0.01, 0.001),
                     out=paste0("Regression_results/FE/workhome_racialseg_COMBINED.html"))

## robustness with acs data ##
# Following code generates Supplementary Annex 4, Table 4.3 Regression Models using ACS 5 year survey derived estimates of residential racial segregation 
formula_input<-paste0("multigroupd_race_wac~", 
                      paste( paste("multigroupd_race_rac_acs*",chosen_interactvariables_ss, sep=""), collapse="+"),"+", paste(chosen_mainvariables_ss, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m5<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")

formula_input<-paste0("multigroupd_race_wac~multigroupd_education_rac+", 
                      paste( paste("multigroupd_race_rac_acs*",chosen_interactvariables_ss, sep=""), collapse="+"),"+", paste(chosen_mainvariables_ss, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m6<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")

formula_input<-paste0("multigroupd_race_wac~multigroupd_education_wac+", 
                      paste( paste("multigroupd_race_rac_acs*",chosen_interactvariables_ss, sep=""), collapse="+"),"+", paste(chosen_mainvariables_ss, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m7<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")

preferredordertermFULL<-c("multigroupd_race_rac_acs", "multigroupd_education_rac", "multigroupd_education_wac",
                          transportvariables2, paste0("multigroupd_race_rac_acs:", transportvariables2),
                          MSA_jobvariables, paste0("multigroupd_race_rac_acs:", MSA_jobvariables),
                          socialnetworkvariables, paste0("multigroupd_race_rac_acs:", socialnetworkvariables),
                          MSA_poppctvariables )

preferredorderterms<-c("multigroupd_race_rac_acs", "multigroupd_education_rac", "multigroupd_education_wac",
                       chosen_mainvariables_ss, 
                       unlist(lapply(chosen_interactvariables_ss, function(x) paste0("multigroupd_race_rac_acs:", x))),
                       MSA_poppctvariables)
preferredorderterms<-preferredordertermFULL[preferredordertermFULL%in%preferredorderterms]

variablelabels<-unlist(lapply(preferredorderterms, function(x) label_dict[x]))
stargazer::stargazer(m5, m6, m7, type="html", single.row=T,
                     order= paste0("^", preferredorderterms , "$"),
                     covariate.labels = variablelabels,
                     dep.var.labels="Workplace Racial Segregation Score",
                     star.cutoffs=c(0.05, 0.01, 0.001),
                     out=paste0("Regression_results/FE/workhome_racialseg_COMBINED_ACS.html"))

##robustness with density data
## Following code generates Supplementary Annex 4, Table 4.1: Two-way fixed effects regression results from models using alternative measure of density  
chosen_mainvariables_ss2<-c("pct_workers_AgriMiningUtilCons","percent_white",
                           "all_bua","pct_workers_TradeCommerceTransport")
chosen_interactvariables_ss2<-c("all_bua","pct_workers_TradeCommerceTransport","percent_white")
formula_input<-paste0("multigroupd_race_wac~", 
                      paste( paste("multigroupd_race_rac*",chosen_interactvariables_ss2, sep=""), 
                             collapse="+"),"+", paste(chosen_mainvariables_ss2, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m5<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")
formula_input<-paste0("multigroupd_race_wac~multigroupd_education_rac+", 
                      paste( paste("multigroupd_race_rac*",chosen_interactvariables_ss2, sep=""), 
                             collapse="+"),"+", paste(chosen_mainvariables_ss2, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m6<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")

formula_input<-paste0("multigroupd_race_wac~multigroupd_education_wac+", 
                      paste( paste("multigroupd_race_rac*",chosen_interactvariables_ss2, sep=""), 
                             collapse="+"),"+", paste(chosen_mainvariables_ss2, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m7<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")

preferredordertermFULL<-c("multigroupd_race_rac", "multigroupd_education_rac", "multigroupd_education_wac",
                          transportvariables2, "all_bua", paste0("multigroupd_race_rac:", transportvariables2), 
                          "multigroupd_race_rac:all_bua",
                          MSA_jobvariables, paste0("multigroupd_race_rac:", MSA_jobvariables),
                          socialnetworkvariables, paste0("multigroupd_race_rac:", socialnetworkvariables),
                          MSA_poppctvariables )

preferredorderterms<-c("multigroupd_race_rac", "multigroupd_education_rac", "multigroupd_education_wac",
                       chosen_mainvariables_ss2, 
                       unlist(lapply(chosen_interactvariables_ss2, function(x) paste0("multigroupd_race_rac:", x))),
                       MSA_poppctvariables)
preferredorderterms<-preferredordertermFULL[preferredordertermFULL%in%preferredorderterms]

variablelabels<-unlist(lapply(preferredorderterms, function(x) label_dict[x]))
stargazer::stargazer(m5, m6, m7, type="html", single.row=T,
                     order= paste0("^", preferredorderterms , "$"),
                     covariate.labels = variablelabels,
                     dep.var.labels="Workplace Racial Segregation Score",
                     star.cutoffs=c(0.05, 0.01, 0.001),
                     out=paste0("Regression_results/FE/workhome_racialseg_TransportRobustness.html"))

### Following code generates Supplementary Annex 3, Figure 3.2: Interaction Effects for Significant Moderators of the relationship between workhood and neighborhood racial segregation
wut<-unlist(lapply(chosen_interactvariables_ss, function(x) label_dict[x])) 
huh <- msa_ss_full %>%
  filter(year %in% c(2011, 2018)) %>%
  group_by(CBSAFP) %>%
  arrange(year) %>%
  summarize(across(all_of(c(chosen_interactvariables_ss )), ~ .[year == 2018] - .[year == 2011])) %>%
  ungroup()%>%
  summarize(across(all_of(c(chosen_interactvariables_ss)), 
                   list(mean = ~ round(mean(.),2), max = ~ round(max(.),2), min = ~ round(min(.),2))))  

n<-length(chosen_interactvariables_ss)
plotlist<-list()
for(i in 1:n){
  averageresiseg=round(mean(msa_ss_full[msa_ss_full$year==2011,"multigroupd_race_rac"]),2)
  summaryvalues<-huh[c(1,2,3)+(i-1)*3]
  plotlist[[i]]<-plot_model(m2, type = "pred", terms = c("multigroupd_race_rac",
                                                         paste0(chosen_interactvariables_ss[i]," [", paste(summaryvalues, collapse=","),"]")),
                            legend.title="Min, Mean, Max",
                            axis.title=c("Residential Racial Segregation","Workhood Racial Segregation")) +
    geom_vline(xintercept=averageresiseg, linetype = "dashed", color = "darkgrey", linewidth=0.5)+
    annotate("text", x = averageresiseg-0.12, y = Inf, label = "Mean Resi Seg, 2011", vjust = 1.5, color = "black", size=3)+
    ggplot2::labs(title = paste0("Interaction with ",wut[[i]]))
  
}
title <- textGrob("Plotting Interactions with Residential Segregation", gp = gpar(fontsize = 20, fontface = "bold"))
combined_plot <- grid.arrange(
  title,
  arrangeGrob(plotlist[[1]], plotlist[[2]],plotlist[[3]],ncol = 3),
  ncol = 1,
  heights = c(0.1, 1)  # Adjust the heights as needed
)
#pull out m4 results for the final model for main paper
m4race<-m4

### run combined FE model FOR EDUCATION
### Following code generates Supplementary Annex 3, Table 3.2: Two-way fixed effects regression results, showing the associations between workhood class segregation, residential racial segregation, potential moderators and other covariates

MSA_poppctvariables<-c(unlist(lapply(c("age55_older", "income1250_less", "education_college", "female"),
                                     function(x) paste0("pct_workers_",x))),"percent_suburban", "GDP_percapita") #"education_college", "education_lesshighschool" " choose either to avoid multicollinearity
chosen_interactvariables_all <- sub(".*:", "", sigvars_education[grepl(":", sigvars_education)])
chosen_mainvariables_all<-unique(c(sigvars_education[!grepl(":", sigvars_education)],chosen_interactvariables_all))

formula_input<-paste0("multigroupd_education_wac~", 
                      paste( paste("multigroupd_education_rac*",chosen_interactvariables_all, sep=""), collapse="+"),"+", 
                      paste(chosen_mainvariables_all, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m1<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")

coefficients_table <- summary(m1)$coefficients
# Extract the names of the variables that are significant
significant_variables <- rownames(coefficients_table)[coefficients_table[, "Pr(>|t|)"] < 0.05]
significant_variables<-significant_variables[significant_variables%in%c(chosen_mainvariables_all, paste0("multigroupd_education_rac:",chosen_interactvariables_all))]

chosen_interactvariables_ss <- sub(".*:", "", significant_variables[grepl(":", significant_variables)])
chosen_mainvariables_ss<-unique(c(significant_variables[!grepl(":", significant_variables)],chosen_interactvariables_ss))

formula_input<-paste0("multigroupd_education_wac~", 
                      paste( paste("multigroupd_education_rac*",chosen_interactvariables_ss, sep=""), collapse="+"),"+", paste(chosen_mainvariables_ss, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m2<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")

formula_input<-paste0("multigroupd_education_wac~multigroupd_race_rac+", 
                      paste( paste("multigroupd_education_rac*",chosen_interactvariables_ss, sep=""), collapse="+"),"+", paste(chosen_mainvariables_ss, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m3<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")

formula_input<-paste0("multigroupd_education_wac~multigroupd_race_wac+", 
                      paste( paste("multigroupd_education_rac*",chosen_interactvariables_ss, sep=""), collapse="+"),"+", paste(chosen_mainvariables_ss, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m4<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")

preferredordertermFULL<-c("multigroupd_education_rac", "multigroupd_race_rac", "multigroupd_race_wac",
                          transportvariables2, paste0("multigroupd_education_rac:", transportvariables2),
                          MSA_jobvariables, paste0("multigroupd_education_rac:", MSA_jobvariables),
                          socialnetworkvariables, paste0("multigroupd_education_rac:", socialnetworkvariables),
                          MSA_poppctvariables )

preferredorderterms<-c("multigroupd_education_rac", "multigroupd_race_rac", "multigroupd_race_wac",
                       chosen_mainvariables_all, 
                       unlist(lapply(chosen_interactvariables_all, function(x) paste0("multigroupd_education_rac:", x))),
                       MSA_poppctvariables)
preferredorderterms<-preferredordertermFULL[preferredordertermFULL%in%preferredorderterms]

variablelabels<-unlist(lapply(preferredorderterms, function(x) label_dict[x]))
stargazer::stargazer(m1, m2, m3, m4, type="html", single.row=T,
                     order= paste0("^", preferredorderterms , "$"),
                     covariate.labels = variablelabels,
                     dep.var.labels="Workplace Educational Class Segregation Score",
                     star.cutoffs=c(0.05, 0.01, 0.001),
                     out=paste0("Regression_results/FE/workhome_educationseg_COMBINED.html"))

## robustness with acs data
# Following code generates Supplementary Annex 4, Table 4.4 Regression Models using ACS 5 year survey derived estimates of residential class segregation 
formula_input<-paste0("multigroupd_education_wac~", 
                      paste( paste("multigroupd_education_rac_acs*",chosen_interactvariables_ss, sep=""), collapse="+"),"+", paste(chosen_mainvariables_ss, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m5<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")

formula_input<-paste0("multigroupd_education_wac~multigroupd_race_rac+", 
                      paste( paste("multigroupd_education_rac_acs*",chosen_interactvariables_ss, sep=""), collapse="+"),"+", paste(chosen_mainvariables_ss, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m6<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")

formula_input<-paste0("multigroupd_education_wac~multigroupd_race_wac+", 
                      paste( paste("multigroupd_education_rac_acs*",chosen_interactvariables_ss, sep=""), collapse="+"),"+", paste(chosen_mainvariables_ss, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m7<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")

preferredordertermFULL<-c("multigroupd_education_rac_acs", "multigroupd_race_rac", "multigroupd_race_wac",
                          transportvariables2, paste0("multigroupd_education_rac_acs:", transportvariables2),
                          MSA_jobvariables, paste0("multigroupd_education_rac_acs:", MSA_jobvariables),
                          socialnetworkvariables, paste0("multigroupd_education_rac_acs:", socialnetworkvariables),
                          MSA_poppctvariables )

preferredorderterms<-c("multigroupd_education_rac_acs", "multigroupd_race_rac", "multigroupd_race_wac",
                       chosen_mainvariables_ss, 
                       unlist(lapply(chosen_interactvariables_ss, function(x) paste0("multigroupd_education_rac_acs:", x))),
                       MSA_poppctvariables)

preferredorderterms<-preferredordertermFULL[preferredordertermFULL%in%preferredorderterms]

variablelabels<-unlist(lapply(preferredorderterms, function(x) label_dict[x]))
stargazer::stargazer(m5, m6, m7, type="html", single.row=T,
                     order= paste0("^", preferredorderterms , "$"),
                     covariate.labels = variablelabels,
                     dep.var.labels="Workhood Educational Class Segregation Score",
                     star.cutoffs=c(0.05, 0.01, 0.001),
                     out=paste0("Regression_results/FE/workhome_educationseg_COMBINED_ACS.html"))

                      
## robustness with transit variables 
# Following code generates Supplementary Annex 4, Table 4.2 Two-way fixed effects regression results from models using alternative measure of transportation infrastructure provision 
chosen_mainvariables_ss2<-c("percent_publictransport","distance_k","pct_workers_Entertainment","all_bua","delta_nonwhite_white_unemployed")
chosen_interactvariables_ss2<-c("percent_publictransport", "distance_k", "pct_workers_Entertainment")


formula_input<-paste0("multigroupd_education_wac~", 
                      paste( paste("multigroupd_education_rac*",chosen_interactvariables_ss2, sep=""), 
                             collapse="+"),"+", paste(chosen_mainvariables_ss2, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m5<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")

formula_input<-paste0("multigroupd_education_wac~multigroupd_race_rac+", 
                      paste( paste("multigroupd_education_rac*",chosen_interactvariables_ss2, sep=""), 
                             collapse="+"),"+", paste(chosen_mainvariables_ss2, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m6<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")

formula_input<-paste0("multigroupd_education_wac~multigroupd_race_wac+", 
                      paste( paste("multigroupd_education_rac*",chosen_interactvariables_ss2, sep=""), 
                             collapse="+"),"+", paste(chosen_mainvariables_ss2, collapse="+"),"+",
                      paste(MSA_poppctvariables, collapse="+"))
m7<-plm(as.formula(formula_input), data=pdata, model = "within", effect = "twoways")

preferredordertermFULL<-c("multigroupd_education_rac", "multigroupd_race_rac", "multigroupd_race_wac",
                          transportvariables2, "all_bua", "distance_k", "mean_local_gridness", paste0("multigroupd_education_rac:", transportvariables2), 
                          "multigroupd_education_rac:distance_k", "multigroupd_education_rac:mean_local_gridness",
                          MSA_jobvariables, paste0("multigroupd_education_rac:", MSA_jobvariables),
                          socialnetworkvariables, paste0("multigroupd_education_rac:", socialnetworkvariables),
                          MSA_poppctvariables )

preferredorderterms<-c("multigroupd_education_rac", "multigroupd_race_rac", "multigroupd_race_wac",
                       chosen_mainvariables_ss2, 
                       unlist(lapply(chosen_interactvariables_ss2, function(x) paste0("multigroupd_education_rac:", x))),
                       MSA_poppctvariables)
preferredorderterms<-preferredordertermFULL[preferredordertermFULL%in%preferredorderterms]

variablelabels<-unlist(lapply(preferredorderterms, function(x) label_dict[x]))
stargazer::stargazer(m5, m6, m7, type="html", single.row=T,
                     order= paste0("^", preferredorderterms , "$"),
                     covariate.labels = variablelabels,
                     dep.var.labels="Workplace Racial Segregation Score",
                     star.cutoffs=c(0.05, 0.01, 0.001),
                     out=paste0("Regression_results/FE/workhome_educationseg_TransportRobustness.html"))

### Following code generates Supplementary Annex, Figure 3.3: Interaction Effects for Significant Moderators of the relationship between workhood and neighborhood class segregation
huh <- msa_ss_full %>%
  filter(year %in% c(2011, 2018)) %>%
  group_by(CBSAFP) %>%
  arrange(year) %>%
  summarize(across(all_of(c(chosen_interactvariables_ss,"multigroupd_education_rac")), ~ .[year == 2018] - .[year == 2011])) %>%
  ungroup()%>%
  summarize(across(all_of(c(chosen_interactvariables_ss,"multigroupd_education_rac")), 
                   list(mean = ~ round(mean(.),2), max = ~ round(max(.),2), min = ~ round(min(.),2))))  

wut<-unlist(lapply(chosen_interactvariables_ss, function(x) label_dict[x])) 
n<-length(chosen_interactvariables_ss)
plotlist<-list()
for(i in 1:n){
  summaryvalues<-huh[c(1,2,3)+(i-1)*3]
  averageresiseg=round(mean(msa_ss_full[msa_ss_full$year==2011,"multigroupd_education_rac"]),2)
  plotlist[[i]]<-plot_model(m2, type = "pred", terms = c("multigroupd_education_rac",
                                                         paste0(chosen_interactvariables_ss[i]," [", paste(summaryvalues, collapse=","),"]")),
                            legend.title="Min, Mean, Max",
                            axis.title=c("Residential Class Segregation","Workhood Class Segregation"))+
    geom_vline(xintercept=averageresiseg, linetype = "dashed", color = "darkgrey", linewidth=0.5)+
    annotate("text", x = averageresiseg-0.03, y = Inf, label = "Mean Resi Seg, 2011", vjust = 1.5, color = "black", size=3)+
    ggplot2::labs(title = paste0("Interaction with ",wut[[i]]))

}
title <- textGrob("Plotting Interactions with Residential Segregation", gp = gpar(fontsize = 20, fontface = "bold"))
combined_plot <- grid.arrange(
  title,
  arrangeGrob(plotlist[[1]], plotlist[[2]],plotlist[[3]],plotlist[[4]],ncol = 2),
  ncol = 1,
  heights = c(0.1, 1)  # Adjust the heights as needed
)

### Following code generates Table 1,  Main Paper #####
preferredordertermFULL<-c( "multigroupd_race_rac","multigroupd_education_wac", "multigroupd_education_rac", "multigroupd_race_wac",
                            transportvariables2, paste0("multigroupd_race_rac:", transportvariables2), paste0("multigroupd_education_rac:", transportvariables2),
                            MSA_jobvariables, paste0("multigroupd_race_rac:", MSA_jobvariables), paste0("multigroupd_education_rac:", MSA_jobvariables),
                            socialnetworkvariables,paste0("multigroupd_race_rac:", socialnetworkvariables), paste0("multigroupd_education_rac:", socialnetworkvariables),
                            MSA_poppctvariables )

preferredorderterms<-unique(c(names(coefficients(m4race)), names(coefficients(m4))))
preferredorderterms<-preferredordertermFULL[preferredordertermFULL%in%preferredorderterms]

variablelabels<-unlist(lapply(preferredorderterms, function(x) label_dict[x]))

stargazer::stargazer(m4race, m4, type="html", single.row=T,
                     order= paste0("^", preferredorderterms , "$"),
                     covariate.labels = variablelabels,
                     dep.var.labels="Workplace Segregation Score",
                     star.cutoffs=c(0.05, 0.01, 0.001),
                     out=paste0("Regression_results/FE/workhome_class_race_finaltable.html"))

