#Diagnostic performance and clinical implications of rapid SARS-CoV-2 antigen testing in Mexico using real-world nationwide COVID-19 registry data
#Data Analysis: Omar Yaxmehen Bello-Chavolla (oyaxbell@yahoo.com.mx); Neftali Eduardo Antonio-Villa (nefoantonio@hotmail.com)
#Latest version of Analysis 01-Jan-2021
#Any question regarding analysis contact Omar Yaxmehen Bello-Chavolla

#### Database management####
library(readr); library(tidyverse); library(survival); library(mediation); library(ggpubr); library(rms); library(psych); library(smoothHR);library(survey); library(caret)
library(survminer); library(haven); library(rsq); library(ResourceSelection); library(ggsci);library(timereg); library(coxme); library(MatchIt); library(timeROC);library(officer)
library(pROC);librafactor(covid, labels=c("SARS-CoV-2 (-), SARS-CoV-2 (+)"))ry(sf); library(rgdal); library(viridis); library(ggsci); library(ggmap); library(scales); library(jtools); library(cowplot);library(npROCRegression)
library(ggstance); library(flextable); library(simPH); library(ggthemes); library(lme4); library(lmerTest); library(prismatic); library(readxl); library(zoo); library(OptimalCutpoints)
library(ggstatsplot);library(Rmpfr);library(ggrepel)
setwd("~/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/COVID - Antígeno")

#file<-paste0(if(format(Sys.time(), "%H:%M:%S")>"19:00:00"){
#  print(format(Sys.Date(),"%y%m%d"))
#  } else {
#    print(format(Sys.Date()-1,"%y%m%d"))
#    },"COVID19MEXICO.csv")

temp <- tempfile()
#download.file("http://datosabiertos.salud.gob.mx/gobmx/salud/datos_abiertos/datos_abiertos_covid19.zip",temp)
download.file("http://datosabiertos.salud.gob.mx/gobmx/salud/datos_abiertos/historicos/12/datos_abiertos_covid19_31.12.2020.zip",temp)
#covid <- read.csv(unz(temp, file),header=TRUE,sep = ",", encoding = "UTF-8")
covid <- read.csv(unz(temp, "201231COVID19MEXICO.csv"),header=TRUE,sep = ",", encoding = "UTF-8")
unlink(temp)

covid$id<-paste0(str_pad(covid$ENTIDAD_RES, 2,pad = "0"),str_pad(covid$MUNICIPIO_RES,3, pad="0"))
covid1<-covid[,c(14:15,18:30,32:35)]
covid<-covid[,-c(14:15,18:30,32:35)]
covid1[covid1==2]<-0
covid1[covid1==97]<-NA;covid1[covid1==98]<-NA;covid1[covid1==99]<-NA
covid<-as.data.frame(cbind(covid, covid1))
covid$EMBARAZO[is.na(covid$EMBARAZO)]<-0
covid$FECHA_DEF[covid$FECHA_DEF=="9999-99-99"]<-NA
covid$FECHA_DEF<-as.Date(covid$FECHA_DEF)
covid$TIPO_PACIENTE[covid$TIPO_PACIENTE==1]<-0;covid$TIPO_PACIENTE[covid$TIPO_PACIENTE==2]<-1
covid$muestra<-NULL; covid$muestra[covid$RESULTADO==1]<-1;covid$muestra[covid$RESULTADO!=1]<-0
covid$edad60<-NULL;covid$edad60[covid$EDAD>=60]<-1;covid$edad60[covid$EDAD<60]<-0
covid$edad40<-NULL;covid$edad40[covid$EDAD>=40]<-0;covid$edad40[covid$EDAD<40]<-1
covid$diabetes_40<-NULL;covid$diabetes_40[covid$DIABETES==1 & covid$edad40==1]<-1;covid$diabetes_40[covid$DIABETES!=1 & covid$edad40!=1]<-0
covid$Mortalidad<-NULL; covid$Mortalidad[is.na(covid$FECHA_DEF)]<-0;covid$Mortalidad[is.na(covid$FECHA_DEF)==FALSE]<-1
covid$FECHA_DEF[is.na(covid$FECHA_DEF)]<-as.Date(covid$FECHA_SINTOMAS)+30
covid$FU_time<-as.numeric(as.Date(covid$FECHA_DEF)-as.Date(covid$FECHA_SINTOMAS))
covid$FU_time[covid$FU_time<0]<-1
covid$Latencia<-as.numeric(as.Date(covid$FECHA_INGRESO)-as.Date(covid$FECHA_SINTOMAS))
covid$comorb<-covid$DIABETES+covid$OBESIDAD+covid$EPOC+covid$ASMA+covid$INMUSUPR+covid$HIPERTENSION+
  covid$CARDIOVASCULAR+covid$RENAL_CRONICA+covid$TABAQUISMO+covid$OTRA_COM
covid$comorb_d<-covid$OBESIDAD+covid$EPOC+covid$ASMA+covid$INMUSUPR+covid$HIPERTENSION+
  covid$CARDIOVASCULAR+covid$RENAL_CRONICA+covid$TABAQUISMO+covid$OTRA_COM
covid$comorb_ob<-covid$DIABETES+covid$EPOC+covid$ASMA+covid$INMUSUPR+covid$HIPERTENSION+
  covid$CARDIOVASCULAR+covid$RENAL_CRONICA+covid$TABAQUISMO+covid$OTRA_COM
covid$comorb_ckd<-covid$DIABETES+covid$OBESIDAD+covid$EPOC+covid$ASMA+covid$INMUSUPR+covid$HIPERTENSION+
  covid$CARDIOVASCULAR+covid$TABAQUISMO+covid$OTRA_COM
covid$comorb_cvd<-covid$DIABETES+covid$OBESIDAD+covid$EPOC+covid$ASMA+covid$INMUSUPR+covid$HIPERTENSION+
  covid$RENAL_CRONICA+covid$TABAQUISMO+covid$OTRA_COM
covid$comorb_copd<-covid$DIABETES+covid$OBESIDAD+covid$ASMA+covid$INMUSUPR+covid$HIPERTENSION+
  covid$CARDIOVASCULAR+covid$RENAL_CRONICA+covid$TABAQUISMO+covid$OTRA_COM
covid$comorb_asma<-covid$DIABETES+covid$OBESIDAD+covid$EPOC+covid$INMUSUPR+covid$HIPERTENSION+
  covid$CARDIOVASCULAR+covid$RENAL_CRONICA+covid$TABAQUISMO+covid$OTRA_COM
covid$comorb_immun<-covid$DIABETES+covid$OBESIDAD+covid$EPOC+covid$ASMA+covid$HIPERTENSION+
  covid$CARDIOVASCULAR+covid$RENAL_CRONICA+covid$TABAQUISMO+covid$OTRA_COM
covid$comorb_has<-covid$DIABETES+covid$OBESIDAD+covid$EPOC+covid$ASMA+covid$INMUSUPR+
  covid$CARDIOVASCULAR+covid$RENAL_CRONICA+covid$TABAQUISMO+covid$OTRA_COM
covid$comorb_tab<-covid$DIABETES+covid$OBESIDAD+covid$EPOC+covid$ASMA+covid$INMUSUPR+covid$HIPERTENSION+
  covid$CARDIOVASCULAR+covid$RENAL_CRONICA+covid$OTRA_COM

covid$comorb_dic[covid$comorb>0]<-1;covid$comorb_dic[covid$comorb==0]<-0
covid$diab_ob<-2*covid$DIABETES+covid$OBESIDAD
covid$HABLA_LENGUA_INDIG[is.na(covid$HABLA_LENGUA_INDIG)]<-0
covid$priv[covid$SECTOR==9]<-1; covid$priv[!covid$SECTOR==9]<-0
covid<- covid %>% mutate(comorb_cat = comorb %>% 
                           cut(c(0, 1, 2, Inf)))

comorb_levels <- paste0(c("0", "1", ">2"))

covid <- covid %>% 
  mutate(comorb_cat = comorb_cat %>% lvls_revalue(comorb_levels))

covid_et<-covid%>%filter(HABLA_LENGUA_INDIG==1, covid==1)

covid<- covid %>% mutate(age_cat = EDAD %>% 
                           cut(c(-Inf, 17, 30, 50, 70,Inf)))

age_levels <- paste0(c("0-17", "18-30", "31-50","51-70",">70"))

covid <- covid %>% 
  mutate(age_cat = age_cat %>% lvls_revalue(age_levels))


covid$comorb_d_dic<-2*(covid$comorb_d>0)+covid$DIABETES
covid$comorb_ob_dic<-2*(covid$comorb_ob>0)+covid$OBESIDAD
covid$comorb_ckd_dic<-2*(covid$comorb_ckd>0)+covid$RENAL_CRONICA
covid$comorb_has_dic<-2*(covid$comorb_has>0)+covid$HIPERTENSION
covid$comorb_asma_dic<-2*(covid$comorb_asma>0)+covid$ASMA
covid$comorb_cvd_dic<-2*(covid$comorb_cvd>0)+covid$CARDIOVASCULAR
covid$comorb_copd_dic<-2*(covid$comorb_copd>0)+covid$EPOC
covid$comorb_tab_dic<-2*(covid$comorb_tab>0)+covid$TABAQUISMO
covid$comodb_immun_dic<-2*(covid$comorb_immun>0)+covid$INMUSUPR
covid$tiempo<-as.numeric(as.Date(covid$FECHA_INGRESO)-as.Date(covid$FECHA_SINTOMAS))

load("df_mx.rda")
df_mx$ENTIDAD_UM<-df_mx$CLAVE_ENT
pop_mun<-df_mx %>% filter(AÑO==2020)%>%group_by(ENTIDAD_UM)%>%summarise(pop=sum(POB))


covid$muestra<-ifelse(covid$CLASIFICACION_FINAL==3, 1, 0)
covid<-covid %>% left_join(pop_mun, by="ENTIDAD_UM")

covid$muestra<-covid$TOMA_MUESTRA_ANTIGENO+2*covid$TOMA_MUESTRA_LAB
covid$FECHA_INGRESO<-as.Date(covid$FECHA_INGRESO)
covid$FECHA_SINTOMAS<-as.Date(covid$FECHA_SINTOMAS)
covid$covid<-ifelse(covid$CLASIFICACION_FINAL==3, 1, 0)
#### Main table ####
covid_pcr<-covid%>%filter(muestra==2)
covid_ag<-covid%>%filter(muestra==1)
covid_both<-covid%>%filter(muestra==3)

# Age
pcr0<-paste0(round(mean(covid_pcr$EDAD), 1),"±",round(sd(covid_pcr$EDAD),1))
ag0<-paste0(round(mean(covid_ag$EDAD), 1),"±",round(sd(covid_ag$EDAD),1))
both0<-paste0(round(mean(covid_both$EDAD), 1),"±",round(sd(covid_both$EDAD),1))

#Sexo
pcr1<-table(covid_pcr$SEXO)[2];pcr1.1<-(table(covid_pcr$SEXO)/nrow(covid_pcr))[2]
ag1<-table(covid_ag$SEXO)[2];ag1.1<-(table(covid_ag$SEXO)/nrow(covid_ag))[2]
both1<-table(covid_both$SEXO)[2];both1.1<-(table(covid_both$SEXO)/nrow(covid_both))[2]

# Confirmed COVID-19
pcrcov<-table(covid_pcr$muestra)[2];pcrcov.1<-(table(covid_pcr$muestra)/nrow(covid_pcr))[2]
agcov<-table(covid_ag$muestra)[2];agcov.1<-(table(covid_ag$muestra)/nrow(covid_ag))[2]
bothcov<-table(covid_both$muestra)[2];bothcov.1<-(table(covid_both$muestra)/nrow(covid_both))[2]

#Diabetes
pcr2<-table(covid_pcr$DIABETES)[2];pcr2.1<-(table(covid_pcr$DIABETES)/nrow(covid_pcr))[2]
ag2<-table(covid_ag$DIABETES)[2];ag2.1<-(table(covid_ag$DIABETES)/nrow(covid_ag))[2]
both2<-table(covid_both$DIABETES)[2];both2.1<-(table(covid_both$DIABETES)/nrow(covid_both))[2]

#EPOC
pcr3<-table(covid_pcr$EPOC)[2];pcr3.1<-(table(covid_pcr$EPOC)/nrow(covid_pcr))[2]
ag3<-table(covid_ag$EPOC)[2];ag3.1<-(table(covid_ag$EPOC)/nrow(covid_ag))[2]
both3<-table(covid_both$EPOC)[2];both3.1<-(table(covid_both$EPOC)/nrow(covid_both))[2]

#ASMA
pcr4<-table(covid_pcr$ASMA)[2];pcr4.1<-(table(covid_pcr$ASMA)/nrow(covid_pcr))[2]
ag4<-table(covid_ag$ASMA)[2];ag4.1<-(table(covid_ag$ASMA)/nrow(covid_ag))[2]
both4<-table(covid_both$ASMA)[2];both4.1<-(table(covid_both$ASMA)/nrow(covid_both))[2]

#INMUSUPR
pcr5<-table(covid_pcr$INMUSUPR)[2];pcr5.1<-(table(covid_pcr$INMUSUPR)/nrow(covid_pcr))[2]
ag5<-table(covid_ag$INMUSUPR)[2];ag5.1<-(table(covid_ag$INMUSUPR)/nrow(covid_ag))[2]
both5<-table(covid_both$INMUSUPR)[2];both5.1<-(table(covid_both$INMUSUPR)/nrow(covid_both))[2]

#HIPERTENSION
pcr6<-table(covid_pcr$HIPERTENSION)[2];pcr6.1<-(table(covid_pcr$HIPERTENSION)/nrow(covid_pcr))[2]
ag6<-table(covid_ag$HIPERTENSION)[2];ag6.1<-(table(covid_ag$HIPERTENSION)/nrow(covid_ag))[2]
both6<-table(covid_both$HIPERTENSION)[2];both6.1<-(table(covid_both$HIPERTENSION)/nrow(covid_both))[2]

#OTRA_COM
pcr7<-table(covid_pcr$OTRA_COM)[2];pcr7.1<-(table(covid_pcr$OTRA_COM)/nrow(covid_pcr))[2]
ag7<-table(covid_ag$OTRA_COM)[2];ag7.1<-(table(covid_ag$OTRA_COM)/nrow(covid_ag))[2]
both7<-table(covid_both$OTRA_COM)[2];both7.1<-(table(covid_both$OTRA_COM)/nrow(covid_both))[2]

#CARDIOVASCULAR
pcr8<-table(covid_pcr$CARDIOVASCULAR)[2];pcr8.1<-(table(covid_pcr$CARDIOVASCULAR)/nrow(covid_pcr))[2]
ag8<-table(covid_ag$CARDIOVASCULAR)[2];ag8.1<-(table(covid_ag$CARDIOVASCULAR)/nrow(covid_ag))[2]
both8<-table(covid_both$CARDIOVASCULAR)[2];both8.1<-(table(covid_both$CARDIOVASCULAR)/nrow(covid_both))[2]

#OBESIDAD
pcr9<-table(covid_pcr$OBESIDAD)[2];pcr9.1<-(table(covid_pcr$OBESIDAD)/nrow(covid_pcr))[2]
ag9<-table(covid_ag$OBESIDAD)[2];ag9.1<-(table(covid_ag$OBESIDAD)/nrow(covid_ag))[2]
both9<-table(covid_both$OBESIDAD)[2];both9.1<-(table(covid_both$OBESIDAD)/nrow(covid_both))[2]

#RENAL_CRONICA
pcr10<-table(covid_pcr$RENAL_CRONICA)[2];pcr10.1<-(table(covid_pcr$RENAL_CRONICA)/nrow(covid_pcr))[2]
ag10<-table(covid_ag$RENAL_CRONICA)[2];ag10.1<-(table(covid_ag$RENAL_CRONICA)/nrow(covid_ag))[2]
both10<-table(covid_both$RENAL_CRONICA)[2];both10.1<-(table(covid_both$RENAL_CRONICA)/nrow(covid_both))[2]

#TABAQUISMO
pcr11<-table(covid_pcr$TABAQUISMO)[2];pcr11.1<-(table(covid_pcr$TABAQUISMO)/nrow(covid_pcr))[2]
ag11<-table(covid_ag$TABAQUISMO)[2];ag11.1<-(table(covid_ag$TABAQUISMO)/nrow(covid_ag))[2]
both11<-table(covid_both$TABAQUISMO)[2];both11.1<-(table(covid_both$TABAQUISMO)/nrow(covid_both))[2]

#NEUMONIA
pcr12<-table(covid_pcr$NEUMONIA)[2];pcr12.1<-(table(covid_pcr$NEUMONIA)/nrow(covid_pcr))[2]
ag12<-table(covid_ag$NEUMONIA)[2];ag12.1<-(table(covid_ag$NEUMONIA)/nrow(covid_ag))[2]
both12<-table(covid_both$NEUMONIA)[2];both12.1<-(table(covid_both$NEUMONIA)/nrow(covid_both))[2]

#HOSPITALIZATION
pcr13<-table(covid_pcr$TIPO_PACIENTE)[2];pcr13.1<-(table(covid_pcr$TIPO_PACIENTE)/nrow(covid_pcr))[2]
ag13<-table(covid_ag$TIPO_PACIENTE)[2];ag13.1<-(table(covid_ag$TIPO_PACIENTE)/nrow(covid_ag))[2]
both13<-table(covid_both$TIPO_PACIENTE)[2];both13.1<-(table(covid_both$TIPO_PACIENTE)/nrow(covid_both))[2]

#ICU
pcr14<-table(covid_pcr$UCI)[2];pcr14.1<-(table(covid_pcr$UCI)/nrow(covid_pcr))[2]
ag14<-table(covid_ag$UCI)[2];ag14.1<-(table(covid_ag$UCI)/nrow(covid_ag))[2]
both14<-table(covid_both$UCI)[2];both14.1<-(table(covid_both$UCI)/nrow(covid_both))[2]

#ICU
pcr17<-table(covid_pcr$INTUBADO)[2];pcr17.1<-(table(covid_pcr$INTUBADO)/nrow(covid_pcr))[2]
ag17<-table(covid_ag$INTUBADO)[2];ag17.1<-(table(covid_ag$INTUBADO)/nrow(covid_ag))[2]
both17<-table(covid_both$INTUBADO)[2];both17.1<-(table(covid_both$INTUBADO)/nrow(covid_both))[2]

#DEATH
pcr15<-table(covid_pcr$Mortalidad)[2];pcr15.1<-(table(covid_pcr$Mortalidad)/nrow(covid_pcr))[2]
ag15<-table(covid_ag$Mortalidad)[2];ag15.1<-(table(covid_ag$Mortalidad)/nrow(covid_ag))[2]
both15<-table(covid_both$Mortalidad)[2];both15.1<-(table(covid_both$Mortalidad)/nrow(covid_both))[2]

#Days to assessment
pcr16<-paste0(round(median(covid_pcr$tiempo), 1)," (",round(quantile(covid_pcr$tiempo)[2],1),"-",round(quantile(covid_pcr$tiempo)[4],1),")")
ag16<-paste0(round(median(covid_ag$tiempo), 1)," (",round(quantile(covid_ag$tiempo)[2],1),"-",round(quantile(covid_ag$tiempo)[4],1),")")
both16<-paste0(round(median(covid_both$tiempo), 1)," (",round(quantile(covid_both$tiempo)[2],1),"-",round(quantile(covid_both$tiempo)[4],1),")")

pcr_1<-paste0(pcr1," ","(",round(pcr1.1*100,1),")");ag_1<-paste0(ag1," ","(",round(ag1.1*100,1),")");both_1<-paste0(both1," ","(",round(both1.1*100,1),")")
pcr_2<-paste0(pcr2," ","(",round(pcr2.1*100,1),")");ag_2<-paste0(ag2," ","(",round(ag2.1*100,1),")");both_2<-paste0(both2," ","(",round(both2.1*100,1),")")
pcr_3<-paste0(pcr3," ","(",round(pcr3.1*100,1),")");ag_3<-paste0(ag3," ","(",round(ag3.1*100,1),")");both_3<-paste0(both3," ","(",round(both3.1*100,1),")")
pcr_4<-paste0(pcr4," ","(",round(pcr4.1*100,1),")");ag_4<-paste0(ag4," ","(",round(ag4.1*100,1),")");both_4<-paste0(both4," ","(",round(both4.1*100,1),")")
pcr_5<-paste0(pcr5," ","(",round(pcr5.1*100,1),")");ag_5<-paste0(ag5," ","(",round(ag5.1*100,1),")");both_5<-paste0(both5," ","(",round(both5.1*100,1),")")
pcr_6<-paste0(pcr6," ","(",round(pcr6.1*100,1),")");ag_6<-paste0(ag6," ","(",round(ag6.1*100,1),")");both_6<-paste0(both6," ","(",round(both6.1*100,1),")")
pcr_7<-paste0(pcr7," ","(",round(pcr7.1*100,1),")");ag_7<-paste0(ag7," ","(",round(ag7.1*100,1),")");both_7<-paste0(both7," ","(",round(both7.1*100,1),")")
pcr_8<-paste0(pcr8," ","(",round(pcr8.1*100,1),")");ag_8<-paste0(ag8," ","(",round(ag8.1*100,1),")");both_8<-paste0(both8," ","(",round(both8.1*100,1),")")
pcr_9<-paste0(pcr9," ","(",round(pcr9.1*100,1),")");ag_9<-paste0(ag9," ","(",round(ag9.1*100,1),")");both_9<-paste0(both9," ","(",round(both9.1*100,1),")")
pcr_10<-paste0(pcr10," ","(",round(pcr10.1*100,1),")");ag_10<-paste0(ag10," ","(",round(ag10.1*100,1),")");both_10<-paste0(both10," ","(",round(both10.1*100,1),")")
pcr_11<-paste0(pcr11," ","(",round(pcr11.1*100,1),")");ag_11<-paste0(ag11," ","(",round(ag11.1*100,1),")");both_11<-paste0(both11," ","(",round(both11.1*100,1),")")
pcr_12<-paste0(pcr12," ","(",round(pcr12.1*100,1),")");ag_12<-paste0(ag12," ","(",round(ag12.1*100,1),")");both_12<-paste0(both12," ","(",round(both12.1*100,1),")")
pcr_13<-paste0(pcr13," ","(",round(pcr13.1*100,1),")");ag_13<-paste0(ag13," ","(",round(ag13.1*100,1),")");both_13<-paste0(both13," ","(",round(both13.1*100,1),")")
pcr_14<-paste0(pcr14," ","(",round(pcr14.1*100,1),")");ag_14<-paste0(ag14," ","(",round(ag14.1*100,1),")");both_14<-paste0(both14," ","(",round(both14.1*100,1),")")
pcr_15<-paste0(pcr15," ","(",round(pcr15.1*100,1),")");ag_15<-paste0(ag15," ","(",round(ag15.1*100,1),")");both_15<-paste0(both15," ","(",round(both15.1*100,1),")")
pcr_cov<-paste0(pcrcov," ","(",round(pcrcov.1*100,1),")");ag_cov<-paste0(agcov," ","(",round(agcov.1*100,1),")");both_cov<-paste0(bothcov," ","(",round(bothcov.1*100,1),")")
pcr_17<-paste0(pcr15," ","(",round(pcr17.1*100,1),")");ag_17<-paste0(ag17," ","(",round(ag17.1*100,1),")");both_17<-paste0(both17," ","(",round(both17.1*100,1),")")


covid_tests<-covid %>% filter(muestra!=0)
P_0<-summary(aov(covid_tests$EDAD~covid_tests$muestra))[[1]][["Pr(>F)"]][1]
P_1<-round(chisq.test(covid_tests$SEXO,covid_tests$muestra)$p.value,3);P_1
P_2<-round(chisq.test(covid_tests$DIABETES,covid_tests$muestra)$p.value,3);P_2
P_3<-round(chisq.test(covid_tests$EPOC,covid_tests$muestra)$p.value,3);P_3
P_4<-round(chisq.test(covid_tests$ASMA,covid_tests$muestra)$p.value,3);P_4
P_5<-round(chisq.test(covid_tests$INMUSUPR,covid_tests$muestra)$p.value,3);P_5
P_6<-round(chisq.test(covid_tests$HIPERTENSION,covid_tests$muestra)$p.value,3);P_6
P_7<-round(chisq.test(covid_tests$OTRA_COM,covid_tests$muestra)$p.value,3);P_7
P_8<-round(chisq.test(covid_tests$CARDIOVASCULAR,covid_tests$muestra)$p.value,3);P_8
P_9<-round(chisq.test(covid_tests$OBESIDAD,covid_tests$muestra)$p.value,3);P_9
P_10<-round(chisq.test(covid_tests$RENAL_CRONICA,covid_tests$muestra)$p.value,3);P_10
P_11<-round(chisq.test(covid_tests$TABAQUISMO,covid_tests$muestra)$p.value,3);P_11
P_12<-round(chisq.test(covid_tests$NEUMONIA,covid_tests$muestra)$p.value,3);P_12
P_13<-round(chisq.test(covid_tests$TIPO_PACIENTE,covid_tests$muestra)$p.value,3);P_13
P_14<-round(chisq.test(covid_tests$UCI,covid_tests$muestra)$p.value,3);P_14
P_15<-round(chisq.test(covid_tests$Mortalidad,covid_tests$muestra)$p.value,3);P_15
P_cov<-round(chisq.test(covid_tests$muestra,covid_tests$muestra)$p.value,3);P_cov
P_16<-kruskal.test(covid_tests$tiempo, covid_tests$muestra)$p.value;P_16
P_17<-round(chisq.test(covid_tests$Mortalidad,covid_tests$muestra)$p.value,3);P_1

GGG<-base::format.pval(c(P_0,P_1,P_cov,P_2,P_3,P_4,P_5,P_6,P_7,P_8,P_9,P_10,P_11,P_12,P_13,P_14,P_17,P_15, P_16), eps = .001, digits = 2)

AM_JOV_COVID<-data.frame("bold(Parameter)"=c("Age (years)","Male sex (%)","Confirmed SARS-CoV-2 (%)","Diabetes (%)","COPD (%)","Asthma (%)","Immunosuppression (%)","Hypertension (%)",
                                                 "Other (%)","CVD (%)","Obesity","CKD (%)","Smoking (%)","Pneumonia(%)",
                                                 "Hospitalization(%)","ICU admission (%)", "Intubation (%)","Death (%)", "Time to assessment* (days)"),
                         "RT-PCR test"=c(pcr0, pcr_1, pcr_cov,pcr_2, pcr_3, pcr_4, pcr_5, pcr_6, pcr_7, pcr_8, pcr_9, pcr_10, pcr_11, pcr_12, pcr_13, pcr_14, pcr_17,pcr_15, pcr16),
                         "Rapid Ag-T"=c(ag0,ag_1, ag_cov,ag_2, ag_3, ag_4, ag_5, ag_6, ag_7, ag_8, ag_9, ag_10, ag_11, ag_12, ag_13, ag_14, ag_17,ag_15, ag16),
                         "Both tests"=c(both0, both_1, both_cov,both_2, both_3, both_4, both_5, both_6, both_7, both_8, both_9, both_10, both_11, both_12, both_13, both_14, both_17,both_15, both16),
                         "p-value"=GGG)

colnames(AM_JOV_COVID)<-c("Parameters",paste0("RT-PCR test \nn=",nrow(covid_pcr)),
                          paste0("Rapig Ag test \nn=",nrow(covid_ag)),paste0("Both tests \nn=",nrow(covid_both)),"p-value")

AM_JOV_COVID<-flextable(AM_JOV_COVID,cwidth = 0.5*ncol(AM_JOV_COVID)) %>% autofit()


doc <- read_docx() %>%
  body_add_flextable(value = AM_JOV_COVID, split = TRUE) %>%
  body_end_section_landscape() %>% # a landscape section is ending here
  print( target = "table1.docx" )

#### Population analysis ####

tests<-covid %>% filter(FECHA_SINTOMAS>"2020-11-05") %>%
  group_by(ENTIDAD_UM,TOMA_MUESTRA_ANTIGENO) %>%
  summarise(pruebas=n(), pop=median(pop))%>% mutate(rate=(pruebas/pop)*100000) 

covid_tests<- covid %>% group_by(ENTIDAD_UM,) %>%
  summarise(tests=n(), pop=median(pop)) %>% 
  mutate(rate1=(tests/pop)*100000, ENTIDAD_UM=as.factor(ENTIDAD_UM))
covid_tests$state<-c("AGS", "BCN", "BCS", "CAM", "CHP", "CHH", "COA","COL", 
                     "CMX", "DUR", "GUA", "GRO", "HID", "JAL", "MEX", "MIC", 
                     "MOR", "NAY", "NLE", "OAX", "PUE", "QUE", "ROO", "SLP",
                     "SIN", "SON", "TAB", "TAM", "TLA", "VER", "YUC", "ZAC")
covid$cdmx<-factor(ifelse(covid$ENTIDAD_UM==9,1 ,0), labels = c("Rest of Mexico", "Mexico City"))

f1a<-covid %>% filter(FECHA_SINTOMAS>"2020-11-05")%>%
  group_by(ENTIDAD_UM) %>%
  summarise(pruebas=sum(TOMA_MUESTRA_ANTIGENO), pop=median(pop))%>% 
  mutate(rate=(pruebas/pop)*100000) %>%
  left_join(covid_tests%>%dplyr::select(ENTIDAD_UM, state))%>%
  mutate(state=as.factor(state))%>%
  mutate(state=forcats::fct_reorder(state, rate))%>%
  filter(rate>=1) %>%
  ggplot(aes(x=state, y=rate))+
  geom_col(width = 0.7)+
  coord_flip()+theme_bw()+labs(fill="Type of test")+
  ylab("Rapid antigen tests per 100,000 hab. (log scale)")+xlab("Mexican state")+
  scale_y_log10()

f1b<-covid %>% group_by(FECHA_SINTOMAS, cdmx) %>%
  summarise(pcr=sum(TOMA_MUESTRA_LAB), antigeno=sum(TOMA_MUESTRA_ANTIGENO))%>%
  mutate(ratio=(antigeno/(pcr+antigeno))*100) %>% filter(antigeno!=0) %>% group_by(cdmx)%>%
  mutate(ratio2=rollmean(ratio, 7, fill=NA))%>% ungroup()%>%
  ggplot(aes(x=FECHA_SINTOMAS, y=ratio2, col=cdmx))+geom_point()+geom_line()+scale_color_jama()+
  theme_classic()+ylab("Assessment with rapid antigen tests (%) 7-day rolling average")+
  xlab("Date at symptom onset")+labs(col="Location")+theme(legend.position = "top")

covid$covid1<-factor(covid$covid, labels=c("SARS-CoV-2 (-)", "SARS-CoV-2 (+)"))
covid$muestra1<-factor(covid$muestra, labels = c("Not tested", "Ag-T", "RT-PCR", "Both"))

f1c<-covid %>% filter(CLASIFICACION_FINAL==3)%>%
  mutate(test=factor(TOMA_MUESTRA_ANTIGENO, labels = c("RT-PCR", "Rapid test")))%>%
  filter(FECHA_SINTOMAS>"2020-11-04") %>%
  ggplot(aes(x=FECHA_SINTOMAS, fill=test))+
  geom_histogram(col="black", binwidth = 1)+
  ylab("New confirmed COVID-19 cases")+
  xlab("Symptom onset (date)")+
  theme_classic()+facet_wrap(~cdmx, scales = "free")+
  geom_vline(xintercept = 5, size = 1, colour = "#FF3721",linetype = "dashed")+
  labs(fill="Test")+theme(legend.position="top")+scale_fill_jama()

f1top<-ggarrange(f1a, f1b, labels = LETTERS[1:2])
fig1<-ggarrange(f1top,f1c, nrow = 2, ncol=1, labels = c("", "C"))

ggsave(fig1,filename = "Figure1.jpg", 
       width = 35, 
       height = 30,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)

hist<-covid %>%filter(tiempo<30) %>%
  ggplot(aes(x=tiempo,fill=covid1))+
  geom_histogram(col="black", binwidth = 1)+
  ylab("Tested cases")+
  xlab("Time from symptom onset to evaluation (days)")+
  theme_bw()+facet_wrap(~muestra1, scales = "free")+
  geom_vline(xintercept = 7, size = 1, colour = "#FF3721",linetype = "dashed")+
  labs(fill="Test")+theme(legend.position="top")+scale_fill_jama()


#### Diagnostic performance ####
covid0<-covid %>% filter(muestra==3)

covid1<-covid %>% filter(RESULTADO_LAB<3 & !is.na(RESULTADO_ANTIGENO))
nrow(covid0)
caret::confusionMatrix(factor(covid1$RESULTADO_ANTIGENO),factor(covid1$RESULTADO_LAB))
k1<-vcd::Kappa(table(covid1$RESULTADO_LAB, covid1$RESULTADO_ANTIGENO))
cbind(k1[1]$Unweighted[1], confint(k1))

pROC::roc(covid1$RESULTADO_LAB, covid1$RESULTADO_ANTIGENO, ci=T)

covid1$tiempo7<-ifelse(covid1$tiempo>7, 1, 0)

p1 <- optimal.cutpoints(X = "RESULTADO_ANTIGENO", status = "RESULTADO_LAB", tag.healthy = 0, 
                        methods = "Youden", data = covid1, pop.prev = NULL, 
                        control = control.cutpoints(), ci.fit = TRUE, conf.level = 0.95, trace = FALSE)

p2 <- optimal.cutpoints(X = "RESULTADO_ANTIGENO", status = "RESULTADO_LAB", tag.healthy = 0, 
                        methods = "Youden", data = covid1, pop.prev = NULL, categorical.cov = "cdmx",
                        control = control.cutpoints(), ci.fit = TRUE, conf.level = 0.95, trace = FALSE)
p3 <- optimal.cutpoints(X = "RESULTADO_ANTIGENO", status = "RESULTADO_LAB", tag.healthy = 0, 
                        methods = "Youden", data = covid1, pop.prev = NULL, categorical.cov = "TIPO_PACIENTE",
                        control = control.cutpoints(), ci.fit = TRUE, conf.level = 0.95, trace = FALSE)
p4 <- optimal.cutpoints(X = "RESULTADO_ANTIGENO", status = "RESULTADO_LAB", tag.healthy = 0, 
                        methods = "Youden", data = covid1, pop.prev = NULL, categorical.cov = "comorb_dic",
                        control = control.cutpoints(), ci.fit = TRUE, conf.level = 0.95, trace = FALSE)
p5 <- optimal.cutpoints(X = "RESULTADO_ANTIGENO", status = "RESULTADO_LAB", tag.healthy = 0, 
                        methods = "Youden", data = covid1, pop.prev = NULL, categorical.cov = "edad60",
                        control = control.cutpoints(), ci.fit = TRUE, conf.level = 0.95, trace = FALSE)
p6 <- optimal.cutpoints(X = "RESULTADO_ANTIGENO", status = "RESULTADO_LAB", tag.healthy = 0, 
                        methods = "Youden", data = covid1, pop.prev = NULL, categorical.cov = "tiempo7",
                        control = control.cutpoints(), ci.fit = TRUE, conf.level = 0.95, trace = FALSE)

table1<-data.frame("Parameter"=c("Overall", "Mexico City","Rest of Mexico","Outpatient", "Inpatient", "No comorbidities", "≥1 comorbidity", "<60 years", "≥60 years", "<7d from onset", "≥7d from onset"),
                   "False positives/False negatives"=c(paste0(p1[1]$Youden$Global$optimal.cutoff$FP[1],"/",p1[1]$Youden$Global$optimal.cutoff$FN[1]),
                                                       paste0(p2[1]$Youden$`Mexico City`$optimal.cutoff$FP[1],"/",p2[1]$Youden$`Mexico City`$optimal.cutoff$FN[1]),
                                                       paste0(p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$FP[1],"/",p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$FN[1]),
                                                       paste0(round(p3[1]$Youden$`0`$optimal.cutoff$FP[1],0),"/",p3[1]$Youden$`0`$optimal.cutoff$FN[1]),
                                                       paste0(p3[1]$Youden$`1`$optimal.cutoff$FP[1],"/",p3[1]$Youden$`1`$optimal.cutoff$FN[1]),
                                                       paste0(p4[1]$Youden$`0`$optimal.cutoff$FP[1],"/",p4[1]$Youden$`0`$optimal.cutoff$FN[1]),
                                                       paste0(p4[1]$Youden$`1`$optimal.cutoff$FP[1],"/",p4[1]$Youden$`1`$optimal.cutoff$FN[1]),
                                                       paste0(p5[1]$Youden$`0`$optimal.cutoff$FP[1],"/",p5[1]$Youden$`0`$optimal.cutoff$FN[1]),
                                                       paste0(p5[1]$Youden$`1`$optimal.cutoff$FP[1],"/",p5[1]$Youden$`1`$optimal.cutoff$FN[1]),
                                                       paste0(p6[1]$Youden$`0`$optimal.cutoff$FP[1],"/",p6[1]$Youden$`0`$optimal.cutoff$FN[1]),
                                                       paste0(p6[1]$Youden$`1`$optimal.cutoff$FP[1],"/",p6[1]$Youden$`1`$optimal.cutoff$FN[1])),
                   "AUROC (95%CI)"=c(paste0(round(p1[1]$Youden$Global$measures.acc$AUC[1],3), " (",round(p1[1]$Youden$Global$measures.acc$AUC[2],3), "-",round(p1[1]$Youden$Global$measures.acc$AUC[3],3),")"),
                                     paste0(round(p2[1]$Youden$`Mexico City`$measures.acc$AUC[1],3), " (",round(p2[1]$Youden$`Mexico City`$measures.acc$AUC[2],3), "-",round(p2[1]$Youden$`Mexico City`$measures.acc$AUC[3],3),")"),
                                     paste0(round(p2[1]$Youden$`Rest of Mexico`$measures.acc$AUC[1],3), " (",round(p2[1]$Youden$`Rest of Mexico`$measures.acc$AUC[2],3), "-",round(p2[1]$Youden$`Rest of Mexico`$measures.acc$AUC[3],3),")"),
                                     paste0(round(p3[1]$Youden$`0`$measures.acc$AUC[1],3), " (",round(p3[1]$Youden$`0`$measures.acc$AUC[2],3), "-",round(p3[1]$Youden$`0`$measures.acc$AUC[3],3),")"),
                                     paste0(round(p3[1]$Youden$`1`$measures.acc$AUC[1],3), " (",round(p3[1]$Youden$`1`$measures.acc$AUC[2],3), "-",round(p3[1]$Youden$`1`$measures.acc$AUC[3],3),")"),
                                     paste0(round(p4[1]$Youden$`0`$measures.acc$AUC[1],3), " (",round(p4[1]$Youden$`0`$measures.acc$AUC[2],3), "-",round(p4[1]$Youden$`0`$measures.acc$AUC[3],3),")"),
                                     paste0(round(p4[1]$Youden$`1`$measures.acc$AUC[1],3), " (",round(p4[1]$Youden$`1`$measures.acc$AUC[2],3), "-",round(p4[1]$Youden$`1`$measures.acc$AUC[3],3),")"),
                                     paste0(round(p5[1]$Youden$`0`$measures.acc$AUC[1],3), " (",round(p5[1]$Youden$`0`$measures.acc$AUC[2],3), "-",round(p5[1]$Youden$`0`$measures.acc$AUC[3],3),")"),
                                     paste0(round(p5[1]$Youden$`1`$measures.acc$AUC[1],3), " (",round(p5[1]$Youden$`1`$measures.acc$AUC[2],3), "-",round(p5[1]$Youden$`1`$measures.acc$AUC[3],3),")"),
                                     paste0(round(p6[1]$Youden$`0`$measures.acc$AUC[1],3), " (",round(p6[1]$Youden$`0`$measures.acc$AUC[2],3), "-",round(p6[1]$Youden$`0`$measures.acc$AUC[3],3),")"),
                                     paste0(round(p6[1]$Youden$`1`$measures.acc$AUC[1],3), " (",round(p6[1]$Youden$`1`$measures.acc$AUC[2],3), "-",round(p6[1]$Youden$`1`$measures.acc$AUC[3],3),")")),
                   "Sensitivity (%, 95%CI)"=c(paste0(round(p1[1]$Youden$Global$optimal.cutoff$Se[1],3)*100, " (",round(p1[1]$Youden$Global$optimal.cutoff$Se[2],3)*100, "-",round(p1[1]$Youden$Global$optimal.cutoff$Se[3],3)*100,")"),
                                              paste0(round(p2[1]$Youden$`Mexico City`$optimal.cutoff$Se[1],3)*100, " (",round(p2[1]$Youden$`Mexico City`$optimal.cutoff$Se[2],3)*100, "-",round(p2[1]$Youden$`Mexico City`$optimal.cutoff$Se[3],3)*100,")"),
                                              paste0(round(p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$Se[1],3)*100, " (",round(p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$Se[2],3)*100, "-",round(p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$Se[3],3)*100,")"),
                                              paste0(round(p3[1]$Youden$`0`$optimal.cutoff$Se[1],3)*100, " (",round(p3[1]$Youden$`0`$optimal.cutoff$Se[2],3)*100, "-",round(p3[1]$Youden$`0`$optimal.cutoff$Se[3],3)*100,")"),
                                              paste0(round(p3[1]$Youden$`1`$optimal.cutoff$Se[1],3)*100, " (",round(p3[1]$Youden$`1`$optimal.cutoff$Se[2],3)*100, "-",round(p3[1]$Youden$`1`$optimal.cutoff$Se[3],3)*100,")"),
                                              paste0(round(p4[1]$Youden$`0`$optimal.cutoff$Se[1],3)*100, " (",round(p4[1]$Youden$`0`$optimal.cutoff$Se[2],3)*100, "-",round(p4[1]$Youden$`0`$optimal.cutoff$Se[3],3)*100,")"),
                                              paste0(round(p4[1]$Youden$`1`$optimal.cutoff$Se[1],3)*100, " (",round(p4[1]$Youden$`1`$optimal.cutoff$Se[2],3)*100, "-",round(p4[1]$Youden$`1`$optimal.cutoff$Se[3],3)*100,")"),
                                              paste0(round(p5[1]$Youden$`0`$optimal.cutoff$Se[1],3)*100, " (",round(p5[1]$Youden$`0`$optimal.cutoff$Se[2],3)*100, "-",round(p5[1]$Youden$`0`$optimal.cutoff$Se[3],3)*100,")"),
                                              paste0(round(p5[1]$Youden$`1`$optimal.cutoff$Se[1],3)*100, " (",round(p5[1]$Youden$`1`$optimal.cutoff$Se[2],3)*100, "-",round(p5[1]$Youden$`1`$optimal.cutoff$Se[3],3)*100,")"),
                                              paste0(round(p6[1]$Youden$`0`$optimal.cutoff$Se[1],3)*100, " (",round(p6[1]$Youden$`0`$optimal.cutoff$Se[2],3)*100, "-",round(p6[1]$Youden$`0`$optimal.cutoff$Se[3],3)*100,")"),
                                              paste0(round(p6[1]$Youden$`1`$optimal.cutoff$Se[1],3)*100, " (",round(p6[1]$Youden$`1`$optimal.cutoff$Se[2],3)*100, "-",round(p6[1]$Youden$`1`$optimal.cutoff$Se[3],3)*100,")")),
                   "Specificity (%, 95%CI)"=c(paste0(round(p1[1]$Youden$Global$optimal.cutoff$Sp[1],3)*100, " (",round(p1[1]$Youden$Global$optimal.cutoff$Sp[2],3)*100, "-",round(p1[1]$Youden$Global$optimal.cutoff$Sp[3],3)*100,")"),
                                              paste0(round(p2[1]$Youden$`Mexico City`$optimal.cutoff$Sp[1],3)*100, " (",round(p2[1]$Youden$`Mexico City`$optimal.cutoff$Sp[2],3)*100, "-",round(p2[1]$Youden$`Mexico City`$optimal.cutoff$Sp[3],3)*100,")"),
                                              paste0(round(p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$Sp[1],3)*100, " (",round(p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$Sp[2],3)*100, "-",round(p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$Sp[3],3)*100,")"),
                                              paste0(round(p3[1]$Youden$`0`$optimal.cutoff$Sp[1],3)*100, " (",round(p3[1]$Youden$`0`$optimal.cutoff$Sp[2],3)*100, "-",round(p3[1]$Youden$`0`$optimal.cutoff$Sp[3],3)*100,")"),
                                              paste0(round(p3[1]$Youden$`1`$optimal.cutoff$Sp[1],3)*100, " (",round(p3[1]$Youden$`1`$optimal.cutoff$Sp[2],3)*100, "-",round(p3[1]$Youden$`1`$optimal.cutoff$Sp[3],3)*100,")"),
                                              paste0(round(p4[1]$Youden$`0`$optimal.cutoff$Sp[1],3)*100, " (",round(p4[1]$Youden$`0`$optimal.cutoff$Sp[2],3)*100, "-",round(p4[1]$Youden$`0`$optimal.cutoff$Sp[3],3)*100,")"),
                                              paste0(round(p4[1]$Youden$`1`$optimal.cutoff$Sp[1],3)*100, " (",round(p4[1]$Youden$`1`$optimal.cutoff$Sp[2],3)*100, "-",round(p4[1]$Youden$`1`$optimal.cutoff$Sp[3],3)*100,")"),
                                              paste0(round(p5[1]$Youden$`0`$optimal.cutoff$Sp[1],3)*100, " (",round(p5[1]$Youden$`0`$optimal.cutoff$Sp[2],3)*100, "-",round(p5[1]$Youden$`0`$optimal.cutoff$Sp[3],3)*100,")"),
                                              paste0(round(p5[1]$Youden$`1`$optimal.cutoff$Sp[1],3)*100, " (",round(p5[1]$Youden$`1`$optimal.cutoff$Sp[2],3)*100, "-",round(p5[1]$Youden$`1`$optimal.cutoff$Sp[3],3)*100,")"),
                                              paste0(round(p6[1]$Youden$`0`$optimal.cutoff$Sp[1],3)*100, " (",round(p6[1]$Youden$`0`$optimal.cutoff$Sp[2],3)*100, "-",round(p6[1]$Youden$`0`$optimal.cutoff$Sp[3],3)*100,")"),
                                              paste0(round(p6[1]$Youden$`1`$optimal.cutoff$Sp[1],3)*100, " (",round(p6[1]$Youden$`1`$optimal.cutoff$Sp[2],3)*100, "-",round(p6[1]$Youden$`1`$optimal.cutoff$Sp[3],3)*100,")")),
                   "PPV (%, 95%CI)"=c(paste0(round(p1[1]$Youden$Global$optimal.cutoff$PPV[1],3)*100, " (",round(p1[1]$Youden$Global$optimal.cutoff$PPV[2],3)*100, "-",round(p1[1]$Youden$Global$optimal.cutoff$PPV[3],3)*100,")"),
                                     paste0(round(p2[1]$Youden$`Mexico City`$optimal.cutoff$PPV[1],3)*100, " (",round(p2[1]$Youden$`Mexico City`$optimal.cutoff$PPV[2],3)*100, "-",round(p2[1]$Youden$`Mexico City`$optimal.cutoff$PPV[3],3)*100,")"),
                                     paste0(round(p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$PPV[1],3)*100, " (",round(p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$PPV[2],3)*100, "-",round(p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$PPV[3],3)*100,")"),
                                     paste0(round(p3[1]$Youden$`0`$optimal.cutoff$PPV[1],3)*100, " (",round(p3[1]$Youden$`0`$optimal.cutoff$PPV[2],3)*100, "-",round(p3[1]$Youden$`0`$optimal.cutoff$PPV[3],3)*100,")"),
                                     paste0(round(p3[1]$Youden$`1`$optimal.cutoff$PPV[1],3)*100, " (",round(p3[1]$Youden$`1`$optimal.cutoff$PPV[2],3)*100, "-",round(p3[1]$Youden$`1`$optimal.cutoff$PPV[3],3)*100,")"),
                                     paste0(round(p4[1]$Youden$`0`$optimal.cutoff$PPV[1],3)*100, " (",round(p4[1]$Youden$`0`$optimal.cutoff$PPV[2],3)*100, "-",round(p4[1]$Youden$`0`$optimal.cutoff$PPV[3],3)*100,")"),
                                     paste0(round(p4[1]$Youden$`1`$optimal.cutoff$PPV[1],3)*100, " (",round(p4[1]$Youden$`1`$optimal.cutoff$PPV[2],3)*100, "-",round(p4[1]$Youden$`1`$optimal.cutoff$PPV[3],3)*100,")"),
                                     paste0(round(p5[1]$Youden$`0`$optimal.cutoff$PPV[1],3)*100, " (",round(p5[1]$Youden$`0`$optimal.cutoff$PPV[2],3)*100, "-",round(p5[1]$Youden$`0`$optimal.cutoff$PPV[3],3)*100,")"),
                                     paste0(round(p5[1]$Youden$`1`$optimal.cutoff$PPV[1],3)*100, " (",round(p5[1]$Youden$`1`$optimal.cutoff$PPV[2],3)*100, "-",round(p5[1]$Youden$`1`$optimal.cutoff$PPV[3],3)*100,")"),
                                     paste0(round(p6[1]$Youden$`0`$optimal.cutoff$PPV[1],3)*100, " (",round(p6[1]$Youden$`0`$optimal.cutoff$PPV[2],3)*100, "-",round(p6[1]$Youden$`0`$optimal.cutoff$PPV[3],3)*100,")"),
                                     paste0(round(p6[1]$Youden$`1`$optimal.cutoff$PPV[1],3)*100, " (",round(p6[1]$Youden$`1`$optimal.cutoff$PPV[2],3)*100, "-",round(p6[1]$Youden$`1`$optimal.cutoff$PPV[3],3)*100,")")),
                   "NPV (%, 95%CI)"=c(paste0(round(p1[1]$Youden$Global$optimal.cutoff$NPV[1],3)*100, " (",round(p1[1]$Youden$Global$optimal.cutoff$NPV[2],3)*100, "-",round(p1[1]$Youden$Global$optimal.cutoff$NPV[3],3)*100,")"),
                                      paste0(round(p2[1]$Youden$`Mexico City`$optimal.cutoff$NPV[1],3)*100, " (",round(p2[1]$Youden$`Mexico City`$optimal.cutoff$NPV[2],3)*100, "-",round(p2[1]$Youden$`Mexico City`$optimal.cutoff$NPV[3],3)*100,")"),
                                      paste0(round(p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$NPV[1],3)*100, " (",round(p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$NPV[2],3)*100, "-",round(p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$NPV[3],3)*100,")"),
                                      paste0(round(p3[1]$Youden$`0`$optimal.cutoff$NPV[1],3)*100, " (",round(p3[1]$Youden$`0`$optimal.cutoff$NPV[2],3)*100, "-",round(p3[1]$Youden$`0`$optimal.cutoff$NPV[3],3)*100,")"),
                                      paste0(round(p3[1]$Youden$`1`$optimal.cutoff$NPV[1],3)*100, " (",round(p3[1]$Youden$`1`$optimal.cutoff$NPV[2],3)*100, "-",round(p3[1]$Youden$`1`$optimal.cutoff$NPV[3],3)*100,")"),
                                      paste0(round(p4[1]$Youden$`0`$optimal.cutoff$NPV[1],3)*100, " (",round(p4[1]$Youden$`0`$optimal.cutoff$NPV[2],3)*100, "-",round(p4[1]$Youden$`0`$optimal.cutoff$NPV[3],3)*100,")"),
                                      paste0(round(p4[1]$Youden$`1`$optimal.cutoff$NPV[1],3)*100, " (",round(p4[1]$Youden$`1`$optimal.cutoff$NPV[2],3)*100, "-",round(p4[1]$Youden$`1`$optimal.cutoff$NPV[3],3)*100,")"),
                                      paste0(round(p5[1]$Youden$`0`$optimal.cutoff$NPV[1],3)*100, " (",round(p5[1]$Youden$`0`$optimal.cutoff$NPV[2],3)*100, "-",round(p5[1]$Youden$`0`$optimal.cutoff$NPV[3],3)*100,")"),
                                      paste0(round(p5[1]$Youden$`1`$optimal.cutoff$NPV[1],3)*100, " (",round(p5[1]$Youden$`1`$optimal.cutoff$NPV[2],3)*100, "-",round(p5[1]$Youden$`1`$optimal.cutoff$NPV[3],3)*100,")"),
                                      paste0(round(p6[1]$Youden$`0`$optimal.cutoff$NPV[1],3)*100, " (",round(p6[1]$Youden$`0`$optimal.cutoff$NPV[2],3)*100, "-",round(p6[1]$Youden$`0`$optimal.cutoff$NPV[3],3)*100,")"),
                                      paste0(round(p6[1]$Youden$`1`$optimal.cutoff$NPV[1],3)*100, " (",round(p6[1]$Youden$`1`$optimal.cutoff$NPV[2],3)*100, "-",round(p6[1]$Youden$`1`$optimal.cutoff$NPV[3],3)*100,")")),
                   "LR (+) (95%CI)"=c(paste0(round(p1[1]$Youden$Global$optimal.cutoff$DLR.Positive[1],1), " (",round(p1[1]$Youden$Global$optimal.cutoff$DLR.Positive[2],1), "-",round(p1[1]$Youden$Global$optimal.cutoff$DLR.Positive[3],1),")"),
                                      paste0(round(p2[1]$Youden$`Mexico City`$optimal.cutoff$DLR.Positive[1],1), " (",round(p2[1]$Youden$`Mexico City`$optimal.cutoff$DLR.Positive[2],1), "-",round(p2[1]$Youden$`Mexico City`$optimal.cutoff$DLR.Positive[3],1),")"),
                                      paste0(round(p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$DLR.Positive[1],1), " (",round(p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$DLR.Positive[2],1), "-",round(p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$DLR.Positive[3],1),")"),
                                      paste0(round(p3[1]$Youden$`0`$optimal.cutoff$DLR.Positive[1],1), " (",round(p3[1]$Youden$`0`$optimal.cutoff$DLR.Positive[2],1), "-",round(p3[1]$Youden$`0`$optimal.cutoff$DLR.Positive[3],1),")"),
                                      paste0(round(p3[1]$Youden$`1`$optimal.cutoff$DLR.Positive[1],1), " (",round(p3[1]$Youden$`1`$optimal.cutoff$DLR.Positive[2],1), "-",round(p3[1]$Youden$`1`$optimal.cutoff$DLR.Positive[3],1),")"),
                                      paste0(round(p4[1]$Youden$`0`$optimal.cutoff$DLR.Positive[1],1), " (",round(p4[1]$Youden$`0`$optimal.cutoff$DLR.Positive[2],1), "-",round(p4[1]$Youden$`0`$optimal.cutoff$DLR.Positive[3],1),")"),
                                      paste0(round(p4[1]$Youden$`1`$optimal.cutoff$DLR.Positive[1],1), " (",round(p4[1]$Youden$`1`$optimal.cutoff$DLR.Positive[2],1), "-",round(p4[1]$Youden$`1`$optimal.cutoff$DLR.Positive[3],1),")"),
                                      paste0(round(p5[1]$Youden$`0`$optimal.cutoff$DLR.Positive[1],1), " (",round(p5[1]$Youden$`0`$optimal.cutoff$DLR.Positive[2],1), "-",round(p5[1]$Youden$`0`$optimal.cutoff$DLR.Positive[3],1),")"),
                                      paste0(round(p5[1]$Youden$`1`$optimal.cutoff$DLR.Positive[1],1), " (",round(p5[1]$Youden$`1`$optimal.cutoff$DLR.Positive[2],1), "-",round(p5[1]$Youden$`1`$optimal.cutoff$DLR.Positive[3],1),")"),
                                      paste0(round(p6[1]$Youden$`0`$optimal.cutoff$DLR.Positive[1],1), " (",round(p6[1]$Youden$`0`$optimal.cutoff$DLR.Positive[2],1), "-",round(p6[1]$Youden$`0`$optimal.cutoff$DLR.Positive[3],1),")"),
                                      paste0(round(p6[1]$Youden$`1`$optimal.cutoff$DLR.Positive[1],1), " (",round(p6[1]$Youden$`1`$optimal.cutoff$DLR.Positive[2],1), "-",round(p6[1]$Youden$`1`$optimal.cutoff$DLR.Positive[3],1),")")),
                   "LR (-) (95%CI)"=c(paste0(round(p1[1]$Youden$Global$optimal.cutoff$DLR.Negative[1],2), " (",round(p1[1]$Youden$Global$optimal.cutoff$DLR.Negative[2],2), "-",round(p1[1]$Youden$Global$optimal.cutoff$DLR.Negative[3],2),")"),
                                      paste0(round(p2[1]$Youden$`Mexico City`$optimal.cutoff$DLR.Negative[1],2), " (",round(p2[1]$Youden$`Mexico City`$optimal.cutoff$DLR.Negative[2],2), "-",round(p2[1]$Youden$`Mexico City`$optimal.cutoff$DLR.Negative[3],2),")"),
                                      paste0(round(p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$DLR.Negative[1],2), " (",round(p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$DLR.Negative[2],2), "-",round(p2[1]$Youden$`Rest of Mexico`$optimal.cutoff$DLR.Negative[3],2),")"),
                                      paste0(round(p3[1]$Youden$`0`$optimal.cutoff$DLR.Negative[1],2), " (",round(p3[1]$Youden$`0`$optimal.cutoff$DLR.Negative[2],2), "-",round(p3[1]$Youden$`0`$optimal.cutoff$DLR.Negative[3],2),")"),
                                      paste0(round(p3[1]$Youden$`1`$optimal.cutoff$DLR.Negative[1],2), " (",round(p3[1]$Youden$`1`$optimal.cutoff$DLR.Negative[2],2), "-",round(p3[1]$Youden$`1`$optimal.cutoff$DLR.Negative[3],2),")"),
                                      paste0(round(p4[1]$Youden$`0`$optimal.cutoff$DLR.Negative[1],2), " (",round(p4[1]$Youden$`0`$optimal.cutoff$DLR.Negative[2],2), "-",round(p4[1]$Youden$`0`$optimal.cutoff$DLR.Negative[3],2),")"),
                                      paste0(round(p4[1]$Youden$`1`$optimal.cutoff$DLR.Negative[1],2), " (",round(p4[1]$Youden$`1`$optimal.cutoff$DLR.Negative[2],2), "-",round(p4[1]$Youden$`1`$optimal.cutoff$DLR.Negative[3],2),")"),
                                      paste0(round(p5[1]$Youden$`0`$optimal.cutoff$DLR.Negative[1],2), " (",round(p5[1]$Youden$`0`$optimal.cutoff$DLR.Negative[2],2), "-",round(p5[1]$Youden$`0`$optimal.cutoff$DLR.Negative[3],2),")"),
                                      paste0(round(p5[1]$Youden$`1`$optimal.cutoff$DLR.Negative[1],2), " (",round(p5[1]$Youden$`1`$optimal.cutoff$DLR.Negative[2],2), "-",round(p5[1]$Youden$`1`$optimal.cutoff$DLR.Negative[3],2),")"),
                                      paste0(round(p6[1]$Youden$`0`$optimal.cutoff$DLR.Negative[1],2), " (",round(p6[1]$Youden$`0`$optimal.cutoff$DLR.Negative[2],2), "-",round(p6[1]$Youden$`0`$optimal.cutoff$DLR.Negative[3],2),")"),
                                      paste0(round(p6[1]$Youden$`1`$optimal.cutoff$DLR.Negative[1],2), " (",round(p6[1]$Youden$`1`$optimal.cutoff$DLR.Negative[2],2), "-",round(p6[1]$Youden$`1`$optimal.cutoff$DLR.Negative[3],2),")")))

table1<-`names<-`(table1,c("Parameter","FP/FN","AUROC (95%CI)","Sensitivity (%, 95%CI)","Specificity (%, 95%CI)","PPV (%, 95%CI)","NPV (%, 95%CI)", "LR (+) (95%CI)", "LR (-) (95%CI)"))
table1<-flextable::align(flextable(table1),align = "center",part = "all") 

doc <- read_docx() %>%
  body_add_flextable(value = table1, split = TRUE) %>%
  body_end_section_landscape() %>% # a landscape section is ending here
  print( target = "table2.docx" )

covid1$antigeno<-factor(covid1$RESULTADO_ANTIGENO, labels = c("Negativo", "Positivo"))
covid1$pcr<-factor(covid1$RESULTADO_LAB, labels = c("Negativo", "Positivo"))

#### Time-dependent ROC curves ####
covid1<-covid %>% filter(RESULTADO_LAB<3 & !is.na(RESULTADO_ANTIGENO))

roc1<-timeROC(T=covid1$tiempo,
              delta=covid1$RESULTADO_LAB,marker=covid1$RESULTADO_ANTIGENO,
              cause=1,weighting="cox",
              other_markers=as.matrix(covid1$EDAD, covid1$SEXO),
              times=c(1,3,7,10, 15));roc1

roc2<-timeROC(T=covid1$tiempo,
              delta=covid1$TIPO_PACIENTE,marker=covid1$RESULTADO_LAB,
              cause=1,weighting="cox",
              other_markers=as.matrix(covid1$EDAD, covid1$SEXO),
              times=c(1,3,7,10, 15));roc2

roc3<-timeROC(T=covid1$tiempo,
              delta=covid1$TIPO_PACIENTE,marker=covid1$RESULTADO_ANTIGENO,
              cause=1,weighting="cox",
              other_markers=as.matrix(covid1$EDAD, covid1$SEXO),
              times=c(1,3,7,10, 15));roc3

roc4<-timeROC(T=covid1$tiempo,
              delta=covid1$INTUBADO,marker=covid1$RESULTADO_LAB,
              cause=1,weighting="cox",
              other_markers=as.matrix(covid1$EDAD, covid1$SEXO),
              times=c(1,3,7,10, 15));roc4

roc5<-timeROC(T=covid1$tiempo,
              delta=covid1$INTUBADO,marker=covid1$RESULTADO_ANTIGENO,
              cause=1,weighting="cox",
              other_markers=as.matrix(covid1$EDAD, covid1$SEXO),
              times=c(1,3,7,10, 15));roc5

roc6<-timeROC(T=covid1$tiempo,
              delta=covid1$Mortalidad,marker=covid1$RESULTADO_LAB,
              cause=1,weighting="cox",
              other_markers=as.matrix(covid1$EDAD, covid1$SEXO),
              times=c(1,3,7,10, 15));roc6

roc7<-timeROC(T=covid1$tiempo,
              delta=covid1$Mortalidad,marker=covid1$RESULTADO_ANTIGENO,
              cause=1,weighting="cox",
              other_markers=as.matrix(covid1$EDAD, covid1$SEXO),
              times=c(1,3,7,10, 15));roc7


table2<- data.frame(cols=c("1 day", "3 days", "7 days", "10 days", "15 days"),
                    pcr_vs_ag=c(paste0(round(roc1$AUC[1],3)),paste0(round(roc1$AUC[2],3)),
                                paste0(round(roc1$AUC[3],3)),paste0(round(roc1$AUC[4],3)),paste0(round(roc1$AUC[5],3))),
                    pcr_vs_hosp=c(paste0(round(roc2$AUC[1],3)),paste0(round(roc2$AUC[2],3)),
                                paste0(round(roc2$AUC[3],3)),paste0(round(roc2$AUC[4],3)),paste0(round(roc2$AUC[5],3))),
                    ag_vs_hosp=c(paste0(round(roc3$AUC[1],3)),paste0(round(roc3$AUC[2],3)),
                                paste0(round(roc3$AUC[3],3)),paste0(round(roc1$AUC[4],3)),paste0(round(roc1$AUC[5],3))),
                    pcr_vs_int=c(paste0(round(roc4$AUC[1],3)),paste0(round(roc4$AUC[2],3)),
                                paste0(round(roc4$AUC[3],3)),paste0(round(roc4$AUC[4],3)),paste0(round(roc4$AUC[5],3))),
                    ag_vs_int=c(paste0(round(roc5$AUC[1],3)),paste0(round(roc5$AUC[2],3)),
                                paste0(round(roc5$AUC[3],3)),paste0(round(roc5$AUC[4],3)),paste0(round(roc5$AUC[5],3))),
                    pcr_vs_mort=c(paste0(round(roc6$AUC[1],3)),paste0(round(roc6$AUC[2],3)),
                                paste0(round(roc6$AUC[3],3)),paste0(round(roc6$AUC[4],3)),paste0(round(roc6$AUC[5],3))),
                    ag_vs_mort=c(paste0(round(roc7$AUC[1],3)),paste0(round(roc7$AUC[2],3)),
                                paste0(round(roc7$AUC[3],3)),paste0(round(roc7$AUC[4],3)),paste0(round(roc7$AUC[5],3))))

table2<-`names<-`(table2,c("Time AUROC","PCR vs. Ag-T","RT-PCR \nHospitalization","Ag-T \nHospitalization",
                           "RT-PCR \nintubation","Ag-T \nIntubation","RT-PCR \nMortality", "Ag-T \nMortality"))
table2<-flextable::align(flextable(table2),align = "center",part = "all") 

doc <- read_docx() %>%
  body_add_flextable(value = table2, split = TRUE) %>%
  body_end_section_landscape() %>% # a landscape section is ending here
  print( target = "table3.docx" )


#### Discordant test models ####
covid1<-covid %>% filter(RESULTADO_LAB<3 & !is.na(RESULTADO_ANTIGENO))

covid1$false_status<-ifelse(covid1$RESULTADO_LAB==1, 1, 0)
covid1$false_status[covid1$RESULTADO_ANTIGENO==0 & covid1$RESULTADO_LAB==1]<-2
covid1$false_status[covid1$RESULTADO_ANTIGENO==1 & covid1$RESULTADO_LAB==0]<-3
covid1$false_status<-factor(covid1$false_status, labels = c("True negative", "True positive", "False negative", "False positive"))
covid1$false_positive<-ifelse(covid1$false_status=="False positive", 1, 0)
covid1$false_negative<-ifelse(covid1$false_status=="False negative", 1, 0)
covid1$discordant<-ifelse((covid1$false_status=="False negative" |covid1$false_status=="False positive"), 1, 0)

covid1.1<-covid1%>%filter(false_status=="True positive" |  false_status=="False negative")
covid1.2<-covid1%>%filter(false_status=="True negative" |  false_status=="False positive")

#False Negative Predictors model

m1<-glmer(false_negative~scale(EDAD)+(tiempo>=7)+SEXO+DIABETES+HIPERTENSION+EPOC+INMUSUPR+CARDIOVASCULAR+
            TABAQUISMO+OBESIDAD+ASMA+RENAL_CRONICA+(1|id), data=covid1.1, family="binomial")
summ(m1, exp=T, ci=T)

Fig1<-ggstatsplot::ggcoefstats(
  x = m1,
  exp=T,
  conf.int = TRUE,
  vline = F,
  exclude.intercept = T,
  vline.args = list(size = 1, linetype = "dashed"),
  package = "RColorBrewer",
  palette = "Dark2",
  ggtheme = ggthemes::theme_hc(),
  title = "False Negative Predictors for Rapid Antigen Test in Mexico",
  xlab = "Random Effect Logistic Regression Model; OR (95% CI)",
  ylab = "",
  mean.label.size = 3)+
  scale_x_log10(limits = c(-1, 10))+
  geom_vline(aes(xintercept = 1), size = 0.5, linetype = "dashed")+
  ggplot2::scale_y_discrete(labels = c("Age", "≥7 days Since Symptoms Onset","Male Sex", "Diabetes","Hypertension","COPD","Inmunosupression","CVD","Smoking Status",
                                       "Obesity","Asthma","CKD"))
Fig1

#False Positive Predictors model

m2<-glmer(false_positive~scale(EDAD)+(tiempo>=7)+SEXO+DIABETES+HIPERTENSION+EPOC+INMUSUPR+CARDIOVASCULAR+
            TABAQUISMO+OBESIDAD+ASMA+RENAL_CRONICA+(1|id), data=covid1.2, family="binomial")
summ(m2, exp=T, ci=T)

Fig2<-ggstatsplot::ggcoefstats(
  x = m2,
  exp=T,
  conf.int = TRUE,
  vline = F,
  exclude.intercept = T,
  vline.args = list(size = 1, linetype = "dashed"),
  package = "RColorBrewer",
  palette = "Dark2",
  ggtheme = ggthemes::theme_hc(),
  title = "False Positive Predictors for Rapid Antigen Test in Mexico",
  xlab = "Random Effect Logistic Regression Model; OR (95% CI)",
  ylab = "",
  mean.label.size = 3)+
  scale_x_log10(limits = c(-1, 10))+
  geom_vline(aes(xintercept = 1), size = 0.5, linetype = "dashed")+
  ggplot2::scale_y_discrete(labels = c("Age", "≥7 days Since Symptoms Onset","Male Sex", "Diabetes","Hypertension","COPD","Inmunosupression","CVD","Smoking Status",
                                       "Obesity","Asthma","CKD"))

Figure_1<-ggarrange(Fig1,Fig2, ncol = 2, nrow = 1, labels = c("A","B"))

ggsave(Figure_1,
       filename = "Figure2.jpg", 
       width = 55, 
       height = 25,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)





#### Risk for confusion matrix categories ####

covid2<-covid%>%filter(RESULTADO_LAB<=2)%>%filter(!is.na(RESULTADO_ANTIGENO))
covid2$resultado_final<-NULL;
covid2$resultado_final[covid2$RESULTADO_ANTIGENO==0 & covid2$RESULTADO_LAB==0]<-1
covid2$resultado_final[covid2$RESULTADO_ANTIGENO==1 & covid2$RESULTADO_LAB==0]<-2
covid2$resultado_final[covid2$RESULTADO_ANTIGENO==0 & covid2$RESULTADO_LAB==1]<-3
covid2$resultado_final[covid2$RESULTADO_ANTIGENO==1 & covid2$RESULTADO_LAB==1]<-4
covid2$resultado_final<-factor(covid2$resultado_final, labels = c("True-Negative", "False-Positive",
                                                                  "False-Negative","True-Positive"))


#Modelo de Hospitalizacion

m1<-coxme::coxme(formula = Surv(tiempo, TIPO_PACIENTE) ~ resultado_final+priv+comorb+SEXO+scale(EDAD)+(1|id), data = covid2)
summary(m1)
exp(confint(m1))

Fig1<-ggstatsplot::ggcoefstats(
  x = m1,
  exp=T,
  conf.int = TRUE,
  vline = F,
  vline.args = list(size = 1, linetype = "dashed"),
  package = "RColorBrewer",
  palette = "Dark2",
  ggtheme = ggthemes::theme_hc(),
  title = "Risk for Hospitalization in Subjects Tested for PCR/Rapid Antigen Test in Mexico",
  xlab = "Cox Proportional Hazards Regression Model with Frailty penalty function; HR (95% CI)",
  ylab = "")+
  scale_x_log10(limits = c(-1, 10))+
  ggplot2::scale_y_discrete(labels = c("False-Positive","False-Negative","True-Positive","Private HC","Number of Comorbidities","Male Sex","Age"))+
  geom_vline(aes(xintercept = 1), size = 0.5, linetype = "dashed")

#Modelo de Intubacion

m2 <-lme4::glmer(formula = INTUBADO ~ resultado_final+priv+comorb+SEXO+scale(EDAD)+(1|id),data = covid2, family = "binomial")
summary(m2)
summ(m2,exp=T,confint=T)

Fig2<-ggstatsplot::ggcoefstats(
  x = m2,
  exp=T,
  conf.int = TRUE,
  vline = F,
  exclude.intercept = T,
  vline.args = list(size = 1, linetype = "dashed"),
  package = "RColorBrewer",
  palette = "Dark2",
  ggtheme = ggthemes::theme_hc(),
  title = "Mechanical Ventilation Support Odds in Subjects Tested for PCR/Rapid Antigen Test in Mexico",
  xlab = "Random Effect Logistic Regression Model; OR (95% CI)",
  ylab = "")+
  scale_x_log10(limits = c(-1, 10))+
  ggplot2::scale_y_discrete(labels = c("False-Positive","False-Negative","True-Positive","Private HC","Number of Comorbidities","Male Sex","Age"))+
  geom_vline(aes(xintercept = 1), size = 0.5, linetype = "dashed")

#Severe Outcome
covid2$severe_outcome<-NULL;
covid2$severe_outcome[covid2$INTUBADO==1 | covid2$UCI==1 | covid2$Mortalidad==1]<-1
covid2$severe_outcome<-na.tools::na.replace(covid2$severe_outcome,0)

m3<-coxme::coxme(formula = Surv(tiempo, severe_outcome) ~ resultado_final+priv+comorb+SEXO+scale(EDAD)+(1|id), data = covid2)
summary(m3)
exp(confint(m3))

Fig3<-ggstatsplot::ggcoefstats(
  x = m3,
  exp=T,
  conf.int = TRUE,
  vline = F,
  package = "RColorBrewer",
  palette = "Dark2",
  ggtheme = ggthemes::theme_hc(),
  title = "Risk for Severe Outcome in Subjects Tested for PCR/Rapid Antigen Test in Mexico",
  xlab = "Cox Proportional Hazards Regression Model with Frailty penalty function; HR (95% CI)",
  ylab = "")+
  scale_x_log10(limits = c(-1, 10))+
  ggplot2::scale_y_discrete(labels = c("False-Positive","False-Negative","True-Positive","Private HC","Number of Comorbidities","Male Sex","Age"))+
  geom_vline(aes(xintercept = 1), size = 0.5, linetype = "dashed")

#Modelo de Mortalidad
m4<-coxme::coxme(formula = Surv(tiempo, Mortalidad) ~ resultado_final+priv+comorb+SEXO+scale(EDAD)+(1|id), data = covid2)
summary(m4)
exp(confint(m4))

Fig4<-ggstatsplot::ggcoefstats(
  x = m4,
  exp=T,
  conf.int = TRUE,
  package = "RColorBrewer",
  vline = F,
  palette = "Dark2",
  ggtheme = ggthemes::theme_hc(),
  title = "Risk for Letality in Subjects Tested for PCR/Rapid Antigen Test in Mexico",
  xlab = "Cox Proportional Hazards Regression Model with Frailty penalty function; HR (95% CI)",
  ylab = "")+
  scale_x_log10(limits = c(-1, 10))+
  ggplot2::scale_y_discrete(labels = c("False-Positive","False-Negative","True-Positive","Private HC","Number of Comorbidities","Male Sex","Age"))+
  geom_vline(aes(xintercept = 1), size = 0.5, linetype = "dashed")

Figure_2<-ggarrange(Fig1,Fig2,Fig3,Fig4, ncol = 2, nrow = 2, labels = c("A","B","C","D"))
Figure_2

ggsave(Figure_2,
       filename = "Figure3.jpg", 
       width = 55, 
       height = 25,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)

#### Predictors of outcomes in positive Ag tests #####

covid3<-covid%>%filter(is.na(RESULTADO_LAB))%>%filter(!is.na(RESULTADO_ANTIGENO))

#Modelo de Hospitalización
m1<-coxme::coxme(formula = Surv(tiempo, TIPO_PACIENTE) ~ scale(EDAD)+SEXO+DIABETES+HIPERTENSION+EPOC+INMUSUPR+CARDIOVASCULAR+
                   TABAQUISMO+OBESIDAD+ASMA+RENAL_CRONICA+(1|id), data = covid3[covid3$RESULTADO_ANTIGENO==1,])
summary(m1)
exp(confint(m1))

Fig1<-ggstatsplot::ggcoefstats(
  x = m1,
  exp=T,
  conf.int = TRUE,
  package = "RColorBrewer",
  vline = F,
  palette = "Dark2",
  ggtheme = ggthemes::theme_hc(),
  title = "Risk for Hospitalization in Subjects with Positive Rapid Antigen Test",
  xlab = "Random Effects Cox Proportional Hazards Regression Model; HR (95% CI)",
  ylab = "")+
  scale_x_log10(limits = c(-1, 20))+
  ggplot2::scale_y_discrete(labels = c("Age", "Male Sex", "Diabetes", "Hypertension", "COPD", "Immunosuppression","CVD","Smoking Status","Obesity","Asthma","CKD"))+
  geom_vline(aes(xintercept = 1), size = 0.5, linetype = "dashed")

#Modelo de Intubacion

m2 <- lme4::glmer(
  formula = INTUBADO~scale(EDAD)+SEXO+DIABETES+HIPERTENSION+EPOC+CARDIOVASCULAR+
    TABAQUISMO+OBESIDAD+RENAL_CRONICA+(1|id),
  data = covid3[covid3$RESULTADO_ANTIGENO==1,],
  family = "binomial")
summary(m2)
summ(m2,exp=T,confint=T)

Fig2<-ggstatsplot::ggcoefstats(
  x = m2,
  exp=T,
  conf.int = TRUE,
  vline = F,
  exclude.intercept = T,
  vline.args = list(size = 1, linetype = "dashed"),
  package = "RColorBrewer",
  palette = "Dark2",
  ggtheme = ggthemes::theme_hc(),
  title = "Mechanical Ventilation Support Odds in Subjects with Positive Rapid Antigen Test",
  xlab = "Random Effect Logistic Regression Model; OR (95% CI)",
  ylab = "")+
  scale_x_log10(limits = c(-1, 20))+
  ggplot2::scale_y_discrete(labels = c("Age", "Male Sex", "Diabetes", "Hypertension", "COPD", "Immunosuppression","CVD","Smoking Status","Obesity","Asthma","CKD"))+
  geom_vline(aes(xintercept = 1), size = 0.5, linetype = "dashed")


#Severe Outcome
covid3$severe_outcome<-NULL;
covid3$severe_outcome[covid3$INTUBADO==1 | covid3$UCI==1 | covid3$Mortalidad==1]<-1
covid3$severe_outcome<-na.tools::na.replace(covid3$severe_outcome,0)

m3<-coxme::coxme(formula = Surv(tiempo, severe_outcome) ~ scale(EDAD)+SEXO+DIABETES+HIPERTENSION+EPOC+INMUSUPR+CARDIOVASCULAR+
                   TABAQUISMO+OBESIDAD+ASMA+RENAL_CRONICA+(1|id), data = covid3[covid3$RESULTADO_ANTIGENO==1,])
summary(m3)
exp(confint(m3))

Fig3<-ggstatsplot::ggcoefstats(
  x = m3,
  exp=T,
  conf.int = TRUE,
  vline = F,
  package = "RColorBrewer",
  palette = "Dark2",
  ggtheme = ggthemes::theme_hc(),
  title = "Risk for Severe Outcome in Subjects with Positive Rapid Antigen Test",
  xlab = "Random Effects Cox Proportional Hazards Regression Model; HR (95% CI)",
  ylab = "")+
  scale_x_log10(limits = c(-1, 20))+
  ggplot2::scale_y_discrete(labels = c("Age", "Male Sex", "Diabetes", "Hypertension", "COPD", "Immunosuppression","CVD","Smoking Status","Obesity","Asthma","CKD"))+
  geom_vline(aes(xintercept = 1), size = 0.5, linetype = "dashed")
Fig3
#Modelo de Mortalidad
m4<-coxme::coxme(formula = Surv(tiempo, Mortalidad) ~ scale(EDAD)+SEXO+DIABETES+HIPERTENSION+EPOC+INMUSUPR+CARDIOVASCULAR+
                   TABAQUISMO+OBESIDAD+ASMA+RENAL_CRONICA+(1|id), data = covid3[covid3$RESULTADO_ANTIGENO==1,])
summary(m4)
exp(confint(m4))

Fig4<-ggstatsplot::ggcoefstats(
  x = m4,
  exp=T,
  conf.int = TRUE,
  package = "RColorBrewer",
  vline = F,
  palette = "Dark2",
  ggtheme = ggthemes::theme_hc(),
  title = "Risk for Letality in Subjects with Positive Rapid Antigen Test",
  xlab = "Random Effect Logistic Regression Model; OR (95% CI)",
  ylab = "")+
  scale_x_log10(limits = c(-1, 10))+
  ggplot2::scale_y_discrete(labels = c("Age", "Male Sex", "Diabetes", "Hypertension", "COPD", "Immunosuppression","CVD","Smoking Status","Obesity","Asthma","CKD"))+
  geom_vline(aes(xintercept = 1), size = 0.5, linetype = "dashed")

Figure_3<-ggarrange(Fig1,Fig2,Fig3,Fig4, ncol = 2, nrow = 2, labels = c("A","B","C","D"))
Figure_3

ggsave(Figure_3,
       filename = "Figure4.jpg", 
       width = 55, 
       height = 35,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)



#### Supplementary table ####
covid$muestra<-covid$TOMA_MUESTRA_ANTIGENO+2*covid$TOMA_MUESTRA_LAB

covid4<-covid %>% filter(covid==1, muestra!=0, muestra !=3)
covid_pcr<-covid4%>%filter(muestra==1)
covid_ag<-covid4%>%filter(muestra==2)

# Age
pcr0<-paste0(round(mean(covid_pcr$EDAD), 1),"±",round(sd(covid_pcr$EDAD),1))
ag0<-paste0(round(mean(covid_ag$EDAD), 1),"±",round(sd(covid_ag$EDAD),1))

#Sexo
pcr1<-table(covid_pcr$SEXO)[2];pcr1.1<-(table(covid_pcr$SEXO)/nrow(covid_pcr))[2]
ag1<-table(covid_ag$SEXO)[2];ag1.1<-(table(covid_ag$SEXO)/nrow(covid_ag))[2]

#Diabetes
pcr2<-table(covid_pcr$DIABETES)[2];pcr2.1<-(table(covid_pcr$DIABETES)/nrow(covid_pcr))[2]
ag2<-table(covid_ag$DIABETES)[2];ag2.1<-(table(covid_ag$DIABETES)/nrow(covid_ag))[2]

#EPOC
pcr3<-table(covid_pcr$EPOC)[2];pcr3.1<-(table(covid_pcr$EPOC)/nrow(covid_pcr))[2]
ag3<-table(covid_ag$EPOC)[2];ag3.1<-(table(covid_ag$EPOC)/nrow(covid_ag))[2]

#ASMA
pcr4<-table(covid_pcr$ASMA)[2];pcr4.1<-(table(covid_pcr$ASMA)/nrow(covid_pcr))[2]
ag4<-table(covid_ag$ASMA)[2];ag4.1<-(table(covid_ag$ASMA)/nrow(covid_ag))[2]

#INMUSUPR
pcr5<-table(covid_pcr$INMUSUPR)[2];pcr5.1<-(table(covid_pcr$INMUSUPR)/nrow(covid_pcr))[2]
ag5<-table(covid_ag$INMUSUPR)[2];ag5.1<-(table(covid_ag$INMUSUPR)/nrow(covid_ag))[2]

#HIPERTENSION
pcr6<-table(covid_pcr$HIPERTENSION)[2];pcr6.1<-(table(covid_pcr$HIPERTENSION)/nrow(covid_pcr))[2]
ag6<-table(covid_ag$HIPERTENSION)[2];ag6.1<-(table(covid_ag$HIPERTENSION)/nrow(covid_ag))[2]

#OTRA_COM
pcr7<-table(covid_pcr$OTRA_COM)[2];pcr7.1<-(table(covid_pcr$OTRA_COM)/nrow(covid_pcr))[2]
ag7<-table(covid_ag$OTRA_COM)[2];ag7.1<-(table(covid_ag$OTRA_COM)/nrow(covid_ag))[2]

#CARDIOVASCULAR
pcr8<-table(covid_pcr$CARDIOVASCULAR)[2];pcr8.1<-(table(covid_pcr$CARDIOVASCULAR)/nrow(covid_pcr))[2]
ag8<-table(covid_ag$CARDIOVASCULAR)[2];ag8.1<-(table(covid_ag$CARDIOVASCULAR)/nrow(covid_ag))[2]

#OBESIDAD
pcr9<-table(covid_pcr$OBESIDAD)[2];pcr9.1<-(table(covid_pcr$OBESIDAD)/nrow(covid_pcr))[2]
ag9<-table(covid_ag$OBESIDAD)[2];ag9.1<-(table(covid_ag$OBESIDAD)/nrow(covid_ag))[2]

#RENAL_CRONICA
pcr10<-table(covid_pcr$RENAL_CRONICA)[2];pcr10.1<-(table(covid_pcr$RENAL_CRONICA)/nrow(covid_pcr))[2]
ag10<-table(covid_ag$RENAL_CRONICA)[2];ag10.1<-(table(covid_ag$RENAL_CRONICA)/nrow(covid_ag))[2]

#TABAQUISMO
pcr11<-table(covid_pcr$TABAQUISMO)[2];pcr11.1<-(table(covid_pcr$TABAQUISMO)/nrow(covid_pcr))[2]
ag11<-table(covid_ag$TABAQUISMO)[2];ag11.1<-(table(covid_ag$TABAQUISMO)/nrow(covid_ag))[2]

#NEUMONIA
pcr12<-table(covid_pcr$NEUMONIA)[2];pcr12.1<-(table(covid_pcr$NEUMONIA)/nrow(covid_pcr))[2]
ag12<-table(covid_ag$NEUMONIA)[2];ag12.1<-(table(covid_ag$NEUMONIA)/nrow(covid_ag))[2]

#HOSPITALIZATION
pcr13<-table(covid_pcr$TIPO_PACIENTE)[2];pcr13.1<-(table(covid_pcr$TIPO_PACIENTE)/nrow(covid_pcr))[2]
ag13<-table(covid_ag$TIPO_PACIENTE)[2];ag13.1<-(table(covid_ag$TIPO_PACIENTE)/nrow(covid_ag))[2]

#ICU
pcr14<-table(covid_pcr$UCI)[2];pcr14.1<-(table(covid_pcr$UCI)/nrow(covid_pcr))[2]
ag14<-table(covid_ag$UCI)[2];ag14.1<-(table(covid_ag$UCI)/nrow(covid_ag))[2]

#DEATH
pcr15<-table(covid_pcr$Mortalidad)[2];pcr15.1<-(table(covid_pcr$Mortalidad)/nrow(covid_pcr))[2]
ag15<-table(covid_ag$Mortalidad)[2];ag15.1<-(table(covid_ag$Mortalidad)/nrow(covid_ag))[2]

#Days to assessment
pcr16<-paste0(round(median(covid_pcr$tiempo), 1)," (",round(quantile(covid_pcr$tiempo)[2],1),"-",round(quantile(covid_pcr$tiempo)[4],1),")")
ag16<-paste0(round(median(covid_ag$tiempo), 1)," (",round(quantile(covid_ag$tiempo)[2],1),"-",round(quantile(covid_ag$tiempo)[4],1),")")

pcr_1<-paste0(pcr1," ","(",round(pcr1.1*100,1),")");ag_1<-paste0(ag1," ","(",round(ag1.1*100,1),")")
pcr_2<-paste0(pcr2," ","(",round(pcr2.1*100,1),")");ag_2<-paste0(ag2," ","(",round(ag2.1*100,1),")")
pcr_3<-paste0(pcr3," ","(",round(pcr3.1*100,1),")");ag_3<-paste0(ag3," ","(",round(ag3.1*100,1),")")
pcr_4<-paste0(pcr4," ","(",round(pcr4.1*100,1),")");ag_4<-paste0(ag4," ","(",round(ag4.1*100,1),")")
pcr_5<-paste0(pcr5," ","(",round(pcr5.1*100,1),")");ag_5<-paste0(ag5," ","(",round(ag5.1*100,1),")")
pcr_6<-paste0(pcr6," ","(",round(pcr6.1*100,1),")");ag_6<-paste0(ag6," ","(",round(ag6.1*100,1),")")
pcr_7<-paste0(pcr7," ","(",round(pcr7.1*100,1),")");ag_7<-paste0(ag7," ","(",round(ag7.1*100,1),")")
pcr_8<-paste0(pcr8," ","(",round(pcr8.1*100,1),")");ag_8<-paste0(ag8," ","(",round(ag8.1*100,1),")")
pcr_9<-paste0(pcr9," ","(",round(pcr9.1*100,1),")");ag_9<-paste0(ag9," ","(",round(ag9.1*100,1),")")
pcr_10<-paste0(pcr10," ","(",round(pcr10.1*100,1),")");ag_10<-paste0(ag10," ","(",round(ag10.1*100,1),")")
pcr_11<-paste0(pcr11," ","(",round(pcr11.1*100,1),")");ag_11<-paste0(ag11," ","(",round(ag11.1*100,1),")")
pcr_12<-paste0(pcr12," ","(",round(pcr12.1*100,1),")");ag_12<-paste0(ag12," ","(",round(ag12.1*100,1),")")
pcr_13<-paste0(pcr13," ","(",round(pcr13.1*100,1),")");ag_13<-paste0(ag13," ","(",round(ag13.1*100,1),")")
pcr_14<-paste0(pcr14," ","(",round(pcr14.1*100,1),")");ag_14<-paste0(ag14," ","(",round(ag14.1*100,1),")")
pcr_15<-paste0(pcr15," ","(",round(pcr15.1*100,1),")");ag_15<-paste0(ag15," ","(",round(ag15.1*100,1),")")


covid_tests<-covid4
P_0<-round(t.test(covid_tests$EDAD~covid_tests$muestra)$p.value,3)
P_1<-round(chisq.test(covid_tests$SEXO,covid_tests$muestra)$p.value,3);P_1
P_2<-round(chisq.test(covid_tests$DIABETES,covid_tests$muestra)$p.value,3);P_2
P_3<-round(chisq.test(covid_tests$EPOC,covid_tests$muestra)$p.value,3);P_3
P_4<-round(chisq.test(covid_tests$ASMA,covid_tests$muestra)$p.value,3);P_4
P_5<-round(chisq.test(covid_tests$INMUSUPR,covid_tests$muestra)$p.value,3);P_5
P_6<-round(chisq.test(covid_tests$HIPERTENSION,covid_tests$muestra)$p.value,3);P_6
P_7<-round(chisq.test(covid_tests$OTRA_COM,covid_tests$muestra)$p.value,3);P_7
P_8<-round(chisq.test(covid_tests$CARDIOVASCULAR,covid_tests$muestra)$p.value,3);P_8
P_9<-round(chisq.test(covid_tests$OBESIDAD,covid_tests$muestra)$p.value,3);P_9
P_10<-round(chisq.test(covid_tests$RENAL_CRONICA,covid_tests$muestra)$p.value,3);P_10
P_11<-round(chisq.test(covid_tests$TABAQUISMO,covid_tests$muestra)$p.value,3);P_11
P_12<-round(chisq.test(covid_tests$NEUMONIA,covid_tests$muestra)$p.value,3);P_12
P_13<-round(chisq.test(covid_tests$TIPO_PACIENTE,covid_tests$muestra)$p.value,3);P_13
P_14<-round(chisq.test(covid_tests$UCI,covid_tests$muestra)$p.value,3);P_14
P_15<-round(chisq.test(covid_tests$Mortalidad,covid_tests$muestra)$p.value,3);P_15
P_16<-wilcox.test(covid_tests$tiempo, covid_tests$muestra)$p.value;P_16

GGG<-base::format.pval(c(P_0,P_1,P_2,P_3,P_4,P_5,P_6,P_7,P_8,P_9,P_10,P_11,P_12,P_13,P_14,P_15, P_16), eps = .001, digits = 2)

AM_JOV_COVID<-data.frame("bold(Parameter)"=c("Age (years)","Male sex (%)","Diabetes (%)","COPD (%)","Asthma (%)","Immunosuppression (%)","Hypertension (%)",
                                             "Other (%)","CVD (%)","Obesity","CKD (%)","Smoking (%)","Pneumonia(%)",
                                             "Hospitalization(%)","ICU admission (%)", "Death (%)", "Time to assessment* (days)"),
                         "RT-PCR test"=c(pcr0, pcr_1,pcr_2, pcr_3, pcr_4, pcr_5, pcr_6, pcr_7, pcr_8, pcr_9, pcr_10, pcr_11, pcr_12, pcr_13, pcr_14, pcr_15, pcr16),
                         "Rapid Ag-T"=c(ag0,ag_1,ag_2, ag_3, ag_4, ag_5, ag_6, ag_7, ag_8, ag_9, ag_10, ag_11, ag_12, ag_13, ag_14, ag_15, ag16),
                         "p-value"=GGG)

colnames(AM_JOV_COVID)<-c("Parameters",paste0("RT-PCR positive SARS-CoV-2 \nn=",nrow(covid_pcr)),
                          paste0("Ag-T positive SARS-CoV2- \nn=",nrow(covid_ag)),"p-value")

AM_JOV_COVID<-flextable(AM_JOV_COVID,cwidth = 0.5*ncol(AM_JOV_COVID)) %>% autofit()


doc <- read_docx() %>%
  body_add_flextable(value = AM_JOV_COVID, split = TRUE) %>%
  body_end_section_landscape() %>% # a landscape section is ending here
  print( target = "SuppTable1.docx" )



