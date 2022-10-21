# Databricks notebook source
library(data.table)
# library(sapbio.ds.ml)
library(ggplot2)
# library(lme4)
library(MASS)
# library(broom.mixed)

# COMMAND ----------

library(ordinal)

# COMMAND ----------

# MAGIC %md # load data

# COMMAND ----------

# df=fread('/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/data/std_log_impute_raw_data_w_pheno_091922.csv')

df=fread('/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/data/raw_data_w_pheno_091922.csv')

df[1:3,1:4]

# COMMAND ----------

df[is.na(df)]=0
pheno_list=colnames(df)[!startsWith(colnames(df), 'rLC')]
mtb_list=colnames(df)[startsWith(colnames(df),'rLC')]
# pheno_list
print(dim(df))
print(length(mtb_list))
mdata <- melt(df, id=c(pheno_list))
print(unique(mdata$Sarnat_T1))
print(unique(mdata$Sarnat_T2))

# COMMAND ----------


mdata$mtb_id=as.character(unlist(mdata$variable))
mdata=as.data.frame(mdata)
mdata$mtb_id=as.character(unlist(mdata$variable))
names(mdata)[names(mdata) == "Study ID (Neonate)"] <- "StudyID"
names(mdata)[names(mdata) == "Study ID (Neonate)"] <- "StudyID"
names(mdata)[names(mdata) == "Overall highest Thompson score"] <- "WorstThompson"
names(mdata)[names(mdata) == "Overall worst Sarnat grade"] <- "WorstSarnat"
names(mdata)[names(mdata) == "Sample type and description"] <- "SampleType"
names(mdata)[names(mdata) == "First time-point worst grade recorded"] <- "WorstGradeTime"
names(mdata)[names(mdata) == "Unique sample ID_x"] <- "USmpleID"
# names(mdata)[names(mdata) == "variable"] <- "mtb"
names(mdata)[names(mdata) == "value"] <- "level"

mdata[mdata$Sarnat_T1=="Moderate (2)",'Sarnat_T1']=2
mdata[mdata$Sarnat_T1=="Severe (3)",'Sarnat_T1']=3
mdata[mdata$Sarnat_T1=="Mild (1)",'Sarnat_T1']=1
mdata[mdata$Sarnat_T1=='-','Sarnat_T1']='-1'
mdata[mdata$Sarnat_T1=='','Sarnat_T1']='-1'

mdata[mdata$Sarnat_T2=="Moderate (2)",'Sarnat_T2']=2
mdata[mdata$Sarnat_T2=="Severe (3)",'Sarnat_T2']=3
mdata[mdata$Sarnat_T2=="Mild (1)",'Sarnat_T2']=1
mdata[mdata$Sarnat_T2=="Normal (0)",'Sarnat_T2']='0'
mdata[mdata$Sarnat_T2=="Normal",'Sarnat_T2']='0'
mdata[mdata$Sarnat_T2=='-','Sarnat_T2']='-1'
mdata[mdata$Sarnat_T2=='','Sarnat_T2']='-1'

mdata[mdata$WorstSarnat=="Moderate",'WorstSarnat']=2
mdata[mdata$WorstSarnat=="Severe",'WorstSarnat']=3
mdata[mdata$WorstSarnat=="Mild",'WorstSarnat']=1
mdata[mdata$WorstSarnat=="Not done",'WorstSarnat']='-1'
mdata[mdata$WorstSarnat=="Data not available on source",'WorstSarnat']='-1'
mdata[mdata$WorstSarnat=='','WorstSarnat']='-1'

mdata$Visit=mdata$Timepoint
mdata[mdata$Visit=="T1",'Visit']=1
mdata[mdata$Visit=="T2",'Visit']=2


mdata=setDT(mdata) 
mdata$Sarnat_T1=as.numeric(unlist(mdata$Sarnat_T1))
mdata$Sarnat_T2=as.numeric(unlist(mdata$Sarnat_T2))
print(unique(mdata$Sarnat_T1))
print(unique(mdata$Sarnat_T2))


# COMMAND ----------

# MAGIC %md # load results

# COMMAND ----------

# MAGIC %md ## severity linear mixed

# COMMAND ----------

# severity_results = readRDS('/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/lmm_rstudio_results_w_warningMsg_msg.RData')
# std_severity_results = readRDS('/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/std_lmm_rstudio_results_w_warningMsg_msg.RData')
std_severity_results = readRDS('/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/results/lmm_std_log_impute_raw_Res_lmm_timepoint_cooled.rdata')

# severity_output=dcast(severity_results,formula = mtb_id+warnings_full+msg~ term, value.var = c('estimate','std.error','statistic','p.value'))
# severity_std_output=dcast(std_severity_results,formula = mtb_id+warnings_full+msg~ term, value.var = c('estimate','std.error','statistic','p.value'))
severity_std_log_output=dcast(std_severity_results,formula = mtb_id + conv_messages + warnings + gradient + d_messages + Hessian ~ term, value.var = c('estimate','std.error','statistic','p.value'))

# print(dim(severity_output))
# print(dim(severity_std_output))
# print(dim(severity_std_log_output))

# COMMAND ----------

std_severity_results |> colnames()

# COMMAND ----------

severity_std_log_output$conv_messages |> unique()

# COMMAND ----------

severity_std_log_output$d_messages |> unique()

# COMMAND ----------

# MAGIC %md ## sarnat time 1 time 2

# COMMAND ----------

results_t1 = readRDS('/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/results/t1_ME_clm_std_log_impute_raw_data_w_pheno_092822.rdata')
results_t2 = readRDS('/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/results/t2_ME_clm_std_log_impute_raw_data_w_pheno_092822.rdata')

severity_t1=dcast(data = results_t1,formula = mtb_id + warnings_full ~ term, value.var = c('estimate','std.error','statistic','p.value'))
severity_t2=dcast(data = results_t2,formula = mtb_id + warnings_full ~ term, value.var = c('estimate','std.error','statistic','p.value'))
# severity_t2

# COMMAND ----------

# c_results_t1 = readRDS('/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/sexGA_t1_std_ordinal_clm_rstudio_results_w_warningMsg_msg.RData')
# c_results_t2 = readRDS('dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/sexGA_t2_std_ordinal_clm_rstudio_results_w_warningMsg_msg.RData')


# c_results_t1
# covariate_severity_t1=dcast(covariate_results_t1,formula = mtb_id+warnings_full ~ term, value.var = c('estimate','std.error','statistic','p.value'))
# covariate_severity_t2=dcast(covariate_results_t2,formula = mtb_id+warnings_full ~ term, value.var = c('estimate','std.error','statistic','p.value'))

# COMMAND ----------

# MAGIC %md ## peak severity ordinal

# COMMAND ----------

results_peak = readRDS('/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/results/peak_ME_clm_std_log_impute_raw_data_w_pheno_092822.rdata')

# COMMAND ----------

severity_peak=dcast(results_peak,formula = mtb_id+warnings_full ~ term, value.var = c('estimate','std.error','statistic','p.value'))

# COMMAND ----------

# MAGIC %md ## thomspon

# COMMAND ----------

tmp_results_t1 = readRDS('/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/results/t1_thompson_ME_clm_std_log_impute_raw_data_w_pheno_092822.rdata')
tmp_results_t2 = readRDS('/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/results/t2_thompson_ME_clm_std_log_impute_raw_data_w_pheno_092822.rdata')
tmp_results_t2
tmp_severity_t1=dcast(tmp_results_t1,formula = mtb_id+warnings_full ~ term, value.var = c('estimate','std.error','statistic','p.value'))
tmp_severity_t2=dcast(tmp_results_t2,formula = mtb_id+warnings_full ~ term, value.var = c('estimate','std.error','statistic','p.value'))

# COMMAND ----------

# MAGIC %md # apgar

# COMMAND ----------

Apgar_m1_lm_data = readRDS('/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/results/Apgar_m1_lm_std_log_impute_raw_data_w_pheno_092822.rdata')
Apgar_m5_lm_data = readRDS('/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/results/Apgar_m5_lm_std_log_impute_raw_data_w_pheno_092822.rdata')
Apgar_m10_lm_data = readRDS('/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/results/Apgar_m10_lm_std_log_impute_raw_data_w_pheno_092822.rdata')


Apgar_m1_severity_t1=dcast(Apgar_m1_lm_data,formula = mtb_id+warnings_full ~ term, value.var = c('estimate','std.error','statistic','p.value'))

Apgar_m5_severity_t1=dcast(Apgar_m5_lm_data,formula = mtb_id+warnings_full ~ term, value.var = c('estimate','std.error','statistic','p.value'))

Apgar_m10_severity_t1=dcast(Apgar_m10_lm_data,formula = mtb_id+warnings_full ~ term, value.var = c('estimate','std.error','statistic','p.value'))


# COMMAND ----------

# MAGIC %md # severity plots

# COMMAND ----------

# print(dim(severity_std_output[severity_std_output$warnings_full!='None',]))
print(dim(severity_std_log_output[severity_std_log_output$warnings_full!='No',]))

# COMMAND ----------

data=as.data.frame(severity_std_log_output)
# data=as.data.frame(output)
# data=data[data$warnings_full=='None',]
data['logp']=-log10(data$p.value_Sarnat)
data['Type']='Below threshold'
data[data$logp>3 & data$estimate_Sarnat>0 ,'Type']='Elevated above threshold'
data[data$logp>3 & data$estimate_Sarnat<0 ,'Type']='Reduced above threshold'

# COMMAND ----------

myplot=ggplot(data, aes(x=estimate_Sarnat, y=logp,color=Type))+geom_point(alpha = 0.4 )
myplot
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlim(-2, 2) + 
scale_color_manual(values=c('grey','#039fbe','#cf1578'))+ 
theme(legend.text = element_text(size = 15))+
 theme(axis.title.y = element_text(size = 20))  + theme(axis.title.x = element_text(size = 20)) +
 theme(axis.text.y = element_text(size = 15)) + theme(axis.text.x = element_text(size = 15)) 

# COMMAND ----------

p_th=3
print(length(data[data$logp>p_th & data$estimate_Sarnat<0 ,'mtb_id']))
print(length(data[data$logp>p_th & data$estimate_Sarnat>0 ,'mtb_id']))
reduced_lmm=data[data$logp>p_th & data$estimate_Sarnat<0 ,'mtb_id']
elevated_lmm=data[data$logp>p_th & data$estimate_Sarnat>0 ,'mtb_id']
lmm=data[data$logp>3,'mtb_id']

data[data$logp>3 & data$estimate_Sarnat<0 & data$warnings_full=='None','mtb_id']

# COMMAND ----------

# mtb='rLC_neg_mtb_1293276' # Optimization stopped
# mtb='rLC_pos_mtb_629943' # best down
mtb='rLC_pos_mtb_869671' # best down
# mtb='rLC_neg_mtb_1375224' # 2nd up
# mtb='rLC_neg_mtb_3724398' # best up
# mtb='rLC_neg_mtb_3755909' #peak
# mtb = 'rLC_neg_mtb_1583136' #warning molecule


# mtb = 'rLC_neg_mtb_1583136' #warning molecule
# mtb = 'rLC_neg_mtb_1583136' #warning molecule
# mtb = 'rLC_neg_mtb_1583136' #warning molecule

test_data2=mdata[mtb_id %in% mtb] 
test_data2=as.data.frame(test_data2)
V1sit1 = test_data2[test_data2$Timepoint == 'T1',c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat_T1','Visit')]
colnames(V1sit1)=c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat','Visit')
V1sit2 = test_data2[test_data2$Timepoint == 'T2',c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat_T2','Visit')]
colnames(V1sit2)=c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat','Visit')
df1 = rbind(V1sit1,V1sit2)
df1=df1[df1$Sarnat>=1,]
print(dim(df1))
df1$Sarnat=as.factor(df1$Sarnat)


# COMMAND ----------

myplot=ggplot(df1, aes(x=Sarnat, y=level, fill=Patient_ID)) +
# myplot=ggplot(df1, aes(x=Sarnat, y=level, fill=Timepoint)) +
  geom_boxplot()
# Add dots
myplot + geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5,
                 position=position_dodge(1))+

# scale_color_manual(values=c('#039fbe','#cf1578'))+
ggtitle(mtb)+ 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black")) +
# # # # #   xlim(-1000, 1000) + 
 theme(legend.position="none") +

# theme(legend.text = element_text(size = 20))+
 theme(axis.title.y = element_text(size = 25))  + theme(axis.title.x = element_text(size = 25)) +
 theme(axis.text.y = element_text(size = 20)) + theme(axis.text.x = element_text(size =20)) + theme(
plot.title = element_text(color="black", size=25, face="bold"))


# COMMAND ----------

# MAGIC %md # log severity plot

# COMMAND ----------

severity_std_log_output

# COMMAND ----------

logp_th=-log10(7.774227e-06)
data=as.data.frame(severity_std_log_output)
# data=as.data.frame(output)
# data=data[data$warnings_full=='None',]
data['logp']=-log10(data$p.value_Sarnat)
data['Type']='Below threshold'
data[data$logp>logp_th & data$estimate_Sarnat>0 ,'Type']='Elevated above threshold'
data[data$logp>logp_th & data$estimate_Sarnat<0 ,'Type']='Reduced above threshold'

print(length(data[data$logp>logp_th & data$estimate_Sarnat<0 ,'mtb_id']))
print(length(data[data$logp>logp_th & data$estimate_Sarnat>0 ,'mtb_id']))
reduced_log_lmm=data[data$logp>logp_th & data$estimate_Sarnat<0 ,'mtb_id']
elevated_log_lmm=data[data$logp>logp_th & data$estimate_Sarnat>0 ,'mtb_id']
log_lmm=data[data$logp>logp_th,'mtb_id']

myplot=ggplot(data, aes(x=estimate_Sarnat, y=logp,color=Type))+geom_point(alpha = 0.4 )
myplot
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
#   xlim(-2, 2) + 
scale_color_manual(values=c('grey','#039fbe','#cf1578'))+ 
theme(legend.text = element_text(size = 15))+
 theme(axis.title.y = element_text(size = 20))  + theme(axis.title.x = element_text(size = 20)) +
 theme(axis.text.y = element_text(size = 15)) + theme(axis.text.x = element_text(size = 15)) 

# COMMAND ----------

# library(tidyverse)

data=as.data.frame(severity_std_log_output)
t2_p_threshold=-log10(0.05/7)
# data=as.data.frame(output)
# data=data[data$warnings_full=='None',]
data['logp']=-log10(data$p.value_Sarnat)
data['Type']='Below threshold'
# data[data$logp>t2_p_threshold & data$mtb_id == 'rLC_neg_mtb_3755909' ,'Type']='Peak above threshold'
# data[data$mtb_id  == 'rLC_neg_mtb_3755909','Type']='Peak'
# peak_data=data[data$mtb_id  == 'rLC_neg_mtb_3755909',]

data[data$logp>t2_p_threshold & data$mtb_id %in% t2 ,'T2 levated above lmm threshold']
data[data$mtb_id %in% t2,'Type']='T2 elevated'
t2_data=data[data$mtb_id %in% t2,]

ggplot(data, aes(x=estimate_Sarnat, y=logp)) +
  theme_classic() +
  geom_point(
    mapping = aes(colour = Type),
    shape = 1,
#     size = 3,
#     stroke = 2,
    alpha = 5 / 6
    ) + xlim(-2, 2) +
scale_color_manual(values=c('grey','#cf1578'))+
  geom_point(t2_data, 
    mapping = aes(x=estimate_Sarnat, y=logp,fill = Type, colour = NA),
#     size = 2.88,
    alpha = 5 /6,
    shape = 21
  ) + geom_hline(yintercept=-log10(0.05/7))


# COMMAND ----------

logp_th=-log10(7.774227e-06)
print(length(data[data$logp>logp_th & data$estimate_Sarnat<0 ,'mtb_id']))
print(length(data[data$logp>logp_th & data$estimate_Sarnat>0 ,'mtb_id']))
reduced_log_lmm=data[data$logp>logp_th & data$estimate_Sarnat<0 ,'mtb_id']
elevated_log_lmm=data[data$logp>logp_th & data$estimate_Sarnat>0 ,'mtb_id']
log_lmm=data[data$logp>logp_th,'mtb_id']

# data[data$logp>3 & data$estimate_Sarnat<0 & data$warnings_full!='No','mtb_id']
data[data$logp>logp_th & data$estimate_Sarnat>0 ,c('mtb_id','estimate_Sarnat','p.value_Sarnat')]
# data[data$logp>3 & data$estimate_Sarnat<0, c('mtb_id','estimate_Sarnat','p.value_Sarnat')]
log_lmm

# COMMAND ----------

# MAGIC %md ## lmm box plot

# COMMAND ----------


# mtb='rLC_neg_mtb_3466281' # best down

# mtb='rLC_pos_mtb_3306281' # best up
# mtb='rLC_neg_mtb_3630952' # best down warning
# mtb = 'rLC_neg_mtb_3920887' #best up warning molecule

# mtb='rLC_neg_mtb_1339474' #overlap
# mtb='rLC_pos_mtb_3273673'
# the 7 overlap
mtb='rLC_neg_mtb_1823778'
# mtb='rLC_neg_mtb_1825916'
# mtb='rLC_neg_mtb_1826191'
# mtb='rLC_neg_mtb_1846759'
# mtb='rLC_neg_mtb_1602111'
# mtb='rLC_neg_mtb_1640306'
# mtb='rLC_neg_mtb_1871859'


# peak
mtb='rLC_neg_mtb_3755909'

# the 7 lmm
mtb='rLC_neg_mtb_1846759'
mtb='rLC_neg_mtb_1825916'
mtb='rLC_pos_mtb_3273673'
mtb='rLC_pos_mtb_3306281'
mtb='rLC_neg_mtb_1823778'
mtb='rLC_neg_mtb_1826191'
mtb='rLC_pos_mtb_3144569'




test_data2=mdata[mtb_id %in% mtb] 
test_data2=as.data.frame(test_data2)
V1sit1 = test_data2[test_data2$Timepoint == 'T1',c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat_T1','Visit')]
colnames(V1sit1)=c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat','Visit')
V1sit2 = test_data2[test_data2$Timepoint == 'T2',c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat_T2','Visit')]
colnames(V1sit2)=c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat','Visit')
df1 = rbind(V1sit1,V1sit2)
df1=df1[df1$Sarnat>=1,]
print(dim(df1))
df1$Sarnat=as.factor(df1$Sarnat)

# COMMAND ----------

myplot=ggplot(df1, aes(x=Sarnat, y=level, fill=Timepoint)) +
geom_boxplot(position = position_dodge(1), outlier.shape = NA)+
      geom_point(aes(fill =  factor(Timepoint),  group = factor(Timepoint)), 
#                  color="black", 
                 alpha  =0.5, size=3, 
                 position = position_jitterdodge(jitter.width = .1, dodge.width = 1)) +
#       scale_shape_manual (values = c(21,23) ) +
      guides(fill = guide_legend(override.aes = list(shape = NA) ) )

myplot+
# ggtitle(mtb)+ xlab("Dose (mg)") +
ylab(mtb)+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black")) +
# # # # #   xlim(-1000, 1000) + 

theme(legend.text = element_text(size = 20))+
 theme(axis.title.y = element_text(size = 25))  + theme(axis.title.x = element_text(size = 25)) +
 theme(axis.text.y = element_text(size = 20)) + theme(axis.text.x = element_text(size =20)) + theme(
plot.title = element_text(color="black", size=25, face="bold"))


# COMMAND ----------

# MAGIC %md # t1

# COMMAND ----------

data=as.data.frame(severity_t1)
t1_p_threshold=-log10(3.07662e-05)
# data=as.data.frame(output)
# data=data[data$warnings_full=='None',]
data['logp']=-log10(data$p.value_level)
data['Type']='Below threshold'
data[data$logp>t1_p_threshold & data$estimate_level>0 ,'Type']='Elevated above threshold'
data[data$logp>t1_p_threshold & data$estimate_level<0 ,'Type']='Reduced above threshold'


myplot=ggplot(data, aes(x=estimate_level, y=logp,color=Type))+geom_point(alpha = 0.4 )
myplot
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlim(-2, 2) + 
scale_color_manual(values=c('grey','#039fbe','#cf1578'))+ 
theme(legend.text = element_text(size = 15))+
 theme(axis.title.y = element_text(size = 20))  + theme(axis.title.x = element_text(size = 20)) +
 theme(axis.text.y = element_text(size = 15)) + theme(axis.text.x = element_text(size = 15)) 

# COMMAND ----------

# library(tidyverse)

data=as.data.frame(severity_t1)
t1_p_threshold=-log10(0.05/7)
# data=as.data.frame(output)
# data=data[data$warnings_full=='None',]
data['logp']=-log10(data$p.value_level)
data['Type']='Below threshold'
data[data$logp>t1_p_threshold & data$mtb_id %in% t2 ,'Type']='T2 elevated above threshold'
data[data$mtb_id %in% t2,'Type']='T2 elevated'
t2_data=data[data$mtb_id %in% t2,]

# data[data$logp>t1_p_threshold & data$mtb_id %in% log_lmm ,'Type']='Lmm reduced above threshold'
# data[data$mtb_id %in% log_lmm,'Type']='lmm Reduced'
# t2_data=data[data$mtb_id %in% log_lmm,]



ggplot(data, aes(x=estimate_level, y=logp)) +
  theme_classic() +
  geom_point(
    mapping = aes(colour = Type),
    shape = 1,
#     size = 3,
#     stroke = 2,
    alpha = 5 / 6
    ) + xlim(-2, 2) +
scale_color_manual(values=c('grey','#cf1578'))+
  geom_point(t2_data, 
    mapping = aes(x=estimate_level, y=logp,fill = Type, colour = NA),
#     size = 2.88,
    alpha = 5 /6,
    shape = 21
  ) + geom_hline(yintercept=-log10(0.05/7))
# +
#   scale_fill_gradient(low = "blue", high = "yellow")

# COMMAND ----------

# library(tidyverse)

data=as.data.frame(severity_t1)
t1_p_threshold=-log10(0.05/7)
# data=as.data.frame(output)
# data=data[data$warnings_full=='None',]
data['logp']=-log10(data$p.value_level)
data['Type']='Below threshold'
data[data$logp>t1_p_threshold & data$mtb_id == 'rLC_neg_mtb_3755909' ,'Type']='peak above threshold'
data[data$mtb_id  == 'rLC_neg_mtb_3755909','Type']='Peak'
t2_data=data[data$mtb_id  == 'rLC_neg_mtb_3755909',]

ggplot(data, aes(x=estimate_level, y=logp)) +
  theme_classic() +
  geom_point(
    mapping = aes(colour = Type),
    shape = 1,
#     size = 3,
#     stroke = 2,
    alpha = 5 / 6
    ) + xlim(-2, 2) +
scale_color_manual(values=c('grey','#cf1578'))+
  geom_point(t2_data, 
    mapping = aes(x=estimate_level, y=logp,fill = Type, colour = NA),
#     size = 2.88,
    alpha = 5 /6,
    shape = 21
  ) 


# COMMAND ----------

print(length(data[data$logp>t1_p_threshold & data$estimate_level<0 ,'mtb_id']))
print(length(data[data$logp>t1_p_threshold & data$estimate_level>0 ,'mtb_id']))
reduced_t1=data[data$logp>t1_p_threshold & data$estimate_level<0 ,'mtb_id']
elevated_t1=data[data$logp>t1_p_threshold & data$estimate_level>0 ,'mtb_id']
t1=data[data$logp>t1_p_threshold  ,'mtb_id']
data[data$logp>t1_p_threshold  ,c('mtb_id','estimate_level','p.value_level')]

# COMMAND ----------

v1=data[data$mtb_id=='rLC_neg_mtb_3638997',]
v1

# COMMAND ----------

# MAGIC %md ## t1 box plot

# COMMAND ----------

# mtb='rLC_neg_mtb_3258421' # best down

mtb='rLC_neg_mtb_5118761' # best up
# mtb='rLC_neg_mtb_294302'

# mtb='rLC_neg_mtb_3630952'
# mtb='rLC_pos_mtb_867621' # estimate reduced outlier
# mtb='rLC_pos_mtb_1035331' # estimate increased outlier

test_data2=mdata[mtb_id %in% mtb] 
test_data2=as.data.frame(test_data2)
V1sit1 = test_data2[test_data2$Timepoint == 'T1',c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat_T1','Visit')]
colnames(V1sit1)=c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat','Visit')
V1sit2 = test_data2[test_data2$Timepoint == 'T2',c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat_T2','Visit')]
colnames(V1sit2)=c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat','Visit')
df1 = rbind(V1sit1,V1sit2)
V1sit1=V1sit1[V1sit1$Sarnat>=0,]
print(dim(V1sit1))
# V1sit1$Sarnat=as.factor(V1sit1$Sarnat)

# COMMAND ----------

V1sit1

# COMMAND ----------

polr(factor(Sarnat)~level,data=V1sit1,Hess=TRUE) 

# COMMAND ----------

myplot=ggplot(V1sit1, aes(x=factor(Sarnat), y=level, fill=Timepoint)) +
  geom_boxplot()
# Add dots
myplot + geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5,
                 position=position_dodge(1))+

scale_color_manual(values=c('#039fbe','#cf1578'))+ ylab(mtb)+ 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black")) +
# # # # #   xlim(-1000, 1000) + 

theme(legend.text = element_text(size = 20))+
 theme(axis.title.y = element_text(size = 25))  + theme(axis.title.x = element_text(size = 25)) +
 theme(axis.text.y = element_text(size = 20)) + theme(axis.text.x = element_text(size =20)) + theme(
plot.title = element_text(color="black", size=25, face="bold"))


# COMMAND ----------

# MAGIC %md # t2

# COMMAND ----------

t2_logP_th=-log10(2.405066e-05)

data=as.data.frame(severity_t2)
# data=as.data.frame(output)
data=data[data$warnings_full=='None',]
data['logp']=-log10(data$p.value_level)
data['Type']='Below threshold'
data[data$logp>t2_logP_th & data$estimate_level>0 ,'Type']='Elevated above threshold'
data[data$logp>t2_logP_th & data$estimate_level<0 ,'Type']='Reduced above threshold'


myplot=ggplot(data, aes(x=estimate_level, y=logp,color=Type))+geom_point(alpha = 0.4 )
myplot
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlim(-2, 2) + 
scale_color_manual(values=c('grey','#039fbe','#cf1578'))+ 
theme(legend.text = element_text(size = 15))+
 theme(axis.title.y = element_text(size = 20))  + theme(axis.title.x = element_text(size = 20)) +
 theme(axis.text.y = element_text(size = 15)) + theme(axis.text.x = element_text(size = 15)) 

# COMMAND ----------

# library(tidyverse)

data=as.data.frame(severity_t2)
t2_p_threshold=-log10(0.05/7)
# data=as.data.frame(output)
# data=data[data$warnings_full=='None',]
data['logp']=-log10(data$p.value_level)
data['Type']='Below threshold'
# data[data$logp>t2_p_threshold & data$mtb_id == 'rLC_neg_mtb_3755909' ,'Type']='peak above threshold'
# data[data$mtb_id  == 'rLC_neg_mtb_3755909','Type']='Peak'
# peak_data=data[data$mtb_id  == 'rLC_neg_mtb_3755909',]

data[data$logp>t2_p_threshold & data$mtb_id %in% log_lmm ,'Type']='Lmm elevated above threshold'
data[data$mtb_id %in% log_lmm,'Type']='lmm elevated'
t2_data=data[data$mtb_id %in% log_lmm,]


ggplot(data, aes(x=estimate_level, y=logp)) +
  theme_classic() +
  geom_point(
    mapping = aes(colour = Type),
    shape = 1,
#     size = 3,
#     stroke = 2,
    alpha = 5 / 6
    ) + xlim(-2, 2) +
scale_color_manual(values=c('grey','#cf1578'))+
  geom_point(t2_data, 
    mapping = aes(x=estimate_level, y=logp,fill = Type, colour = NA),
#     size = 2.88,
    alpha = 5 /6,
    shape = 21
  )  + geom_hline(yintercept=-log10(0.05/7))


# COMMAND ----------

t2_logP_th=-log10(2.405066e-05)
print(length(data[data$logp>t2_logP_th & data$estimate_level<0 ,'mtb_id']))
print(length(data[data$logp>t2_logP_th & data$estimate_level>0 ,'mtb_id']))
reduced_t2=data[data$logp>t2_logP_th & data$estimate_level<0 ,'mtb_id']
elevated_t2=data[data$logp>t2_logP_th & data$estimate_level>0 ,'mtb_id']
t2=data[data$logp>t2_logP_th ,'mtb_id']

data[data$logp>t2_logP_th  ,c('mtb_id','estimate_level','p.value_level')]
# data[order(data$p.value_level ),c('mtb_id','estimate_level','logp')]
t2

# COMMAND ----------

# MAGIC %md ## t2 box plot

# COMMAND ----------

# mtb='rLC_neg_mtb_3367429' # best down

# mtb='rLC_neg_mtb_1825916' # best up
mtb='rLC_neg_mtb_1846759'

#the 7 t2
mtb='rLC_neg_mtb_1602111'
mtb='rLC_neg_mtb_1640306'
mtb='rLC_neg_mtb_1823778'
mtb='rLC_neg_mtb_1825916'
mtb='rLC_neg_mtb_1826191'
mtb='rLC_neg_mtb_1846759'
mtb='rLC_neg_mtb_1871859'

test_data2=mdata[mtb_id %in% mtb] 
test_data2=as.data.frame(test_data2)
V1sit1 = test_data2[test_data2$Timepoint == 'T1',c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat_T1','Visit')]
colnames(V1sit1)=c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat','Visit')
V1sit2 = test_data2[test_data2$Timepoint == 'T2',c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat_T2','Visit')]
colnames(V1sit2)=c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat','Visit')
df1 = rbind(V1sit1,V1sit2)
V1sit2=V1sit2[V1sit2$Sarnat>=1,]
print(dim(V1sit2))
V1sit2$Sarnat=as.factor(V1sit2$Sarnat)

# COMMAND ----------

myplot=ggplot(V1sit2, aes(x=Sarnat, y=level, fill=Timepoint)) +
  geom_boxplot()
# Add dots
myplot + geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5,
                 position=position_dodge(1))+

scale_color_manual(values=c('#039fbe','#cf1578'))+ ylab(mtb)+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black")) +
# # # # #   xlim(-1000, 1000) + 

theme(legend.text = element_text(size = 20))+
 theme(axis.title.y = element_text(size = 25))  + theme(axis.title.x = element_text(size = 25)) +
 theme(axis.text.y = element_text(size = 20)) + theme(axis.text.x = element_text(size =20)) + theme(
plot.title = element_text(color="black", size=25, face="bold"))


# COMMAND ----------

# MAGIC %md # peak severity

# COMMAND ----------

peak_th=-log10(6.565949e-05)
# peak_th=3
data=as.data.frame(severity_peak)
# data=as.data.frame(output)
data=data[data$warnings_full=='None',]
data['logp']=-log10(data$p.value_level)
data['Type']='Below threshold'
data[data$logp>peak_th & data$estimate_level>0 ,'Type']='Elevated above threshold'
data[data$logp>peak_th & data$estimate_level<0 ,'Type']='Reduced above threshold'


myplot=ggplot(data, aes(x=estimate_level, y=logp,color=Type))+geom_point(alpha = 0.4 )
myplot
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlim(-3, 3) + 
scale_color_manual(values=c('grey','#039fbe','#cf1578'))+ 
theme(legend.text = element_text(size = 15))+
 theme(axis.title.y = element_text(size = 20))  + theme(axis.title.x = element_text(size = 20)) +
 theme(axis.text.y = element_text(size = 15)) + theme(axis.text.x = element_text(size = 15)) 

# COMMAND ----------

print(length(data[data$logp>peak_th & data$estimate_level<0 ,'mtb_id']))
print(length(data[data$logp>peak_th & data$estimate_level>0 ,'mtb_id']))
reduced_peak=data[data$logp>peak_th & data$estimate_level<0 ,'mtb_id']
elevated_peak=data[data$logp>peak_th & data$estimate_level>0 ,'mtb_id']
peak=data[data$logp>peak_th,'mtb_id']
data[data$logp>peak_th  ,c('mtb_id','estimate_level','p.value_level')]

# COMMAND ----------

df_worsk_peak=mdata[(mdata$WorstSarnat > 0)&(mdata$Timepoint=='T1') ,c("Timepoint", 'Sarnat_T1','Patient_ID','level','mtb_id','Visit','GA','WorstSarnat','WorstGradeTime')]
colnames(df_worsk_peak)=c("Timepoint",'Sarnat_T1', 'Patient_ID','level','mtb_id','Visit','GA','WorstSarnat','WorstGradeTime')
df_worsk_peak=as.data.frame(df_worsk_peak)

df_worsk_peak$worstime=0
df_worsk_peak[df_worsk_peak$WorstGradeTime=="<6H",'Worstime']=3
df_worsk_peak[df_worsk_peak$WorstGradeTime=="6H",'Worstime']=6
df_worsk_peak[df_worsk_peak$WorstGradeTime=="24H",'Worstime']=24
df_worsk_peak[df_worsk_peak$WorstGradeTime=="48H",'Worstime']=48
df_worsk_peak[df_worsk_peak$WorstGradeTime=="72H",'Worstime']=72
df_worsk_peak[df_worsk_peak$WorstGradeTime=="96H",'Worstime']=96
df_worsk_peak[df_worsk_peak$WorstGradeTime=="D6",'Worstime']=144
df_worsk_peak[df_worsk_peak$WorstGradeTime=="D7",'Worstime']=168
df_worsk_peak[df_worsk_peak$WorstGradeTime=="D8",'Worstime']=192
df_worsk_peak[df_worsk_peak$WorstGradeTime=="D10",'Worstime']=240

unique(setDT(df_worsk_peak), by = c("Patient_ID", "mtb_id"))
print(dim(df_worsk_peak))
df_worsk_peak$WorstSarnat=as.factor(df_worsk_peak$WorstSarnat)
df_worsk_peak=setDT(df_worsk_peak)

# COMMAND ----------

# mtb='rLC_neg_mtb_3886754' #elevated
mtb='rLC_neg_mtb_3755909' #reduced

peak_mtb=df_worsk_peak[mtb_id %in% mtb]

# COMMAND ----------

myplot=ggplot(peak_mtb, aes(x=WorstSarnat, y=level,color=Sarnat_T1))+geom_point(alpha =0.3,size = 5)
myplot
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
# #   xlim(-1000, 1000) + 
# scale_color_manual(values=c('#039fbe','#cf1578'))+
ggtitle(mtb)+
theme(legend.text = element_text(size = 20))+
 theme(axis.title.y = element_text(size = 25))  + theme(axis.title.x = element_text(size = 25)) +
 theme(axis.text.y = element_text(size = 20)) + theme(axis.text.x = element_text(size =20)) + theme(
plot.title = element_text(color="black", size=25, face="bold"))

# COMMAND ----------

# MAGIC %md ## peak box plot

# COMMAND ----------


mtb='rLC_neg_mtb_3755909' # best down

# mtb='rLC_neg_mtb_3886754' # best up
# mtb='rLC_neg_mtb_294302' #overlap with general
df_worsk_peak=mdata[(mdata$WorstSarnat > 0)&(mdata$Timepoint=='T1') ,c("Timepoint", 'Sarnat_T1','Patient_ID','level','mtb_id','Visit','GA','WorstSarnat','WorstGradeTime')]
colnames(df_worsk_peak)=c("Timepoint",'Sarnat_T1', 'Patient_ID','level','mtb_id','Visit','GA','WorstSarnat','WorstGradeTime')
df_worsk_peak=as.data.frame(df_worsk_peak)

df_worsk_peak$worstime=0
df_worsk_peak[df_worsk_peak$WorstGradeTime=="<6H",'Worstime']=3
df_worsk_peak[df_worsk_peak$WorstGradeTime=="6H",'Worstime']=6
df_worsk_peak[df_worsk_peak$WorstGradeTime=="24H",'Worstime']=24
df_worsk_peak[df_worsk_peak$WorstGradeTime=="48H",'Worstime']=48
df_worsk_peak[df_worsk_peak$WorstGradeTime=="72H",'Worstime']=72
df_worsk_peak[df_worsk_peak$WorstGradeTime=="96H",'Worstime']=96
df_worsk_peak[df_worsk_peak$WorstGradeTime=="D6",'Worstime']=144
df_worsk_peak[df_worsk_peak$WorstGradeTime=="D7",'Worstime']=168
df_worsk_peak[df_worsk_peak$WorstGradeTime=="D8",'Worstime']=192
df_worsk_peak[df_worsk_peak$WorstGradeTime=="D10",'Worstime']=240

unique(setDT(df_worsk_peak), by = c("Patient_ID", "mtb_id"))
print(dim(df_worsk_peak))
df_worsk_peak$WorstSarnat=as.factor(df_worsk_peak$WorstSarnat)
df_worsk_peak=setDT(df_worsk_peak)
peak_mtb=df_worsk_peak[mtb_id %in% mtb]

# COMMAND ----------

Peak_Sarnar_Time=as.factor(peak_mtb$Worstime)
Peak_Sarnar_Time=lapply(Peak_Sarnar_Time, paste0, " hr")
peak_mtb$Peak_Sarnar_Time=as.factor(unlist(Peak_Sarnar_Time))
peak_mtb=peak_mtb[order(Worstime,decreasing=TRUE),]

# COMMAND ----------

myplot=ggplot(peak_mtb, aes(x=WorstSarnat, y=level, fill=factor(Worstime))) +
#       stat_boxplot(geom = "errorbar", width=0.5, position = position_dodge(1)) +
      geom_boxplot(position = position_dodge(1), outlier.shape = NA)+
      geom_point(aes(fill =  factor(Worstime),  group = factor(Worstime)), 
#                  color="black", 
                 alpha  =0.5, size=3, 
                 position = position_jitterdodge(jitter.width = .1, dodge.width = 1)) +
#       scale_shape_manual (values = c(21,23) ) +
      guides(fill = guide_legend(override.aes = list(shape = NA) ) )

myplot+ggtitle(mtb)+ 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black")) +
# # # # #   xlim(-1000, 1000) + 

theme(legend.text = element_text(size = 20))+
 theme(axis.title.y = element_text(size = 25))  + theme(axis.title.x = element_text(size = 25)) +
 theme(axis.text.y = element_text(size = 20)) + theme(axis.text.x = element_text(size =20)) + theme(
plot.title = element_text(color="black", size=25, face="bold"))


# COMMAND ----------


myplot=ggplot(peak_mtb, aes(x=WorstSarnat, y=level, fill=factor(Worstime))) +
  geom_boxplot()
# Add dots
myplot + geom_dotplot(binaxis='y', 
                      stackdir='center', 
#                       stackgroups = FALSE,
#                       binwidth = 1, 
#                       binpositions = "all",
#                       fill = factor(peak_mtb$Worstime),
#                       binwidth=1/25,
                      dotsize = 0.5,
                 position=position_dodge(3))

# scale_color_manual(values=c('#039fbe','#cf1578'))+ 
ggtitle(mtb)+ ylab(mtb)+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black")) +
# # # # #   xlim(-1000, 1000) + 
# guides(fill = guide_legend(override.aes = list(shape = NA) ) )

theme(legend.text = element_text(size = 20))+
 theme(axis.title.y = element_text(size = 25))  + theme(axis.title.x = element_text(size = 25)) +
 theme(axis.text.y = element_text(size = 20)) + theme(axis.text.x = element_text(size =20)) + theme(
plot.title = element_text(color="black", size=25, face="bold"))
myplot

# COMMAND ----------


myplot=ggplot(peak_mtb, aes(x=WorstSarnat, y=level,fill='red')) +
  geom_boxplot()
# Add dots


myplot=myplot+ggtitle(mtb)+ 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black")) +
# # # # #   xlim(-1000, 1000) + 
# guides(fill = guide_legend(override.aes = list(shape = NA) ) )

theme(legend.text = element_text(size = 20))+
 theme(axis.title.y = element_text(size = 25))  + theme(axis.title.x = element_text(size = 25)) +
 theme(axis.text.y = element_text(size = 20)) + theme(axis.text.x = element_text(size =20)) + theme(
plot.title = element_text(color="black", size=25, face="bold"))
myplot

# COMMAND ----------

# MAGIC %md # thoms t1

# COMMAND ----------

th_t1_pth=-log10(1.762099e-05)
data=as.data.frame(tmp_severity_t1)
# data=as.data.frame(output)
data=data[data$warnings_full=='None',]
data['logp']=-log10(data$p.value_level)
data['Type']='Below threshold'
data[data$logp>th_t1_pth & data$estimate_level>0 ,'Type']='Elevated above threshold'
data[data$logp>th_t1_pth & data$estimate_level<0 ,'Type']='Reduced above threshold'


myplot=ggplot(data, aes(x=estimate_level, y=logp,color=Type))+geom_point(alpha = 0.4 )
myplot
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
#   xlim(-2, 2) + 
scale_color_manual(values=c('grey','#039fbe','#cf1578'))+ 
theme(legend.text = element_text(size = 15))+
 theme(axis.title.y = element_text(size = 20))  + theme(axis.title.x = element_text(size = 20)) +
 theme(axis.text.y = element_text(size = 15)) + theme(axis.text.x = element_text(size = 15)) 

# COMMAND ----------

print(length(data[data$logp>3 & data$estimate_level<0 ,'mtb_id']))
print(length(data[data$logp>3 & data$estimate_level>0 ,'mtb_id']))
reduced_tmp_t1=data[data$logp>3 & data$estimate_level<0 ,'mtb_id']
elevated_tmp_t1=data[data$logp>3 & data$estimate_level>0 ,'mtb_id']
tmp_t1=data[data$logp>3,'mtb_id']

# COMMAND ----------

# MAGIC %md # thoms t2

# COMMAND ----------

th_t1_pth
data=as.data.frame(tmp_severity_t2)
# data=as.data.frame(output)
data=data[data$warnings_full=='None',]
data['logp']=-log10(data$p.value_level)
data['Type']='Below threshold'
data[data$logp>th_t1_pth & data$estimate_level>0 ,'Type']='Elevated above threshold'
data[data$logp>th_t1_pth & data$estimate_level<0 ,'Type']='Reduced above threshold'


myplot=ggplot(data, aes(x=estimate_level, y=logp,color=Type))+geom_point(alpha = 0.4 )
myplot
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
#   xlim(-2, 2) + 
scale_color_manual(values=c('grey','#039fbe','#cf1578'))+ 
theme(legend.text = element_text(size = 15))+
 theme(axis.title.y = element_text(size = 20))  + theme(axis.title.x = element_text(size = 20)) +
 theme(axis.text.y = element_text(size = 15)) + theme(axis.text.x = element_text(size = 15)) 

# COMMAND ----------

print(length(data[data$logp>th_t1_pth & data$estimate_level<0 ,'mtb_id']))
print(length(data[data$logp>th_t1_pth & data$estimate_level>0 ,'mtb_id']))
reduced_tmp_t2=data[data$logp>th_t1_pth & data$estimate_level<0 ,'mtb_id']
elevated_tmp_t2=data[data$logp>th_t1_pth & data$estimate_level>0 ,'mtb_id']
tmp_t2=data[data$logp>3,'mtb_id']

# COMMAND ----------

# MAGIC %md # apgar

# COMMAND ----------


apgar_p1_th=-log10(2.540426e-06)
apgar_p5_th=apgar_p1_th
apgar_p10_th=apgar_p1_th

data=as.data.frame(Apgar_m1_severity_t1)
# data=as.data.frame(output)
data=data[data$warnings_full=='None',]
data['logp']=-log10(data$p.value_level)
data['Type']='Below threshold'
data[data$logp>apgar_p1_th & data$estimate_level>0 ,'Type']='Elevated above threshold'
data[data$logp>apgar_p1_th & data$estimate_level<0 ,'Type']='Reduced above threshold'


myplot=ggplot(data, aes(x=estimate_level, y=logp,color=Type))+geom_point(alpha = 0.4 )
myplot
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
#   xlim(-2, 2) + 
scale_color_manual(values=c('grey','#039fbe','#cf1578'))+ 
theme(legend.text = element_text(size = 15))+
 theme(axis.title.y = element_text(size = 20))  + theme(axis.title.x = element_text(size = 20)) +
 theme(axis.text.y = element_text(size = 15)) + theme(axis.text.x = element_text(size = 15)) 

# COMMAND ----------

apgar_p1_th=apgar_p1_th
apgar_p5_th=apgar_p1_th
apgar_p10_th=apgar_p1_th

print(length(data[data$logp>apgar_p1_th & data$estimate_level<0 ,'mtb_id']))
print(length(data[data$logp>apgar_p1_th & data$estimate_level>0 ,'mtb_id']))
reduced_apgarm1=data[data$logp>apgar_p1_th & data$estimate_level<0 ,'mtb_id']
elevated_apgarm1=data[data$logp>apgar_p1_th & data$estimate_level>0 ,'mtb_id']
apgarm1=data[data$logp>apgar_p1_th,'mtb_id']

# COMMAND ----------

apgar_p1_th=apgar_p1_th
apgar_p5_th=apgar_p1_th
apgar_p10_th=apgar_p1_th

data=as.data.frame(Apgar_m5_severity_t1)
# data=as.data.frame(output)
data=data[data$warnings_full=='None',]
data['logp']=-log10(data$p.value_level)
data['Type']='Below threshold'
data[data$logp>apgar_p5_th & data$estimate_level>0 ,'Type']='Elevated above threshold'
data[data$logp>apgar_p5_th & data$estimate_level<0 ,'Type']='Reduced above threshold'


myplot=ggplot(data, aes(x=estimate_level, y=logp,color=Type))+geom_point(alpha = 0.4 )
myplot
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
#   xlim(-2, 2) + 
scale_color_manual(values=c('grey','#039fbe','#cf1578'))+ 
theme(legend.text = element_text(size = 15))+
 theme(axis.title.y = element_text(size = 20))  + theme(axis.title.x = element_text(size = 20)) +
 theme(axis.text.y = element_text(size = 15)) + theme(axis.text.x = element_text(size = 15)) 



# COMMAND ----------

print(length(data[data$logp>apgar_p5_th & data$estimate_level<0 ,'mtb_id']))
print(length(data[data$logp>apgar_p5_th & data$estimate_level>0 ,'mtb_id']))
reduced_apgarm5=data[data$logp>apgar_p5_th & data$estimate_level<0 ,'mtb_id']
elevated_apgarm5=data[data$logp>apgar_p5_th & data$estimate_level>0 ,'mtb_id']
apgarm5=data[data$logp>apgar_p1_th,'mtb_id']

# COMMAND ----------

apgar_p10_th=apgar_p1_th
data=as.data.frame(Apgar_m10_severity_t1)
# data=as.data.frame(output)
data=data[data$warnings_full=='None',]
data['logp']=-log10(data$p.value_level)
data['Type']='Below threshold'
data[data$logp>apgar_p10_th & data$estimate_level>0 ,'Type']='Elevated above threshold'
data[data$logp>apgar_p10_th & data$estimate_level<0 ,'Type']='Reduced above threshold'


myplot=ggplot(data, aes(x=estimate_level, y=logp,color=Type))+geom_point(alpha = 0.4 )
myplot
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
#   xlim(-2, 2) + 
scale_color_manual(values=c('grey','#039fbe','#cf1578'))+ 
theme(legend.text = element_text(size = 15))+
 theme(axis.title.y = element_text(size = 20))  + theme(axis.title.x = element_text(size = 20)) +
 theme(axis.text.y = element_text(size = 15)) + theme(axis.text.x = element_text(size = 15)) 



# COMMAND ----------

print(length(data[data$logp>apgar_p10_th & data$estimate_level<0 ,'mtb_id']))
print(length(data[data$logp>apgar_p10_th & data$estimate_level>0 ,'mtb_id']))
reduced_apgarm10=data[data$logp>apgar_p10_th & data$estimate_level<0 ,'mtb_id']
elevated_apgarm10=data[data$logp>apgar_p10_th & data$estimate_level>0 ,'mtb_id']
apgarm10=data[data$logp>3,'mtb_id']

# COMMAND ----------

# MAGIC %md # venn diagram

# COMMAND ----------

t1t2_marker_list <- list(
  T1_markers = t1, 
  T2_markers = t2
  )

# elevated_t1t2_marker_list <- list(
#   elevated_T1_markers = elevated_t1, 
#   elevated_T2_markers = elevated_t2
#   )

# elevated_t1_marker_list <- list(
#   elevated_T1_markers = elevated_t1, 
#   elevated_Thompson_T1_markers = elevated_tmp_t1,
#   elevated_T2_markers = elevated_t2, 
#   elevated_Thompson_T2_markers = elevated_tmp_t2
#   )

# reduced_t1_marker_list <- list(
#   reduced_1_markers = reduced_t1, 
#   reduced_Thompson_T1_markers = reduced_tmp_t1,
#   reduced_T2_markers = reduced_t2, 
#   reduced_Thompson_T2_markers = reduced_tmp_t2
#   )

# reduced_t1t2_marker_list <- list(
#   reduced_T1_markers = reduced_t1, 
#   reduced_T2_markers = reduced_t2
#   )


# all_marker_list <- list(
#   General_markers = log_lmm, 
#   T1_markers = t1, 
#   T2_markers = t2,
#   peak_markers = peak
#   )


# reduced_all_marker_list <- list(
#   General_markers = reduced_log_lmm, 
#   T1_markers = reduced_t1, 
#   T2_markers = reduced_t2,
#   peak_markers = reduced_peak
#   )


# elevated_all_marker_list <- list(
#   General_markers = elevated_log_lmm, 
#   T1_markers = elevated_t1, 
#   T2_markers = elevated_t2,
#   peak_markers = elevated_peak
#   )


all_Sarnart_marker_list <- list(
  General_markers = log_lmm, 
  T1_markers = t1, 
  T2_markers = t2,
  markers = peak
  )

# COMMAND ----------

if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")

# COMMAND ----------

library(ggvenn)

# COMMAND ----------


ggvenn(
  all_marker_list, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
  )

# COMMAND ----------

ggvenn(
  elevated_all_marker_list, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
  )

# COMMAND ----------

Overlap_marker_list <- list(
 T2_markers = t2, 
 General_markers = lmm
  )

ggvenn(
  Overlap_marker_list, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
  )

# COMMAND ----------

elevated_tmp_t1

# COMMAND ----------

ggvenn(
  elevated_t1_marker_list, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
  )

# COMMAND ----------

ggvenn(
  reduced_t1_marker_list, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
  )

# COMMAND ----------

# MAGIC %md # intersect

# COMMAND ----------

# all_marker_list <- list(
#   General_markers = lmm, 
#   T1_markers = t1, 
#   T2_markers = t2,
#   peak_markers = peak
#   )

all_marker_list <- list(
  Reduced_General_markers = reduced_lmm, 
  Elevated_General_markers = elevated_lmm, 
  Good_Reduced_General_markers = reduced_log_lmm, 
  Good_Elevated_General_markers = elevated_log_lmm, 
  Reduced_T1_markers = reduced_t1,
  Elevated_T1_markers = elevated_t1,
  Reduced_T2_markers = reduced_t2,
  Elevated_T2_markers = elevated_t2,
  Reduced_peak_markers = reduced_peak,
  Elevated_peak_markers = elevated_peak
  )

# COMMAND ----------

Reduce(intersect,list(t2,lmm))

# COMMAND ----------

Reduce(intersect,list(elevated_log_lmm,t2))

# COMMAND ----------

log_lmm

# COMMAND ----------

mtb='rLC_neg_mtb_3755909' #peak
test_data2=mdata[mtb_id %in% mtb] 
test_data2=as.data.frame(test_data2)
V1sit1 = test_data2[test_data2$Timepoint == 'T1',c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat_T1','Visit')]
colnames(V1sit1)=c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat','Visit')
V1sit2 = test_data2[test_data2$Timepoint == 'T2',c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat_T2','Visit')]
colnames(V1sit2)=c("Timepoint", 'Patient_ID','level','mtb_id','Sarnat','Visit')
df1 = rbind(V1sit1,V1sit2)
V1sit1=V1sit1[V1sit1$Sarnat>=1,]
df1=df1[df1$Sarnat>=1,]
print(dim(df1))
df1=setDT(df1)

myplot=ggplot(V1sit1, aes(x=Sarnat, y=level))+geom_point(alpha =0.3,size = 5,color='#039fbe')
myplot
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
# #   xlim(-1000, 1000) + 
# scale_color_manual(values=c('#039fbe','#cf1578'))+
ggtitle(mtb)+
theme(legend.text = element_text(size = 20))+
 theme(axis.title.y = element_text(size = 25))  + theme(axis.title.x = element_text(size = 25)) +
 theme(axis.text.y = element_text(size = 20)) + theme(axis.text.x = element_text(size =20)) + theme(
plot.title = element_text(color="black", size=25, face="bold"))

# COMMAND ----------

# Pie Chart with Percentages
slices <- c(dim(V1sit1[V1sit1$Sarnat==1,])[1], dim(V1sit1[V1sit1$Sarnat==2,])[1], dim(V1sit1[V1sit1$Sarnat==3,])[1])
lbls <- c("Mild", "Moderate", "Severe")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, 
   main="Pie Chart of Visit 1")

# COMMAND ----------

slices <- c(dim(V1sit2[V1sit2$Sarnat==1,])[1], dim(V1sit2[V1sit2$Sarnat==2,])[1], dim(V1sit2[V1sit2$Sarnat==3,])[1])
lbls <- c("Mild", "Moderate", "Severe")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, 
   main="Pie Chart of Visit 2")

# COMMAND ----------

slices <- c(dim(df1[df1$Sarnat==1,])[1], dim(df1[df1$Sarnat==2,])[1], dim(df1[df1$Sarnat==3,])[1])
lbls <- c("Mild", "Moderate", "Severe")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, 
   main="Pie Chart of Visits")

# COMMAND ----------

df1

# COMMAND ----------

slices <- c(dim(df1[df1$Timepoint=='T1',])[1], dim(df1[df1$Timepoint=='T2',])[1])
lbls <- c("Visit1", "Visit2")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, 
   main="Pie Chart of Visits")

# COMMAND ----------

# MAGIC %md # upset

# COMMAND ----------

install.packages("UpSetR")
library(UpSetR)

# COMMAND ----------

all_marker_list <- list(
  Reduced_General_markers = reduced_lmm, 
  Elevated_General_markers = elevated_lmm, 
  Good_Reduced_General_markers = reduced_log_lmm, 
  Good_Elevated_General_markers = elevated_log_lmm, 
  Reduced_T1_markers = reduced_t1,
  Elevated_T1_markers = elevated_t1,
  Reduced_T2_markers = reduced_t2,
  Elevated_T2_markers = elevated_t2,
  Reduced_peak_markers = reduced_peak,
  Elevated_peak_markers = elevated_peak
  )

# COMMAND ----------

all_marker_list

# COMMAND ----------

upset(fromList(all_marker_list), order.by = "degree",nsets = 10, number.angles = 0, point.size = 3.5, line.size = 2, 
    mainbar.y.label = "Number of metabolites", sets.x.label = "Number of metabolites",nintersects=3)

# COMMAND ----------


severity_marker_list <- list(
  Reduced_General_markers = reduced_log_lmm, 
  Elevated_General_markers = elevated_log_lmm, 
  Reduced_T1_markers = reduced_t1,
  Elevated_T1_markers = elevated_t1,
  Reduced_T2_markers = reduced_t2,
  Elevated_T2_markers = elevated_t2,
#   Reduced_peak_markers = reduced_peak,
#   Elevated_peak_markers = elevated_peak,
  Elevated_Apgar_1_markers= elevated_apgarm1,
  Reduced_Apgar_1_markers = reduced_apgarm1, 
  Elevated_Apgar_5_markers= elevated_apgarm5,
  Reduced_Apgar_5_markers = reduced_apgarm5, 
  Elevated_Apgar_10_markers= elevated_apgarm10,
  Reduced_Apgar_10_markers = reduced_apgarm10, 
  Elevated_Thompson_t1_markers= elevated_tmp_t1,
  Reduced_Thompson_t1_markers = reduced_tmp_t1, 
  Elevated_Thompson_t2_markers= elevated_tmp_t2,
  Reduced_Thompson_t2_markers = reduced_tmp_t2
  )

# COMMAND ----------

upset(fromList(severity_marker_list), order.by = "degree",nsets = 16, number.angles = 0, point.size = 3.5, line.size = 2, 
    mainbar.y.label = "Number of metabolites", sets.x.label = "Number of metabolites",nintersects=4)

# COMMAND ----------


