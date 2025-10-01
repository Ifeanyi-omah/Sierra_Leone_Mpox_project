#script for running epifilter and creating some kind of epi curve from sl mpox data
setwd("/Users/tbrockfi/projects/sl_mpox/")

library(readxl)
library(data.table)
library(EpiEstim)
library(caTools)
library(ggplot2)
library(seqinr)
library(patchwork)

files.sources = list.files(path = '/Users/tbrockfi/projects/important_scripts/EpiFilter-master/R files/main')
for (i in 1:length(files.sources)) {
  source(paste0(c('/Users/tbrockfi/projects/important_scripts/EpiFilter-master/R files/main/', files.sources[i]), collapse = ''))
}

#genomic data
complete_genome_set = seqinr::read.fasta("working_data/Complete+SLE+sequence+25sept18.fasta")
all_sle_genomes = complete_genome_set[grep("SLE", names(complete_genome_set))]
names(all_sle_genomes)[grep("890_S15", names(all_sle_genomes))] = "890_S15|WesternUrban|SLE|2025-05-13"
write.fasta(all_sle_genomes,
            names = names(all_sle_genomes),
            file.out="working_data/mpox_25sep22_iib_seqs.fasta")

delphy_mut = 7.45/length(getSequence(all_sle_genomes[[1]]))
sl_metadata = data.table(read.csv("raw_data/sl_regions_from_wiki.csv"))
sl_metadata$District = gsub(" ", "", sl_metadata$District)
sequencing_metadata = data.table(sample = sapply(names(all_sle_genomes), function(x){strsplit(x, "\\|")[[1]][1]}, USE.NAMES=FALSE),
                                 district = sapply(names(all_sle_genomes), function(x){strsplit(x, "\\|")[[1]][2]}, USE.NAMES=FALSE),
                                 date = sapply(names(all_sle_genomes), function(x){strsplit(x, "\\|")[[1]][4]}, USE.NAMES=FALSE))
sequencing_metadata$date = as.Date(sequencing_metadata$date, "%Y-%m-%d")
#want to clean this up 
sequencing_metadata[district == "koinadugu"]$district = "Koinadugu"
sequencing_metadata[district == "Portloko"]$district = "PortLoko"
sequencing_metadata[district == "Port-Loko"]$district = "PortLoko"
sequencing_metadata[district == "WesternAreaUrban"]$district = "WesternUrban"
sequencing_metadata[district == "Freetown"]$district = "WesternUrban"
sequencing_metadata[district == "Tongo"]$district = "Kenema"
sequencing_metadata = sl_metadata[,.(District, Province)][sequencing_metadata, on=.(District=district)]
seq_by_date = sequencing_metadata[,.N,date][order(date)]

#epi data
epi_data = data.table(read_excel("raw_data/mpox_25aug5_epi_data_from_allan.xlsx"))
colnames(epi_data) = as.character(unlist(epi_data[1,]))
epi_data = epi_data[-1,]
epi_data[, date := as.Date(paste0("2025-", Days), format = "%Y-%d/%b")]
epi_data = epi_data[-nrow(epi_data),]
epi_data$`Grand total` = as.numeric(epi_data$`Grand total`)
overall_epi_data = epi_data[,.(`Grand total`, date)]
colnames(overall_epi_data) = c("cases","date")

daily_case_and_seq_data = overall_epi_data[seq_by_date, on=.(date)]
daily_case_and_seq_data[, week_start := as.Date(cut(date, breaks = "week", start.on.monday = TRUE))]
weekly_case_and_seq_data = daily_case_and_seq_data[,.(cases = sum(cases), seqs = sum(N)),week_start]
weekly_case_and_seq_data$case_delta = weekly_case_and_seq_data$cases - weekly_case_and_seq_data$seqs
weekly_case_and_seq_data[weekly_case_and_seq_data$case_delta < 0 | is.na(weekly_case_and_seq_data$case_delta),]$case_delta = 0

all_dates = data.table(date = seq(min(epi_data$date), max(epi_data$date), by=1))
all_dates = epi_data[,.(`Grand total`, date)][all_dates, on=.(date)]
colnames(all_dates) = c("cases","date")
all_dates[is.na(cases)]$cases = 0

#epi curve and Rt estiamtes 
#code adapted from 
Iday = all_dates$cases
dates  = all_dates$date
nday = length(dates)
tday = 1:nday
wdist = dgamma(tday, shape = 22, scale = 0.5) #need a source for this
Lday = rep(0, nday) 
for(i in 2:nday){
  Lday[i] = sum(Iday[seq(i-1, 1, -1)]*wdist[1:(i-1)])    
}
Rmin = 0.01; Rmax = 10; eta = 0.1
m = 400; 
pR0 = (1/m)*rep(1, m)
Rgrid = seq(Rmin, Rmax, length.out = m)

Rfilt = epiFilter(Rgrid, m, eta, pR0, nday, Lday[tday], Iday[tday], 0.025)
Ifilt = recursPredict(Rgrid, Rfilt[[4]], Lday[tday], Rfilt[[3]], 0.025)

Rsmooth = epiSmoother(Rgrid, m, Rfilt[[4]], Rfilt[[5]], nday, Rfilt[[6]], 0.025)
Ismooth = recursPredict(Rgrid, Rsmooth[[4]], Lday[tday], Rsmooth[[3]], 0.025)

r_data = data.table(date = dates[-1],
                    r_est = Rfilt[[3]][2:nday],
                    r_lower = Rfilt[[2]][, 2:nday][1,],
                    r_upper = Rfilt[[2]][, 2:nday][2,])
r_data[, week_start := as.Date(cut(date, breaks = "week", start.on.monday = TRUE))]
weekly_r_data = r_data[,.(r_est = mean(r_est), r_lower = mean(r_lower), r_upper = mean(r_upper)),week_start]


dt_long <- data.table::melt(weekly_case_and_seq_data, id.vars = "week_start", measure.vars = c("seqs", "case_delta"),
                            variable.name = "type", value.name = "count")
dt_long[, type := factor(type, levels = c("case_delta","seqs"))]

combined_epi_r_data = dt_long[weekly_r_data, on=.(week_start)]
combined_epi_r_data[is.na(count)]$count = 0


combined_epi_r_data[, week_start := as.Date(week_start)]


bar_data <- combined_epi_r_data[!is.na(type)]
bar_data = bar_data[-c(1:3),]
# Prepare line + ribbon data: one row per week with r_est, r_lower, r_upper
r_line <- combined_epi_r_data[!is.na(r_est), .(
  r_est = unique(r_est),
  r_lower = unique(r_lower),
  r_upper = unique(r_upper)
), by = week_start]
r_line = r_line[-c(1:3, 30:31),]
# Compute a scaling factor for overlaying line on bar plot
scale_factor <- max(bar_data$count, na.rm = TRUE) / max(r_line$r_est, na.rm = TRUE)

epi_curve_plot = ggplot() +
  # Stacked bar chart (count by type)
  geom_bar(data = bar_data, 
           aes(x = week_start, y = count, fill = type), 
           stat = "identity") +
  
  # Shaded ribbon for r_est confidence interval
  geom_ribbon(data = r_line,
              aes(x = week_start, 
                  ymin = r_lower * scale_factor, 
                  ymax = r_upper * scale_factor), 
              fill = "black", alpha = 0.2) +
  
  # r_est line scaled to match bar height
  geom_line(data = r_line, 
            aes(x = week_start, y = r_est * scale_factor), 
            color = "black", size = 1.2) +
  
  geom_point(data = r_line, 
             aes(x = week_start, y = r_est * scale_factor), 
             color = "black", size = 2) +
  
  # Dashed line at r_est = 1 (scaled)
  geom_hline(yintercept = 1 * scale_factor, 
             linetype = "dashed", color = "black") +
  
  # Axis labels
  scale_y_continuous(
    name = "Cases",
    sec.axis = sec_axis(~ . / scale_factor, name = "Rt"),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = c("seqs" = "#9A031E", "case_delta" = "#F08A8B"), labels=c("Cases","Sequences"),name="") +
  labs(
    x = "",
    title = ""
  ) +
  scale_x_date(date_breaks = "1 week", date_labels = "%Y-%m-%d", expand = c(0, 0)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=25,colour = "black", family = "Helvetica"),
        axis.text.x = element_text(size=25,colour = "black", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=25,colour = "black"),
        legend.text=element_text(size=25),
        axis.line = element_line(),
        plot.margin = margin(1,1,1,1, "cm"),
        legend.position = c(0.8, 0.7)) + annotate("segment",
                                                  x = as.Date("2025-03-24"), xend = as.Date("2025-03-24"),     
                                                  y = 610, yend = 0, lty='dotted',  color="grey40",size = 1.2) + 
  annotate("segment",
           x = as.Date("2025-04-28"), xend = as.Date("2025-04-28"), 
           y = 610, yend = 0, lty='dotted',  color="grey40",size = 1.2) + 
  annotate("text", x = as.Date("2025-03-24"), y = 700, label = "Vaccination in W.A.U.\nand W.A.R.", size=25/2.845) + 
  annotate("text", x = as.Date("2025-04-28"), y = 700, label = "Nationwide\nVaccination", size = 25/2.845)

ggsave("/Users/tbrockfi/projects/sl_mpox/plots/epi_curve_plot.pdf",epi_curve_plot,
       height=9,width=18,dpi=1000)


#overdispersion analysis 
cluster_sizes_abbv = data.table(read.csv("/Users/tbrockfi/projects/sl_mpox/working_data/mpox_25sep25_cluster_distrib.csv"))[,-1]
cluster_sizes_abbv[, size_label := as.character(cluster_size)]
cluster_sizes_abbv = rbind(cluster_sizes_abbv, list(5, cluster_sizes_abbv[cluster_size>=5,sum(count),], "5+"))
cluster_sizes_abbv = cluster_sizes_abbv[cluster_size <= 5,]
cluster_sizes_abbv[, size_label := factor(size_label, levels = c("1", "2", "3", "4", "5+"))]
all_seq_r_k = data.table(read.csv("/Users/tbrockfi/projects/sl_mpox/working_data/mpox_25sep25_r_and_k_07.csv"))
half_seq_r_k = data.table(read.csv("/Users/tbrockfi/projects/sl_mpox/working_data/mpox_25sep25_r_and_k_038.csv"))
tenth_seq_r_k = data.table(read.csv("/Users/tbrockfi/projects/sl_mpox/working_data/mpox_25sep25_r_and_k_007.csv"))
r_k_data = data.table(sampling_rate = c(100,100,
                                        50,50,
                                        10,10),
                      param = c("R","k",
                                "R","k",
                                "R","k"),
                      lower = c(all_seq_r_k$lower_95[1], all_seq_r_k$lower_95[2],
                                half_seq_r_k$lower_95[1], half_seq_r_k$lower_95[2],
                                tenth_seq_r_k$lower_95[1], tenth_seq_r_k$lower_95[2]),
                      mean = c(all_seq_r_k$mle_estim[1], all_seq_r_k$mle_estim[2],
                               half_seq_r_k$mle_estim[1], half_seq_r_k$mle_estim[2],
                               tenth_seq_r_k$mle_estim[1], tenth_seq_r_k$mle_estim[2]),
                      upper = c(all_seq_r_k$upper_95[1], all_seq_r_k$upper_95[2],
                                half_seq_r_k$upper_95[1], half_seq_r_k$upper_95[2],
                                tenth_seq_r_k$upper_95[1], tenth_seq_r_k$upper_95[2]))
distro = ggplot(cluster_sizes_abbv, aes(x = size_label, y = count)) +
  geom_bar(stat = "identity", fill = "#1f77b4") +
  labs(x = "Cluster Size", y = "Observations") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=14,colour = "black", family = "Helvetica"),
        axis.text.x = element_text(size=14,colour = "black"),
        axis.text.y = element_text(size=14,colour = "black"),
        legend.text=element_text(size=14),
        axis.line = element_line(),
        plot.margin = margin(1,1,1,1, "cm"),
        legend.position = c(0.5, 1),
        legend.direction = "horizontal")

p_R <- ggplot(r_k_data[param == "R"], aes(x = factor(sampling_rate), y = mean)) +
  geom_point(color = "#1f77b4", size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "#1f77b4") +
  labs(x = "Sampling Rate",y = "Reproductive number") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=14,colour = "black", family = "Helvetica"),
        axis.text.x = element_text(size=14,colour = "black"),
        axis.text.y = element_text(size=14,colour = "black"),
        legend.text=element_text(size=14),
        axis.line = element_line(),
        plot.margin = margin(1,1,1,1, "cm"),
        legend.position = c(0.5, 1),
        legend.direction = "horizontal")

p_k <- ggplot(r_k_data[param == "k"], aes(x = factor(sampling_rate), y = mean)) +
  geom_point(color = "#ff7f0e", size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "#ff7f0e") +
  labs(x = "Sampling Rate",y = "Overdispersion (k)") +theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=14,colour = "black", family = "Helvetica"),
        axis.text.x = element_text(size=14,colour = "black"),
        axis.text.y = element_text(size=14,colour = "black"),
        legend.text=element_text(size=14),
        axis.line = element_line(),
        plot.margin = margin(1,1,1,1, "cm"),
        legend.position = c(0.5, 1),
        legend.direction = "horizontal")
distro + p_R + p_k

### juniper 
combine_results = function(list1, list2){
  all_names = union(names(list1), names(list2))
  master_list = setNames(vector("list", length(all_names)), all_names)
  for(a_name in all_names){
    master_list[[a_name]] = c(list1[[a_name]], list2[[a_name]])
  }
  return(master_list)
}

source("/Users/tbrockfi/projects/important_scripts/important_juniper_helpers_25apr29.R")
load("/Users/tbrockfi/projects/sl_mpox/working_data/tbrockfi-mpox-mpox_25aug25_30k_iter_juniper_downsample.RData")
juniper_res = downsampled_res
juniper_res_summary = juniper0::summarize(juniper_res, burnin=0)
juniper_r = c(quantile(juniper_res_summary$R, 0.025), mean(juniper_res_summary$R), quantile(juniper_res_summary$R, 0.975))
juniper_pi = c(quantile(juniper_res_summary$pi, 0.025), mean(juniper_res_summary$pi), quantile(juniper_res_summary$pi, 0.975))
outbreak_size = rev(length(juniper_res[[3]])/juniper_pi)
juniper_tmrca = juniper_res_summary$time_of_MRCA

load("/Users/tbrockfi/projects/sl_mpox/working_data/mpox_25sep22_15ka_iter_juniper_downsample.RData")
downsampled_15ka = downsampled_res
load("/Users/tbrockfi/projects/sl_mpox/working_data/mpox_25sep22_15kb_iter_juniper_downsample.RData")
downsampled_15kb = downsampled_res

juniper_15ka_summary = juniper0::summarize(downsampled_15ka, burnin=0)
juniper_15kb_summary = juniper0::summarize(downsampled_15kb, burnin=0)

juniper_r = c(juniper_15ka_summary$R, juniper_15kb_summary$R)
juniper_r_cri = c(quantile(juniper_r, 0.025), mean(juniper_r), quantile(juniper_r, 0.975))

juniper_pi = c(juniper_15ka_summary$pi, juniper_15kb_summary$pi)
juniper_pi_cri = c(quantile(juniper_pi, 0.025), mean(juniper_pi), quantile(juniper_pi, 0.975))*100
sampling_cri = c(quantile(juniper_pi, 0.025)*100, mean(juniper_pi)*100, quantile(juniper_pi, 0.975)*100)
outbreak_size_cri = length(all_sle_genomes)/(sampling_cri/100)

juniper_tmrca = data.table(tmrca = c(juniper_15ka_summary$time_of_MRCA, juniper_15kb_summary$time_of_MRCA))
write.csv(juniper_tmrca, "working_data/mpox_25sep24_juniper_tmrca.csv")

juniper_tmrca_cri = c(quantile(juniper_tmrca$tmrca, 0.025,type=1), mean(juniper_tmrca$tmrca), quantile(juniper_tmrca$tmrca, 0.975,type=1))

sampling_data = data.frame(sampling_rate = c(juniper_15kb_summary$pi, juniper_15ka_summary$pi))
sampling_data$outbreak_size = length(all_sle_genomes)/sampling_data$sampling_rate

trans_by_province_15ka = transmission_by_variable_distributions(downsampled_15ka,sequencing_metadata,"Province",name_var = "sample",drop="missing")
trans_by_province_15kb = transmission_by_variable_distributions(downsampled_15kb,sequencing_metadata,"Province",name_var = "sample",drop="missing")
trans_by_province = combine_results(trans_by_province_15ka[[2]], trans_by_province_15kb[[2]])
province_r_data = data.table(province = names(lapply(trans_by_province, mean)),
                             mean = unlist(lapply(trans_by_province, mean)),
                             lower = unlist(lapply(trans_by_province, function(x){quantile(x, 0.025)})),
                             upper = unlist(lapply(trans_by_province, function(x){quantile(x, 0.975)})))
province_r_data_full = rbindlist(lapply(names(trans_by_province), function(nm){data.table(province = nm, value = trans_by_province[[nm]])}))
province_r_data_full[, province := factor(province, levels = names(sort(tapply(value, province, mean))))]
counts_by_province = sequencing_metadata[,.N,Province]

epidemic_size = ggplot(sampling_data, aes(x = outbreak_size)) +
  geom_density(fill = "skyblue", color="#5581B0", alpha = 0.5) +
  labs(x = "Epidemic Size", y = "Density") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=25,colour = "black", family = "Helvetica"),
        axis.text.x = element_text(size=25,colour = "black", angle = 45, vjust = 0.5),
        axis.text.y = element_text(size=25,colour = "black"),
        legend.text=element_text(size=25),
        axis.line = element_line(),
        plot.margin = margin(1,1,1,1, "cm"),
        legend.position = c(0.5, 1),
        legend.direction = "horizontal") +
  scale_y_continuous(breaks=c(0,0.00005,0.0001,0.00015), 
                     labels=c(expression(0),expression(5%*%10^-5),expression(1%*%10^-4),expression(1.5%*%10^-4)), expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0, 0)) + 
  geom_vline(xintercept = median(sampling_data$outbreak_size), color = "black", linetype = "dashed", size = 1)


r_by_province = ggplot(province_r_data_full, aes(x = province, y = value, color = province)) +
  geom_boxplot(outlier.shape = NA) +  # remove default outlier points
  geom_jitter(width = 0.2, alpha = 0.3, size=2) +  # scatter points
  #scale_fill_manual(values = c("Eastern" = "#E15759", "Western" = "#9A031E","North West" = "#F08A8B",  "Southern" = "#F4C2C2",  "Northern" = "#7B002C" )) +
  scale_color_manual(values = c("Western" = "#CFA5A5", "Southern" = "#F08A8B",  "Eastern" = "#E15759","Northern" = "#9A031E","North West" = "#7B002C" ))+
  theme_bw() + 
  theme (panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         text = element_text(size=25,colour = "black", family = "Helvetica"),
         axis.text.x = element_text(size=25,colour = "black", angle = 45, vjust = 0.5), #, hjust = 1),
         axis.text.y = element_text(size=25,colour = "black"),
         legend.text=element_text(size=25),
         axis.line = element_line(),
         plot.margin = margin(1,1,1,1, "cm"),
         legend.position = "none",
         legend.direction = "horizontal") + 
  labs(x = "Province", y = "Reproduction Number") + 
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", size = 1)

juniper_plot = epidemic_size + r_by_province
ggsave("/Users/tbrockfi/projects/sl_mpox/plots/juniper_plot.pdf",juniper_plot,
       height=7,width=15,dpi=100)








