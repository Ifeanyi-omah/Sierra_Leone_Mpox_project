#this version is for upload to vm 
setwd("/home/tbrockfi_broadinstitute_org/mpox/25aug25")
#setwd("/Users/tbrockfi/projects/sl_mpox/working_data")
library(juniper0)
library(data.table)
library(seqinr)
library(TAF)
library(coda)

#function for getting independent samples from juniper
independent_downsample = function(res, burnin_pct = 0.2){
  loglik = res[[1]][(1+length(res[[1]])*burnin_pct):length(res[[1]])]
  ess = floor(unname(coda::effectiveSize(loglik)))
  sample_times = floor(seq(1+length(res[[1]])*burnin_pct, length(res[[1]]), length.out=ess))
  my_res = list()
  my_res[[1]] = res[[1]][sample_times]
  my_res[[2]] = res[[2]][sample_times]
  my_res[[3]] = res[[3]]
  my_res[[4]] = res[[4]]
  my_res[[5]] = res[[5]]
  
  return(my_res)
}

all_seqs = seqinr::read.fasta("mpox_25sep22_iib_seqs.fasta",as.string=TRUE)
#create a juniper alignment and directory 
pos_list = c()
max_pos = nchar(all_seqs[[1]][1])
check_list = 1:(max_pos-1)
for(a_seq_num in 1:length(all_seqs)){ 
  for(current_pos in check_list){
    if(substring(all_seqs[[a_seq_num]][1], current_pos, current_pos+1) == "ga"){
      pos_list = c(pos_list,current_pos)
    }
    if(substring(all_seqs[[a_seq_num]][1], current_pos, current_pos+1) == "tc"){
      pos_list = c(pos_list,current_pos+1)
    }
  }
  pos_list = unique(pos_list)
  check_list = setdiff(check_list, pos_list)
  print(paste(a_seq_num, length(pos_list), sep=": "))
}
pos_list = sort(pos_list)

juniper_seq_list = vector("list", length = length(all_seqs))
for(a_seq_num in 1:length(all_seqs)){
  juniper_seq_list[[a_seq_num]] = getSequence(all_seqs[[a_seq_num]])[pos_list]
}
mkdir("juniper")
write.fasta(juniper_seq_list,
            names = names(all_seqs), 
            file.out="juniper/aligned.fasta")
juniper_metadata = data.table(sample = sapply(names(all_seqs), function(x){strsplit(x, "\\|")[[1]][1]}, USE.NAMES=FALSE),
                              date = sapply(names(all_seqs), function(x){strsplit(x, "\\|")[[1]][4]}, USE.NAMES=FALSE))
write.csv(juniper_metadata, "juniper/metadata.csv", row.names=FALSE)


#run juniper
juniper_mu = 6/length(pos_list)/365
init = initialize(a_g = 11.4, #generation interval of 11.4 days per Marziano et al. (2024)  
                  a_s = 15, #sojourn interval (time between infection date and sample collection date) of 15 days
                  init_mu = juniper_mu,
                  fixed_mu = TRUE, 
                  n_global=15000, 
                  ongoing = TRUE,
                  indir = "juniper")
res = run_mcmc(init, noisy = TRUE)

#analyze juniper results 
downsampled_res = independent_downsample(res)

#save results
save(downsampled_res, file="juniper/mpox_25sep22_15k_iter_juniper_downsample.RData")

