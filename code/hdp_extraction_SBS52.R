library(hdp)
library(ggplot2)
library(RColorBrewer)
library(lsa)
library(lattice)
chlist_file<-commandArgs(T)[1]
hdp_input_path<-commandArgs(T)[2]
output_path<-commandArgs(T)[3]
prefix<-commandArgs(T)[4]
#------------extract HDP results---------
chlist <- vector("list", 10)
for (i in 1:10){
	  chlist[[i]] <- readRDS(paste0(chlist_file,i,".rds"))
}

mut_example_multi <- hdp_multi_chain(chlist)
pdf(paste0(output_path, "/",prefix,"_diagnostic_plots.pdf"), width = 7, height = 7)
par(mfrow = c(2, 2))
lapply(chains(mut_example_multi), plot_lik, bty="L", start = 1000)
lapply(chains(mut_example_multi), plot_numcluster, bty="L")
lapply(chains(mut_example_multi), plot_data_assigned, bty="L")
dev.off()

mut_example_multi <- hdp_extract_components(mut_example_multi)

input_for_hdp = read.table(hdp_input_path, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
channel_names <- colnames(input_for_hdp)
group_factor <- as.factor(gsub(".*\\[(.*)\\].*", "\\1", channel_names))
mut_colours <- c("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#F5ABCC")

pdf(paste0(output_path, "/",prefix,"_sig_plots.pdf"),width=8, height=4)
plot_comp_size(mut_example_multi, bty="L")
plot_comp_distn(mut_example_multi, cat_names=channel_names,
                grouping=group_factor, col=mut_colours,
                col_nonsig="grey80", show_group_labels=TRUE)
dev.off()

hdp_exposures=mut_example_multi@comp_dp_distn[["mean"]][2:dim(mut_example_multi@comp_dp_distn[["mean"]])[1],]
input_for_hdp <- input_for_hdp[apply(input_for_hdp,1,sum)>=1000,]
rownames(hdp_exposures)[(nrow(hdp_exposures)-nrow(input_for_hdp)+1):nrow(hdp_exposures)]=rownames(input_for_hdp)
write.csv(hdp_exposures,paste0(output_path, "/",prefix,"_HDP_exposure.csv"),quote = F)

hdp_sigs=data.frame(t(mut_example_multi@comp_categ_distn[["mean"]][1:dim(mut_example_multi@comp_categ_distn[["mean"]])[1],]))
rownames(hdp_sigs) <- channel_names
write.csv(hdp_sigs,paste0(output_path, "/",prefix,"_HDP_sigs.csv"),quote = F)

