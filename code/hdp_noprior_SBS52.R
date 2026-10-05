options(stringsAsFactors = F)
mut_file<-commandArgs(T)[1]
iter <-as.numeric(commandArgs(T)[2])
out <-as.character(commandArgs(T)[3])


library(hdp)
print(paste0(out,iter,".rds"))
input_for_hdp<-read.table(mut_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)

input_for_hdp_sum <- apply(input_for_hdp,1,sum)
input_for_hdp <- input_for_hdp[apply(input_for_hdp,1,sum)>=1000,]
input_for_hdp_sum <- apply(input_for_hdp,1,sum)

median_muts = 3000
input_for_hdp = round(median_muts*input_for_hdp/input_for_hdp_sum)


ppindex <- c(0, rep(1, nrow(input_for_hdp)))
cpindex <- c(1, rep(2, nrow(input_for_hdp)))


hdp_mut <- hdp_init(ppindex = ppindex, # index of parental node
		    cpindex = cpindex, # index of the CP to use
		    hh = rep(1, 52), # prior is uniform over 52 categories
                    alphaa = rep(1, length(unique(cpindex))), # shape hyperparameters for 2 CPs
                    alphab = rep(5, length(unique(cpindex))))  # rate hyperparameters for 2 CPs

hdp_mut <- hdp_setdata(hdp_mut,  
		       dpindex = 2:numdp(hdp_mut), # index of nodes to add data to
		       input_for_hdp) # input data (mutation counts, sample rows match up with specified dpindex)

hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=10, seed=iter*200)

chlist <- hdp_posterior(hdp_activated,
                        burnin=20000,
			n=100,
			space=1000,
			cpiter=3,
			seed=iter*1e3)
saveRDS(chlist, paste0(out,iter,".rds"))
print(paste0("HDP chain saved to",out,iter,".rds"))

