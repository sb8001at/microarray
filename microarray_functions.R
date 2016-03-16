directory_info <- "http://arayatakao.web.fc2.com/Rscripts/ATH1_info.txt"
source("http://arayatakao.web.fc2.com/Rscripts/mytheme.R") # load mytheme.R

###########################################################################################################
make_testdata <- function(){ # Making data for example analysis
name <- c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3", "D1", "D2", "D3")
probe <- read.table(directory_info, sep="\t", header=T)
probes <- probe[,1]
probes <- probes[1:50]

df <- NULL
df2 <- NULL
vec <- rep(1:4, rep(3,4))

for (i in 1: length(probes)){
	vec1 <- vec + 0.2 * i
	vec1 <- vec1 + rlnorm(length(vec1), 0.5, 0.5) + runif(length(vec1), 0, 0.5)
	df <- rbind(df, vec1)
	vec2 <- NULL
	for(k in 1:12){
		if (runif(1,0,1) > 0.2){
			vec2[k] <- "P"
		} else {
			vec2[k] <- "A"
		}
	}
	df2 <- rbind(df2, as.character(vec2))
}

df <- as.data.frame(df)
df2 <- as.data.frame(df2)
colnames(df) <- name
rownames(df) <- NULL
df <- cbind(probes, df)
colnames(df2) <- name
rownames(df2) <- probes

df3 <- data.frame(
	name,
	factor(c(1,1,1,2,2,2,3,3,3,4,4,4)),
	factor(c(1,2,3,1,2,3,1,2,3,1,2,3)),
	factor(c(1,2,1,2,1,2,1,2,1,2,1,2))
)

df4 <- data.frame(
	name,
	factor(c(0,0,0,0,0,0,0,0,0,1,1,1))
)

write.table(df, "Book1.txt", quote=F, col.names=T, row.names=F, sep="\t")
write.table(df2, "data_flaginfo.txt", quote=F, col.names=T, row.names=T, sep="\t")
write.table(df3, "factors.txt", quote=F, col.names=F, row.names=F, sep="\t")
write.table(df4, "templete.txt", quote=F, col.names=F, row.names=F, sep="\t")
}

###########################################################################################################
# loading libraries used for analysis
require(affy)
library(amap)
library(genefilter)
library(gplots)
library(ggplot2)
library(Heatplus)
library(MASS)
library(multtest)
library(mvtnorm)
library(RColorBrewer)
library(som)
library(stats)
library(tkWidgets)
library(dplyr)
library(tidyr)
library(data.table)
library(SimComp)
library(lawstat)
library(Biobase)
library(GEOquery)
library(stringr)

###########################################################################################################
# Making data.frame of results from soft.gz file
gse_df <- function(in_f, save=T){
  if(str_detect(in_f, "soft.gz")){
    temp <- getGEO(filename = in_f)  
    in_f2 <- str_split(in_f, ".soft")[[1]][1]
  }else{
    temp <- getGEO(in_f, GSEMatrix=F)
    in_f2 <- in_f
  }
  df <- NULL
  for(i in 1:length(names(GSMList(temp)))){
    df <- cbind(df, as.numeric(Table(GSMList(temp)[[i]])[,2]))
  }
  df <- log2(df)
  df <- as.data.frame(df)
  rownames(df) <- Table(GSMList(temp)[[i]])[,1]
  if(save==T){
    t <- strsplit(in_f, "_family.soft.gz")[[1]]
    t <- strsplit(t, ".soft.gz")[[1]]
    write.table(df, paste(t, ".txt", sep=""), col.names=T, row.names=T, quote=F, sep="\t")
  }else{
    return(df)
  }
}

###########################################################################################################

# read tab-separated text file as data.frame
read_m <- function(
	in_f, # input file
	col=T, #output file
	row=1 # set rownames in the text file
){
	d <- read.table(in_f, header=col, row.names=row, sep="\t")
	return(d)
}

###########################################################################################################
# Write variables as files
write_m <- function(
	df, # variables (data.frame)
	out_f, # output file
	row = T, # Save rownames in text
	col = T # Save colnames in text
){
	write.table(df, out_f, quote=F, sep="\t", col.names=col, row.names=row)
}

###########################################################################################################
# extract AGI codes
AGI_m <- function(
	in_f
){
	data <- read_m(in_f)
	df <- rownames(data)
	out_f <- paste(strsplit(in_f, ".txt"), "_agi_list.txt", sep="")
	write_m(df, out_f, row=F)
}

###########################################################################################################
# Correlation analysis with a gene or a vector

corr_m <- function(
	in_f, # input file
	pattern, # AGI code of a gene or the expression pattern
	ts = 0.7, # threshold level of correlation coefficient
	method = "spearman" # methods for calculating correlation coefficient
){
	data <- read.table(in_f, header=T, row.names=1, sep="\t")
  if(is.character(pattern)){
    AGI_data <- unlist(subset(data, row.names(data) == pattern)) 
  } else {
    AGI_data <- pattern
  }

  coefficient <- function(x){
    return(cor(AGI_data, unlist(x), method=method))
  }

  corr_coef <- apply(data,1,coefficient)
  result <- cbind(data, corr_coef)
  result_ts <- subset(result, abs(corr_coef)>ts)
  out_f <- paste(strsplit(in_f, ".txt"), "_", pattern,"_corr.txt", sep="")
	out_f2 <- paste(strsplit(in_f, ".txt"), "_", pattern,"_corr_ts.txt", sep="")
  write.table(result, out_f, sep="\t", col.names=T, row.names=T, quote=F)
	write.table(result_ts, out_f2, sep="\t", col.names=T, row.names=T, quote=F)
  result_ts <- result_ts[, -ncol(result_ts)]
  result_ts <- cbind(AGI=rownames(result_ts), result_ts)
  result_p <- gather(result_ts, sample, value, 2:ncol(result_ts))
  result_p <- cbind(result_p, agi_v=rep(AGI_data, rep(nrow(result_ts), length(AGI_data))))
  if(nrow(result_ts) < 50){
    out_f3 <- paste(strsplit(in_f, ".txt"), "_", pattern,"_corr_ts.pdf", sep="")
    p <- ggplot(result_p, aes(x=agi_v, y=value, color=AGI))+
      geom_point(size=5)+
      facet_wrap(~AGI, scales="free")+
      mytheme()
    ggsave(out_f3, p, width=nrow(result_ts)*1.1, height=nrow(result_ts), limitsize=F)
  }else{
    cat("Too many genes to make graph. Increase ts to decrease candidates.")
  }
}

##########################################################################################################
cutlow_m <- function(in_f, ts=1){
  data <- read.table(in_f, header=T, row.names=1, sep="\t")
  mvalue <- apply(data, 1, mean)
  data <- cbind(data, mvalue)
  data <- subset(data, data$mvalue>ts)
  data$mvalue <- NULL
  out_f <- paste(strsplit(in_f, ".txt"), "_cl",".txt", sep="")
  write.table(data, out_f, col.names=T, row.names=T, quote=F, sep="\t")
}

##########################################################################################################
twogene_plot <- 
  function(
    in_f,
    agi1,
    agi2,
    save=F
  ){
    data <- read.table(in_f, header=T, row.names=1, sep="\t")
    agi1_d <- data%>%subset(rownames(data)==agi1)%>%unlist()
    agi2_d <- data%>%subset(rownames(data)==agi2)%>%unlist()
    coef <- cor(agi1_d, agi2_d)
    coef_text <- paste("R=", sprintf("%1.3f",coef))
    df <- data.frame(agi1_d, agi2_d)
    p <- ggplot(df, aes(x=agi1_d, y=agi2_d))+
      geom_point(size=3)+
      labs(x=agi1, y=agi2)+
      annotate("text", x=max(agi1_d)-0.2, y=min(agi2_d), label=coef_text, color=1, size=5)+
      mytheme()
    if(save==T){
      out_f <- paste(strsplit(in_f, ".txt"), agi1,"_", agi2,".pdf", sep="")
      ggsave(out_f, p, width=5, height=5)  
    }else{
      p
    }
  }

###########################################################################################################

venplot_m <- function(
	in_f,
	lines
){
	data <- read.table(in_f, sep="\t", header=T, fill=T)
	result <- list()
	names <- colnames(data)
	setname <- names[lines]
	
	for(i in 1:length(lines)){
		agiset <- data[, lines[i]]
		agiset <- na.omit(unlist(agiset))
		agiset <- unique(agiset)
		result <- c(result, list(agiset))
	}

	names(result) <- setname
	venn(result)
}

###########################################################################################################

genefilter_m <- function(
	in_f,
	ts1 = 0.25,
	ts_level = 4.5,
	ts_IQR = 0.5
){
	data <- read.table(in_f, sep="\t", header=T, row.names=1)
	f1 <- pOverA(ts1, ts_level)
	f2 <- function(x)(IQR(x) > ts_IQR)
	ff <- filterfun(f1, f2)
	selected <- genefilter(data, ff)
	result <- data[selected,]
	out_f <- paste(strsplit(in_f, ".txt"), "_genefilter.txt", sep="")
	write.table(result, out_f, sep="\t")
}

###########################################################################################################

sammon_m <- function( # sammon
  in_f
){
  data <- read.table(in_f, sep="\t", header=T, row.names=1)
  data.d <- dist(data)
  sam <- sammon(data.d) 
  sam <- as.data.frame(sam$point)
  colnames(sam) <- c("component1", "component2")
  ggplot(sam, aes(x=component1, y=component2, label=rownames(sam)))+geom_text(size=3)+mytheme()
}

###########################################################################################################

relative_m <- function(
	in_f,
  cont
){
	data <- read.table(in_f, sep="\t", header=T, row.names=1)
	result <- 2^(data-apply(data[,cont], 1, mean))
	out_f <- paste(strsplit(in_f, ".txt"), "_r.txt", sep="")
	write.table(result, out_f, sep="\t", quote=F, col.names=T, row.names=T, append=F)
}

###########################################################################################################

panalysis_m <- function(
	in_f
){
	data <- read.table(in_f, sep="\t", header=T, row.names=1)
	data_acp <- acp(data)
	data_acp
	plot2(data_acp)
	plot(data_acp, variables=F, main="scores plot", cex=0.5)
}

###########################################################################################################

cmdscale_m <- function(
	in_f
){
	data <- read.table(in_f, sep="\t", header=T, row.names=1)
	data.d <- dist(data)
	loc <- cmdscale(data.d)
	x <- loc[,1]
	y <- loc[,2]
	plot(x, y, type="n", xlab="", ylab="", main="cmdscale_test")
	text(x, y, rownames(data), cex=0.5)
}

###########################################################################################################

twoway_m <- function(
  in_f,
  facA,
  facB,
  p_ts = 0.01,
  p_i_ts = 0.05,
  q_ts = 0.2
){
  data <- read.table(in_f, sep = "\t", header = T, row.names=1)
  facA <- factor(facA)
  facB <- factor(facB)
  twoway <- function(x){
    x <- unlist(x)
    two_way <- anova(aov(x~facA*facB))
    return(c(unlist(two_way)[13:15], unlist(two_way)[17:19]))
  }
  twoway_result <- as.data.frame(t(apply(data, 1, twoway)))
  two_way_names <- c("F_value_A", "F_value_B", "F_value_int", "p_value_A", "p_value_B", "p_value_int")
  colnames(twoway_result) <- two_way_names
  
  dh_results_A <- p.adjust(twoway_result$p_value_A, "BH")
  dh_results_B <- p.adjust(twoway_result$p_value_B, "BH")
  dh_results_int <- p.adjust(twoway_result$p_value_int, "BH")
  
  data <- cbind(data,twoway_result,BHR_A=dh_results_A, BHR_B=dh_results_B, BHR_int=dh_results_int)
  
  out_f <- paste(strsplit(in_f, ".txt"), "_twoway.txt", sep="")
  out_f2 <- paste(strsplit(in_f, ".txt"), "_twoway_facA_q", as.character(q_ts), ".txt",sep="")
  out_f3 <- paste(strsplit(in_f, ".txt"), "_twoway_facB_q", as.character(q_ts), ".txt",sep="")  
  write.table(data, out_f, sep="\t", quote=F, col.names=T, row.names=F, append=F)
  resultA <- subset(data, as.numeric(p_value_A) < p_ts & as.numeric(BHR_A) < q_ts & as.numeric(p_value_int) > p_i_ts) 
  resultB <- subset(data, as.numeric(p_value_B) < p_ts & as.numeric(BHR_B) < q_ts & as.numeric(p_value_int) > p_i_ts)
  
  write.table(resultA, out_f2, sep="\t", quote=F, col.names=T, row.names=F, append=F)
  write.table(resultB, out_f3, sep="\t", quote=F, col.names=T, row.names=F, append=F)
}

###########################################################################################################

threeway_m <- function(
  in_f,
  facA,
  facB,
  facC
){
  data <- read.table(in_f, sep = "\t", header = T, row.names=1)
  facA <- factor(facA)
  facB <- factor(facB)
  facC <- factor(facC)
  threeway <- function(x){
    x <- unlist(x)
    three_way <- anova(aov(x~facA*facB*facC))
    three_way <- unlist(three_way)[25:40]
    three_way <- three_way[-8]
    three_way <- three_way[-15]
    return(three_way)
  }
  threeway_result <- t(apply(data, 1, threeway))
  three_way_names <- c(
    "F_value_A", "F_value_B", "F_value_C", 
    "F_value_A_B", "F_value_A_C", "F_value_B_C", "F_value_A_B_C",
    "p_value_A", "p_value_B", "p_value_C", 
    "p_value_A_B", "p_value_A_C", "p_value_B_C", 
    "p_value_A_B_C"
  )
  threeway_result <- as.data.frame(threeway_result)
  colnames(threeway_result) <- three_way_names 
  dh_results_A <- p.adjust(threeway_result$p_value_A, "BH")
  dh_results_B <- p.adjust(threeway_result$p_value_B, "BH")
  dh_results_C <- p.adjust(threeway_result$p_value_C, "BH")
  dh_results_AB <- p.adjust(threeway_result$p_value_A_B, "BH")
  dh_results_AC <- p.adjust(threeway_result$p_value_A_C, "BH")
  dh_results_BC <- p.adjust(threeway_result$p_value_B_C, "BH")
  dh_results_ABC <- p.adjust(threeway_result$p_value_A_B_C, "BH")
  dh_results <- data.frame(
    BHR_A=dh_results_A,
    BHR_B=dh_results_B,
    BHR_C=dh_results_C,
    BHR_AB=dh_results_AB,
    BHR_AC=dh_results_AC,
    BHR_BC=dh_results_BC,
    BHR_ABC=dh_results_ABC)
  result <- cbind(data, threeway_result, dh_results)
  out_f <- paste(strsplit(in_f, ".txt"), "_threeway.txt", sep="")
  write.table(result, out_f, sep="\t", quote=F, col.names=T, row.names=T, append=F)
}

###########################################################################################################

kmeans_m <- function(
	in_f,
	par,
	seed = NULL
){
	set.seed(seed) 
	data <- read.table(in_f, sep="\t", header=T, row.names=1)
	data_kmeans <- kmeans(data, par)
	data_kmeans
	clust <- data_kmeans$cluster
	result <- cbind(data, cluster=clust)
	out_f <- paste(strsplit(in_f, ".txt"), "_kmeans.txt", sep="")
	write.table(result, out_f, sep="\t", append=F, quote=F, row.names=T)
}

###########################################################################################################

sort_m <- function(
  in_f,
  sortcol
){
  data <- read.table(in_f, sep="\t", header=T, row.names=1)
  data_s <- data%>%arrange(sortcol)
  data_s <- cbind(AGI=rownames(data), data_s)
  write.table(data_s, in_f, sep="\t", quote=F, col.names=T, row.names=F)
}

###########################################################################################################

ggplot_m <- function(
  in_f,
  cont,
  facA,
  facB=NULL,
  graphtitle = NULL,
  scaling = "free_y",
  relatives = T,
  line = F,
  savefile = F
){
  data <- read.table(in_f, sep="\t", header=T, row.names=1)
  facA <- factor(facA)
  data <- 2^data
  if (relatives == T){
    result <- data/apply(data[,cont], 1, mean)
  }else {
    result <- data
  }
  result <- na.omit(result)
  result <- as.data.frame(result)
  result <- cbind(AGI=rownames(result), result)
  df <- gather(result, exp, value, 2:ncol(result))
  
  if(!is.null(facB)){
    factorA <- rep(facA, rep(nrow(na.omit(data)), length(facA)))
    df <- cbind(df, factorA)
    df <- df%>%group_by(AGI, factorA)%>%summarise(m=mean(value), s=sd(value)/length(value)^0.5)
    p1 <- ggplot(df, aes(x=factorA, y=m, ymin=m-s, ymax=m+s, color=factorA, fill=factorA)) +
      geom_bar(stat="identity", position="dodge") +
      geom_linerange(position=position_dodge(width=0.9)) +
      facet_wrap(~AGI, scales=scaling) +
      theme(title=graphtitle) +
      mytheme()
    p2 <- 
      ggplot(df, aes(x=numeric(factorA), y=m, ymin=m-s, ymax=m+s, color=factorA, fill=factorA)) +
      geom_line(size=1) +
      geom_linerange(size=1) +
      geom_point(size=2.5) +
      facet_wrap(~AGI, scales=scaling) +
      mytheme() +
      expand_limits(y=0)
  } else {
    factorA <- rep(facA, rep(nrow(data), length(facA)))
    factorB <- rep(facB, rep(nrow(data), length(facB)))
    df <- cbind(df, factorA, factorB)
    df <- df%>%group_by(AGI, factorA, factorB)%>%summarize(m=mean(value), s=sd(value)/length(value)^0.5)
    p1 <- ggplot(df, aes(x=factorB, y=m, ymin=m-s, ymax=m+s, color=factorA, fill=factorA)) +
      geom_bar(stat="identity", position="dodge") +
      geom_linerange(position=position_dodge(width=0.9)) +
      facet_wrap(~AGI, scales=scaling) +
      theme(title=graphtitle) +
      mytheme()
    p2 <- 
      ggplot(df, aes(x=factorB, y=m, ymin=m-s, ymax=m+s, color=factorA, fill=factorA)) +
      geom_line(size=1) +
      geom_linerange(size=1) +
      geom_point(size=2.5) +
      facet_wrap(~AGI, scales=scaling) +
      mytheme() +
      expand_limits(y=0)
  }
  out_f <- paste(strsplit(in_f, ".txt"), ".pdf", sep="")
  if(line==F){
    if(savefile == T){
      ggsave(p1, file=out_f, width = nrow(data)/5, height = nrow(data)/6, limitsize=FALSE)
    }else{
      p1
    }
  }else{
    if(savefile == T){
      ggsave(p2, file=out_f, width = nrow(data)/3, height = nrow(data)/4, limitsize=FALSE)
    }else{
      p2
    }
  }
}

###########################################################################################################

coexist_m <- function(
  in_f_A,
  in_f_B,
  param=F
){
  data_A <- read.table(in_f_A, sep="\t")
  data_B <- read.table(in_f_B, sep="\t")
  vec_A <- rownames(data_A)
  vec_B <- rownames(data_B)
  
  if(param==T){
    lt <- list(vec_A, vec_B)
    venn(lt)
  }
  
  length_A <- length(vec_A)
  length_B <- length(vec_B)
  
  if (length_A > length_B){
    leng <- length_A
  }else{
    leng <- length_B
  }
  
  inters <- intersect(vec_A, vec_B)
  dif_A <- setdiff(vec_A, vec_B)
  dif_B <- setdiff(vec_B, vec_A) 
  leng_i <- leng - length(inters)
  leng_dA <- leng - length(dif_A)
  leng_dB <- leng - length(dif_B)
  inters <- c(inters, numeric(leng_i))
  dif_A <- c(dif_A, numeric(leng_dA))
  dif_B <- c(dif_B, numeric(leng_dB))
  df <- data.frame(Intersect=inters, Dif_1 <- dif_A, Dif_2 <- dif_B)
  out_f <- paste(strsplit(in_f_A, ".txt"), "_", strsplit(in_f_B, ".txt"),"_coexist.txt", sep="")
  write.table(df, out_f, sep="\t", quote=F, col.names=T, row.names=F, append=F)
}

###########################################################################################################

agi_data_m <- function(
  in_f, #list of AGI
  in_f2 #data
){
  if(length(grep(".txt", in_f))!=0){
    data <- read.table(in_f, sep="\t")
    agi <- as.character(unique(unlist(data[,1])))
  }else{
    agi <- in_f
  }
  agi <- agi[str_length(agi)>1]
  data2 <- fread(in_f2, sep="\t")
  setnames(data2, colnames(data2),c("AGI", scan(in_f2, nlines=1, sep="\t", what="character")))
  setkey(data2, AGI)
  
  result <- data2[agi]
  result <- na.omit(result)
  out_f <- paste(strsplit(in_f2, ".txt"), "_", strsplit(in_f, ".txt"), "_data.txt", sep="")
  write.table(result, out_f, sep="\t", append=F, quote=F, row.names=F)
}

###########################################################################################################

dunnett_m <- function(
  in_f,
  cont,
  samples
){
  data <- read.table(in_f, sep="\t", header=T, row.names=1)
  dunnett <- function(x){
    df <- data.frame(x,samples)
    vec <- unlist(SimTestDiff(df, grp="samples", resp="x", base=cont))
    dun <- vec[grep("p.val.adj", names(vec))]
    name <- vec[grep("comp.names", names(vec))]
    names(dun) <- name
    return(dun)
  }
  results <- as.data.frame(t(apply(data,1,dunnett)))
  results <- cbind(data, results)
  out_f <- paste(strsplit(in_f, ".txt"), "_dunnett.txt", sep="")
  write.table(results, out_f, col.names=T, row.names=T, quote=F, sep="\t")
}

###########################################################################################################

tukey_m <- function( # Tukey's test
  in_f, # Set input-file name
  samples # vector represent each treatment
){
  data<- read.table(in_f, sep = "\t", header = T, row.names=1)
  samples <- factor(samples)
  tukey <- function(x){
    x <- unlist(x)
    return(TukeyHSD(aov(x~samples))$samples[,4])
  }
  result <- as.data.frame(t(apply(data, 1, tukey)))
  result <- cbind(data, result)
  out_f <- paste(strsplit(in_f, ".txt"), "_tukey.txt", sep="")
  write.table(result, out_f, sep="\t", quote=F, col.names=T, row.names=T, append=F)
}

###########################################################################################################

hcluster_m <- function(
	in_f,
	fig_name = "cluster analysis",
	param = "complete" # ward, single, complete, average, mquitty, median, centroid 
){
	data <- read.table(in_f, header = T, row.names = 1, sep = "\t")
	data.d <- dist(data)
	data.hc <- hclust(data.d, method = param)
	plot(data.hc, hang = -1, main = fig_name)
}

###########################################################################################################

heatmap_m <- function(
	in_f,
	mainname = "",
	xlabel = "",
	ylabel = "",
	plot_fig = F,
	color_h = "by"
){
	data <- read.table(in_f, header=TRUE, row.names=1, sep="\t", quote="")
	redgreen <- colorRampPalette(c("#ff0000", "#000000", "#00ff00"))(100)
	blueyellow <- colorRampPalette(c("#0000ff", "#000000", "#ffff00"))(100)

	if (color_h == "by"){
		color_plot <- blueyellow
	} else if (color_h == "rg"){
		color_plot <- redgreen
	}
	out_f <- paste(strsplit(in_f, ".txt"), "_heatmap.pdf", sep="")
	postscript(out_f, horizontal=FALSE)
	heatplot <- heatmap(as.matrix(data), Rowv =NA, Colv=NA, scale="row", col=color_plot,  
	main=mainname, xlab=xlabel, ylab=ylabel, margin=c(8,10))
	dev.off()

	if (plot_fig == T){
		heatmap(as.matrix(data), Rowv =NA, Colv=NA, scale="row", col=color_plot,
		main=mainname, xlab=xlabel, ylab=ylabel, margin=c(8,10))
	}

}

###########################################################################################################

heatmap2_m <- function(
	in_f,
	group,
	labels,
	scaling ="row",
	plot_fig = F,
	color_h = "by"
){
	data <- read.table(in_f, sep = "\t", header = T, row.names=1)
	data <- subset(data, complete.cases(data))
	genes <- rownames(data)
	genes <- as.vector(genes)
	factors <- rep(1:length(group), group)
	factors <- as.factor(factors)

	df <- NULL

	for(i in 1:nrow(data)){
		row <- data[i,]
		row <- unlist(row)
		ave <- tapply(row, factors, mean)
		df <- rbind(df, ave)
	}

	rownames(df) <- genes
	colnames(df) <- labels

	redgreen <- colorRampPalette(c("#ff0000", "#000000", "#00ff00"))(100)
	blueyellow <- colorRampPalette(c("#0000ff", "#000000", "#ffff00"))(100)

	if (color_h == "by"){
		color_plot <- blueyellow
	} else if (color_h == "rg"){
		color_plot <- redgreen
	}
	out_f <- paste(strsplit(in_f, ".txt"), "_heatmap2.pdf", sep="")
	postscript(out_f, horizontal=FALSE)
	heatplot <- heatmap_2(as.matrix(df), scale = scaling, col=color_plot, legend=1, dist = function(c)Dist(c, method="spearman"))
	dev.off()

	if (plot_fig == T){
		heatplot <- heatmap_2(as.matrix(df), col=color_plot, legend=1, scale = scaling,
		dist = function(c)Dist(c, method="spearman"))
	}
}

###########################################################################################################

cdelete_m <- function(
	df,
	out_f,
	delete_lines
){
	data <- df[,-delete_lines]
	write.table(data, out_f, sep="\t", col.names=T, append=F, quote=F, row.names=T)
}

###########################################################################################################

cZscale_m <- function(
	in_f
){
	data <- read.table(in_f, header=TRUE, row.names=1, sep="\t", quote="")
	data.z <- normalize(data, byrow=FALSE)
	colnames(data.z) <- colnames(data)
	rownames(data.z) <- rownames(data)
	out_f <- paste(strsplit(in_f, ".txt"), "_cZscale.txt", sep="")
	write.table(data.z, out_f, sep="\t", col.names=T, append=F, quote=F, row.names=T)
}

###########################################################################################################

rZscale_m <- function(
	in_f
){
	data <- read.table(in_f, header=TRUE, row.names=1, sep="\t", quote="")
	data.z <- normalize(data, byrow=TRUE)
	colnames(data.z) <- colnames(data)
	rownames(data.z) <- rownames(data)
	out_f <- paste(strsplit(in_f, ".txt"), "_rZscale.txt", sep="")
	write.table(data.z, out_f, sep="\t", col.names=T, append=F, quote=F, row.names=T)
}

###########################################################################################################

genefinder_m <- function(
	in_f1,
	in_f2,
	out_f1,
	out_f2,
	param1 = 50,
	param2 = "zscore",
	param3 = "correlation" #"euclidean", "maximum", "manhattan", "canberra", "correlation", "binary"
){
	data <- read.table(in_f1, header=TRUE, row.names=1, sep="\t", quote="")
	data <- as.matrix(data)
	data_cl <- read.table(in_f2, sep="\t", quote="")
	template <- data_cl[,2] 
	tmp <- rbind(data, template)
	template_posi <- which(rownames(tmp) == "template")
	closeg <- genefinder(tmp, template_posi, param1, scale=param2, method=param3)
	topranked <- tmp[closeg[[1]]$indices,]
	tmp2 <- cbind(rownames(topranked), topranked)
	write.table(tmp2, out_f1, sep="\t", append=F, quote=F, row.names=F)
	write.table(rownames(topranked), out_f2, sep="\t", append=F, quote=F, row.names=F, col.names=F)
}

###########################################################################################################

oneway_m <- function( # conduct one-way ANOVA for microarray data
  in_f, # Set input-file name
  samples, # vector represent each treatment
  p_ts = 0.01, # threshold p-value of 1-way ANOVA
  q_ts = 0.2 # threshold q-value of Benjamin-Hochberg FDR
){
  data <- read.table(in_f, sep="\t", header=T, row.names=1)
  samples <- factor(samples)
  oneway <- function(x){
    x <- unlist(x)
    return(unlist(summary(aov(x~samples)))[9])
  }
  one_way_p <- apply(data, 1, oneway)
  dh_results <- p.adjust(one_way_p, "BH")  
  result <- cbind(data, one_way_p, dh_results)
  out_f <- paste(strsplit(in_f, ".txt"), "_oneway.txt", sep="")
  write.table(result, out_f, sep="\t", quote=F, col.names=T, row.names=T, append=F)
  
  result1 <- subset(result, as.numeric(one_way_p) < p_ts)
  result2 <- subset(result1, as.numeric(dh_results) < q_ts)
  out_f2 <- paste(strsplit(in_f, ".txt"), "_oneway_ts.txt", sep="") 
  write.table(result2, out_f2, sep="\t", quote=F, col.names=T, row.names=T, append=F)
}

###########################################################################################################

t_test_m <- function(
  in_f,
  group,
  var = "t_test",
  ts=1.5,
  xlimit = c(0,15),
  ylimit = c(0,15),
  xcolumn = "",
  ycolumn = "",
  pinXY = c(4,4),
  par1 = 1,
  par2 = 1,
  par3 = 0
){
  data <- read.table(in_f, sep = "\t", header = T, row.names=1)
  fact <- factor(rep(1:2, group))
  f_above <- function(x){return (par1 * par2 * x+log2(ts) + 1 - par1 + par3)} 
  f_below <- function(x){return ((2-par1) * par2 * x - log2(ts) + par1 - 1 + par3)}
  ttest <- function(x){
    if(var == "t_test"){
      t <- unlist(t.test(x[1:group[1]], x[(group[1]+1):sum(group)]))[3]
    } else if (var == "welch"){
      t <- unlist(t.test(x[1:group[1]], x[(group[1]+1):sum(group)], var.equal=F))[3]
    } else if (var == "wilcox"){
      t <- unlist(wilcox.test(x[1:group[1]], x[(group[1]+1):sum(group)]))[2]
    } else if (var == "brunner-munzel"){
      t <- unlist(brunner.munzel.test(x[1:group[1]], x[(group[1]+1):sum(group)]))[6]
    }
  }
  t_test_p <- apply(data, 1, ttest)
  result <- cbind(data, A_average=rowMeans(data[,1:group[1]]), B_average=rowMeans(data[,(group[1]+1):sum(group)]), 
                  Folds_change = 2^(rowMeans(data[,1:group[1]]) - rowMeans(data[,(group[1]+1):sum(group)])), 
                  p_value=t_test_p, q_value=p.adjust(t_test_p, "BH"))
  result1 <- subset(result, f_above(as.numeric(A_average)) < as.numeric(B_average))
  result2 <- subset(result, f_below(as.numeric(A_average)) > as.numeric(B_average))
  result_ts <- rbind(result1, result2)
  out_f <- paste(strsplit(in_f, ".txt"), "_t.test_", var, ".txt", sep="")
  out_f2 <- paste(strsplit(in_f, ".txt"), "_t.test_", var, "2.txt", sep="")
  write.table(result, out_f, sep="\t", quote=F, col.names=T, row.names=T, append=F)
  write.table(result_ts, out_f2, sep="\t", quote=F, col.names=T, row.names=T, append=F)
{
    par(pin = pinXY)
    plot(
      result$A_average, result$B_average, pch=20, 
      col= ifelse((result$B_average > f_above(result$A_average)|result$B_average < f_below(result$A_average)), 2, 5),
      xlab=xcolumn, ylab=ycolumn, xlim=xlimit, ylim=ylimit
    )
    par(new=T, pin = pinXY)
    plot(f_above, ann = F, col=4,  xlim=xlimit, ylim=ylimit)
    par(new=T, pin = pinXY)
    plot(f_below, ann = F, col=4,  xlim=xlimit, ylim=ylimit)
    par(pin = pinXY)
  }
}

##########################################################################################################

probe_to_AGI <- function( 
  in_f,
  in_f2= directory_info
){
  probe_agi <- fread(in_f2, header=TRUE)
  setkey(probe_agi, "Array Element")
  agi_v <- probe_agi[,2]
  unique_agi <- unique(unlist(agi_v))
  
  data <- read.table(in_f, header=TRUE, sep="\t", quote="", row.names=1)
  data_probe <- rownames(data)
  data_probe <- unlist(data_probe)
  
  probe_to_agi <- NULL
  for (i in 1:length(data_probe)){
    var1 <- probe_agi[data_probe[i]]
    probe_to_agi[i] <- var1[,2, with=F]
  }
  
  probe_to_agi <- factor(unlist(probe_to_agi))
  hoge <- cbind(AGI=probe_to_agi, data)
  hoge <- hoge %>%group_by(AGI)%>%summarise_each(funs_("mean"))
  out_f <- paste(strsplit(in_f, ".txt")[[1]], "_agi.txt", sep="")
  write.table(hoge, out_f, sep="\t", append=F, quote=F, col.names=T, row.names=F)
}

##########################################################################################################

flag_m <- function(
	in_f = "data_mas.txt",
	in_f2 = "data_flaginfo.txt",
	out_f = "data_mas_flag.txt",
	values,
	par_marginal = 0
){
	data <- read.table(in_f, sep="\t", header=T, row.names=1)
	flags <- read.table(in_f2, sep="\t", header=T, row.names=1)
	flag_abs_value <- c()
	
	for (i in 1:length(rownames(flags))){
		row <- flags[i,]
		row <- unlist(row)
		absent <- row[row == "A"]
		marginal <- row[row == "M"]
		absent_n <- length(absent)
		marginal_n <- length(marginal)
		abs_value <- absent_n + par_marginal * marginal_n
		flag_abs_value <- append(flag_abs_value, abs_value)
	}
	summary(flag_abs_value) 
	data <- cbind(data, flag=flag_abs_value)
	data <- subset(data, data$flag <= values)
	data$flag <- NULL
	write.table(data, out_f, sep="\t", quote=F, col.names=T, row.names=T, append=F)
}

##########################################################################################################

mas_norm <- function(
	out_f = "data_mas.txt",
	out_f2 = "data_flaginfo.txt",
	tkw=F,
	par=F,
	n.sample
){
	data <- ReadAffy(widget=tkw)
	eset <- mas5(data)

# flag information
	PAcalls <- mas5calls(data)
	PAcalls <- as.data.frame(PAcalls)
	PAcalls <- t(PAcalls)
	PAcalls <- PAcalls[-22811,]

	if (par == T) {
		vec=rep(ceiling(n.sample^0.5),2)
		par(mfrow=vec)
		image(data)
		win.graph()
		layout(matrix(c(1,2,3,3), 2, 2, byrow = T))
		hist(data)
		boxplot(data)
		plotAffyRNAdeg(AffyRNAdeg(data))
		summary(exprs(eset))
		win.graph()
		par(mfrow=vec)
		MAplot(data, cex=1.5)
	}
	exprs(eset)[exprs(eset) < 1] <- 1
	summary(exprs(eset))
	exprs(eset) <- log(exprs(eset), 2)
	id <- rownames(exprs(eset))
	PAcalls <- cbind(ID = id, PAcalls)
	write.exprs(eset, file=out_f)
	write.table(PAcalls, out_f2, sep = "\t", append=F, quote=F, row.names=F) 
}

##########################################################################################################

rma_norm <- function(
	out_f = "data_rma.txt",
	out_f2 = "data_flaginfo.txt",
	tkw = F,
	par = F,
	n.sample
){
	data <- ReadAffy(widget=tkw)
	eset <- rma(data)
	PAcalls <- mas5calls(data)
	PAcalls <- as.data.frame(PAcalls)
	PAcalls <- t(PAcalls)
	PAcalls <- PAcalls[-22811,] 
	if (par == T) {
		vec=rep(ceiling(n.sample^0.5),2)
		par(mfrow=vec)
		image(data)

		win.graph()
		layout(matrix(c(1,2,3,3), 2, 2, byrow = T))
		hist(data)
		boxplot(data)
		plotAffyRNAdeg(AffyRNAdeg(data))
		summary(exprs(eset))
		win.graph()
		par(mfrow=vec)
		MAplot(data, cex=1.5)
	}
	id <- rownames(exprs(eset))
	PAcalls <- cbind(ID = id, PAcalls)
	write.exprs(eset, file=out_f)
	write.table(PAcalls, out_f2, sep = "\t", append=F, quote=F, row.names=F)
}

###########################################################################################################