  # load libs
  if(!require(ggplot2)){install.packages("ggplot2")}
  if(!require(gridExtra)){install.packages("gridExtra")}
  if(!require(stats)){install.packages("stats")}
  if(!require(car)){install.packages("car")}
  if(!require(plyr)){install.packages("plyr")}
  if(!require(scales)){install.packages("scales")}
  
  rm(list = ls())   # clean up working environment
  graphics.off()    # terminate graphics devices
  
  # define function to open and put into data frame
  get_data <- function(indicator_name){
  	# load all scenario & count rows and cols
  	filename <- paste("data4StatsTest_", indicator_name, ".csv", sep="")
  	ghge <- read.csv(filename, sep=";", dec=".", header=FALSE)
  	summary(ghge)
  	class(ghge)
  	class(ghge[,1])
  	
  	# determine data frame size
  	k <- 6 # No of scenarios  
  	n <- length(ghge[1,]) # No of MC realizations
  	N = k*n # Total DF size
  	
  	scenarios <- data.frame(scenario=rep(c("Ref", "Hp", "2ndSC", "2ndEU", "Cp", "All"), each=n), data=rep(as.double(NA),N))
  
  	scenarios$data[1:n] <- as.double(c(ghge[1,], recursive = TRUE))
  	scenarios$data[(n+1):(n*2)] <- as.double(c(ghge[2,], recursive = TRUE))
  	scenarios$data[(2*n+1):(n*3)] <- as.double(c(ghge[3,], recursive = TRUE))
  	scenarios$data[(3*n+1):(n*4)] <- as.double(c(ghge[4,], recursive = TRUE))
  	scenarios$data[(4*n+1):(n*5)] <- as.double(c(ghge[5,], recursive = TRUE))
  	scenarios$data[(5*n+1):(n*6)] <- as.double(c(ghge[6,], recursive = TRUE))
  		return(scenarios)
  }
  
  #############################################################
  # MAIN
  setwd("/home/rber/Work/BR_Stoch_PLUC_emissions/statistical_test")
  k <- 6 # 6scenarios, 6 groups to compare
  n <- 10000 # cheating
  N <- k*n
  indicator <- "TC" # select 'TC', 'SOC' or 'BC' to change the input files 
  
  # extra column with combi scen and year for letters
  tot <- data.frame(scenario=rep("", N), datas=rep(as.double(NA),N), combi=rep("", N), stringsAsFactors = FALSE)
  ghge <- get_data(indicator)
  
  tot$scenario[1:n] <- "1 Ref"
  tot$combi[1:n] <- "1_Ref" #as.character(scenarios_s_1990$combi[1:n])
  tot$datas[1:n] <- ghge$data[1:n]
  
  tot$scenario[(n+1):(2*n)] <- "2 HP"
  tot$combi[(n+1):(2*n)] <- "2_HP" #as.character(scenarios_s_2000$combi[1:n])
  tot$datas[(n+1):(2*n)] <- ghge$data[(n+1):(n*2)]
  
  tot$scenario[(2*n+1):(3*n)] <- "3 2ndSc"
  tot$combi[(2*n+1):(3*n)] <- "3_2ndSc" #as.character(scenarios_s_2009$combi[1:n])
  tot$datas[(2*n+1):(3*n)] <- ghge$data[(2*n+1):(n*3)]
  
  tot$scenario[(3*n+1):(n*4)] = "4 2ndEu"
  tot$combi[(3*n+1):(n*4)] = "4_2ndEu"
  tot$datas[(3*n+1):(n*4)] <- ghge$data[(3*n+1):(n*4)]
  
  tot$scenario[(4*n+1):(n*5)] = "5 CP"
  tot$combi[(4*n+1):(n*5)] = "5_CP"
  tot$datas[(4*n+1):(n*5)] <- ghge$data[(4*n+1):(n*5)]
  
  tot$scenario[(5*n+1) :(n*6)] = "6 All"
  tot$combi[(5*n+1):(n*6)] = "6_All"
  tot$datas[(5*n+1):(n*6)] <- ghge$data[(5*n+1):(n*6)]
  
  # plotting the data
  ggplot(tot, aes(x=scenario, y=datas)) + geom_boxplot()# + facet_wrap(~scenario)
  
  #############################################################
  # Significance tests
  
  if(!require(PMCMR)){install.packages("PMCMR")}
  if(!require(multcompView)){install.packages("multcompView")}
  if(!require(rcompanion)){install.packages("rcompanion")}
  if(!require(FSA)){install.packages("FSA")}
  
  # Kruskal-Wallis with post-hoc tests after Nemenyi
  kk <- kruskal.test(datas ~ as.factor(combi), data = tot) #from stats
  #kruskalTest(datas ~ scenario, data = tot, dist='Chisquare')
  # dist='Chisquare' gives strange(r) results
  pk <- posthoc.kruskal.nemenyi.test(datas ~ as.factor(combi), data = tot)
  kk
  pk
  pt <- pk$p.value
  #pt <- t(pt)

  # try Friedman instead of Kruskal-Wallis
  # does not work because the combination of 'combi' and 'scen' 
  # isn't unique because of MC analysis
  ##agg.data <- aggregate(datas ~ year+scen+combi, data = tot, median)
  ##f <- friedman.test(datas ~ year | scenario, data=agg.data)
  ##f
  ##pf <- posthoc.friedman.nemenyi.test(datas ~ scenario| year, data=agg.data)
  ##pt <- pf$p.value
  ##pt
  
  # https://stats.stackexchange.com/questions/292394/how-to-report-results-based-on-likert-item-group-differences/292540
  pt1 = fullPTable(pt)

  write.table(pt1, file = paste("pvalues_",indicator,".csv", sep=""), append = FALSE, quote = TRUE, sep = " ",
              eol = "\n", na = "NA", dec = ".", row.names = TRUE,
              col.names = TRUE) 
  
  letters <- multcompLetters(pt1,  
                  compare="<",  
                  threshold=0.01,
                  Letters=letters,  
                  reversed = FALSE)
  
  DF = as.data.frame(letters$Letters)
  DF$scenarios <- rownames(DF)
  DF
  
  tot$letter<-with(DF, letters$Letters[match(tot$combi, scenarios)])
  agg.data <- aggregate(datas ~ scenario+combi, data = tot, median)
  agg.data$letter<-with(DF, letters$Letters[match(agg.data$combi, scenarios)])
  agg.data
  
  #----------------------     PLOTS     ----------------------
  
  # five panels
  # graphics.off()    # terminate graphics devices
  # ggplot(tot, aes(x = scenario, y = datas, group=combi, fill=scenario)) +
  #   geom_boxplot(coef = 5) + facet_wrap(~scenario) +
  #   geom_text(data = agg.data, aes(label=letter, vjust=vjust, hjust=+0.5)) +
  #   ylab(expression(atop("GHG emissions", paste('gram CO'[2]*'-eq Mj'[EtOH])))) +
  #   xlab("")
  
  #vjust <- 0
  #hjust <-0.5
  #legpos <- "20"
  
  if(indicator == "TC") {col = "#B2C29D"}
  if(indicator == "SOC") {col = "#DFC284"}
  if(indicator == "BC") {col = "#7BC3BC"}
  
  ghge_min = min(ghge$data)
  ghge_max = max(ghge$data)
  # one panel
  graphics.off()    # terminate graphics devices
  ggplot(tot, aes(x=scenario, y = datas)) +
    stat_boxplot(coef=10, geom = "errorbar", colour = "#585858", size=0.2) + 
    geom_boxplot(coef=10, show.legend=FALSE, fill=col, colour = "#585858", size=0.35) + 
    geom_text(data = agg.data[agg.data$scenario=="1 Ref",], size=5, aes(label=letter, y=max(tot$datas[tot$scenario=="1 Ref"])-1,vjust=0)) +
    geom_text(data = agg.data[agg.data$scenario=="2 HP",], size=5, aes(label=letter, y=max(tot$datas[tot$scenario=="2 HP"])+1,vjust=0)) +
    geom_text(data = agg.data[agg.data$scenario=="3 2ndSc",], size=5, aes(label=letter, y=max(tot$datas[tot$scenario=="3 2ndSc"])+1,vjust=0)) +
    geom_text(data = agg.data[agg.data$scenario=="4 2ndEu",], size=5, aes(label=letter, y=max(tot$datas[tot$scenario=="4 2ndEu"])+1,vjust=0)) +
    geom_text(data = agg.data[agg.data$scenario=="5 CP",], size=5, aes(label=letter, y=max(tot$datas[tot$scenario=="5 CP"])+1,vjust=0)) +
    geom_text(data = agg.data[agg.data$scenario=="6 All",], size=5, aes(label=letter, y=max(tot$datas[tot$scenario=="6 All"])+1,vjust=0)) +
    xlab("") + 
    ylab(expression(atop("GHG emissions", paste('gram CO'[2]*'-eq Mj'[EtOH])))) +
    scale_y_continuous(breaks = seq(-300, 300, 30), limits=c(-30, 120)) +
    scale_x_discrete(labels=c("Ref.", "HP", 
                              expression(atop("2"^{nd}*"SC")), 
                              expression(atop("2"^{nd}*"EU")),
                              "CP", "All")) + 
    theme(plot.title=element_text(hjust = 0.5),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "grey", linetype = "solid", size=0.3, fill = NA),
          panel.grid.major = element_line(colour = "#E6E6E6", size=0.3),
          axis.text = element_text(colour="#585858", size = 13), 
          axis.title.y = element_text(colour="#585858", size = 12))

  
  
  #############################################################
  # OBS: I couldn't figure out how to use savePlot in linux. 
  # Therefore I was manually exporting the plots manually 
 
   # save figure 
  #savePlot(file = 'agb_letters.png', type = "png")
  #savePlot(file = 'agb_letters.pdf', type = "pdf")
  #savePlot(file = 'species_letters.png', type = "png")
  #savePlot(file = 'species_letters.pdf', type = "pdf")
  
  
