library(mendelianRandomization)






m.2011.reverse <- read.table("~/Bioinformatics/Src/R/reverse/final_tables/all_methods_cluster_reverse_2010.txt", header = T)
c.2011.reverse <- read.table("~/Bioinformatics/Src/R/reverse/final_tables/all_calculations_cluster_reverse_2010.txt", header = T)


m.2011.reverse$egger.p.fdr <- 1
m.2011.reverse[which(m.2011.reverse$nSNP >= 3), ]$egger.p.fdr <- p.adjust(p = m.2011.reverse[which(m.2011.reverse$nSNP >= 3), ]$egger.p, method = "fdr", n = length(m.2011.reverse[which(m.2011.reverse$nSNP >= 3), ]$egger.p))

m.2011.reverse$ivw.p.fdr <- 1
m.2011.reverse[which(m.2011.reverse$nSNP >= 3), ]$ivw.p.fdr <- p.adjust(p = m.2011.reverse[which(m.2011.reverse$nSNP >= 3), ]$ivw.p, method = "fdr", n = length(m.2011.reverse[which(m.2011.reverse$nSNP >= 3), ]$ivw.p))




for(phenotype in m.2011.reverse$phenotype){

  m <- m.2011.reverse[which(m.2011.reverse$phenotype == phenotype), ]
  c <- c.2011.reverse[which(c.2011.reverse$phenotype == phenotype), ]
  if(m$ivw.p.fdr < .05 | m$egger.p.fdr < .05){
    if(m$nSNP >= 3){
      p <- mr.plot(By = c$By, Bx = c$Bx, By.se = c$By.se, Bx.se = c$Bx.se, iv = c$iv, iv.se = c$iv.se, ivw = m$ivw, egger = m$egger,
                   egger.i = m$egger.i, chochran.Q = c$cochran.Q, ivw.Q = m$ivw.Q, egger.Q = m$egger.Q, egger.i.Q = m$egger.Q.i, egger.p.fdr = m$egger.p.fdr,
                   ivw.p.fdr = m$ivw.p.fdr, egger.i.p = m$egger.i.p, outcome.name = phenotype, exposure.name = "Celiac 2010")
      ggsave(filename = paste0("~/Bioinformatics/plots/exposure~celiac/2010/mr/", phenotype, ".png"), plot = ggdraw(p))
    }
  }
}


for(phenotype in m.2011.reverse$phenotype){

  m <- m.2011.reverse[which(m.2011.reverse$phenotype == phenotype), ]
  c <- c.2011.reverse[which(c.2011.reverse$phenotype == phenotype), ]
  if(m$ivw.p.fdr < .05 | m$egger.p.fdr < .05){
    if(m$nSNP >= 3){
      f <- mr.forest.plot(SNP = c$SNP, iv = c$iv, iv.se = c$iv.se, chochran.Q = c$cochran.Q, ivw = m$ivw, ivw.se = m$ivw.se, egger = m$egger, egger.se = m$egger.se, ivw.Q = m$ivw.Q,
                          ivw.se.Q = m$ivw.Q.se, egger.Q = m$egger.Q, egger.se.Q = m$egger.Q.se, outcome.name = phenotype, exposure.name = "Celiac 2010")
      ggsave(filename = paste0("~/Bioinformatics/plots/exposure~celiac/2010/forest/", phenotype, ".png"), plot = ggdraw(f))
    }
  }
}



for(phenotype in m.2011.reverse$phenotype){

  m <- m.2011.reverse[which(m.2011.reverse$phenotype == phenotype), ]
  c <- c.2011.reverse[which(c.2011.reverse$phenotype == phenotype), ]
  if(m$ivw.p.fdr < .05 | m$egger.p.fdr < .05){
    if(m$nSNP >= 3){
      if(length(which(is.na(c$beta.iv.se))) == 0){
        p.ivw <- mr.funnel.plot(iv = c$iv, iv.se = c$iv.se, method.estimate.before.Q = m$ivw, method.estimate.after.Q = m$ivw.Q, method.name = "IVW", linetype = "dashed",
                            chochran.Q = c$cochran.Q, outcome.name = phenotype, exposure.name = "Celiac 2011")
        ggsave(filename = paste0("~/Bioinformatics/plots/exposure~celiac/2011/funnel/ivw/", phenotype, ".png"), plot = p.ivw, width = 10, height = 8)

        p.ivw <- mr.funnel.plot(iv = c$iv, iv.se = c$iv.se, method.estimate.before.Q = m$egger, method.estimate.after.Q = m$egger.Q, method.name = "Egger", linetype = "solid",
                                chochran.Q = c$cochran.Q, outcome.name = phenotype, exposure.name = "Celiac 2011")
        ggsave(filename = paste0("~/Bioinformatics/plots/exposure~celiac/2011/funnel/egger/", phenotype, ".png"), plot = p.ivw, width = 10, height = 8)
      }
    }
  }
}

for(phenotype in m.2011.reverse$phenotype){

  m <- m.2011.reverse[which(m.2011.reverse$phenotype == phenotype), ]
  c <- c.2011.reverse[which(c.2011.reverse$phenotype == phenotype), ]
  if(m$beta.ivw.p.fdr < .05 | m$beta.egger.p.fdr < .05){
    if(m$nSNP >= 3){
      if(length(which(is.na(c$iv.se))) == 0){
        f <- mr.qq.plot(By = c$By, Bx = c$Bx, By.se = c$By.se, Bx.se = c$Bx.se, ivw = T, egger = T, chochran.Q = c$cochran.Q,
                        outcome.name = phenotype, exposure.name = "Celiac 2011")

        ggsave(filename = paste0("~/Bioinformatics/plots/2010/qq/", phenotype, ".png"), plot = f)


      }
    }
  }
}


celiac_2011 <- read.table("~/Bioinformatics/Src/R/cluster_results/outcome_without_HLL_region/outcome/celiac_removed/celiac_2011.txt", header = T, sep = " ")
celiac_2010 <- read.table("~/Bioinformatics/Src/R/cluster_results/outcome_without_HLL_region/outcome/celiac_removed/celiac_2010.txt", header = T, sep = " ")

install.packages('R.utils')
library(R.utils)

getAbsolutePath("exe/temp_clump")


find.package(package = "mendelianRandomization")




dir.exists(!paste0(system.file(package="mendelianRandomization", "executables"), "/temp_clump")){
  dir.create(file.path(system.file(package="mendelianRandomization", "executables"), "temp_clump"))
}



m.2011.reverse$egger.p.fdr <- 1
m.2011.reverse[which(m.2011.reverse$nSNP >= 3), ]$egger.p.fdr <- p.adjust(p = m.2011.reverse[which(m.2011.reverse$nSNP >= 3), ]$egger.p, method = "fdr", n = length(m.2011.reverse[which(m.2011.reverse$nSNP >= 3), ]$egger.p))

m.2011.reverse$ivw.p.fdr <- 1
m.2011.reverse[which(m.2011.reverse$nSNP >= 3), ]$ivw.p.fdr <- p.adjust(p = m.2011.reverse[which(m.2011.reverse$nSNP >= 3), ]$ivw.p, method = "fdr", n = length(m.2011.reverse[which(m.2011.reverse$nSNP >= 3), ]$ivw.p))






# m.2010 <- read.table("~/Bioinformatics/Src/R/cluster_results/outcome_without_HLL_region/final_tables/master.table.methods.2010.txt", header = T)
# c.2011 <- read.table("~/Bioinformatics/Src/R/cluster_results/outcome_without_HLL_region/final_tables/master.table.calculations.2011.txt", header = T)


m.2011.reverse <- read.table("~/Bioinformatics/Src/R/reverse/final_tables/all_methods_cluster_reverse_2010.txt", header = T)
c.2011.reverse <- read.table("~/Bioinformatics/Src/R/reverse/final_tables/all_calculations_cluster_reverse_2010.txt", header = T)

m.2011.reverse$egger.p.fdr <- 1
m.2011.reverse[which(m.2011.reverse$nSNP >= 3), ]$egger.p.fdr <- p.adjust(p = m.2011.reverse[which(m.2011.reverse$nSNP >= 3), ]$egger.p, method = "fdr", n = length(m.2011.reverse[which(m.2011.reverse$nSNP >= 3), ]$egger.p))

m.2011.reverse$ivw.p.fdr <- 1
m.2011.reverse[which(m.2011.reverse$nSNP >= 3), ]$ivw.p.fdr <- p.adjust(p = m.2011.reverse[which(m.2011.reverse$nSNP >= 3), ]$ivw.p, method = "fdr", n = length(m.2011.reverse[which(m.2011.reverse$nSNP >= 3), ]$ivw.p))

m.2011 <- m.2011[which(m.2011$beta.ivw.p.fdr < 0.05 | m.2011$beta.egger.p.fdr < 0.05), ]


int <- Reduce(intersect, list(m.2011$phenotype, m.2011.reverse$phenotype))




for(phenotype in int){
  m <- m.2011[which(m.2011$phenotype == phenotype), ]
  c <- c.2011[which(c.2011$phenotype == phenotype), ]
  p1 <- mr.plot(By = c$beta.outcome, Bx = c$beta.exposure, By.se = c$se.outcome, Bx.se = c$se.exposure, iv = c$beta.iv, iv.se = c$beta.iv.se, ivw = m$beta.ivw, egger = m$beta.egger, egger.i = m$beta.egger.i,
                chochran.Q = c$chochrans.q, ivw.Q = m$beta.ivw.Q, egger.Q = m$beta.egger.Q, egger.i.Q = m$beta.egger.i.Q, egger.p.fdr = m$beta.egger.p.fdr, ivw.p.fdr = m$beta.ivw.p.fdr,
                egger.i.p = m$beta.egger.i.p, outcome.name = "Celiac 2010", exposure.name = m$phenotype)

  m.1 <- m.2011.reverse[which(m.2011.reverse$phenotype == phenotype), ]
  c.1 <- c.2011.reverse[which(c.2011.reverse$phenotype == phenotype), ]
  p2 <- mr.plot(By = c.1$By, Bx = c.1$Bx, By.se = c.1$By.se, Bx.se = c.1$Bx.se, iv = c.1$iv, iv.se = c.1$iv.se, ivw = m.1$ivw, egger = m.1$egger, egger.i = m.1$egger.i, chochran.Q = c.1$cochran.Q,
                ivw.Q = m.1$ivw.Q, egger.Q = m.1$egger.Q, egger.i.Q = m.1$egger.Q.i, egger.p.fdr = m.1$egger.p.fdr, ivw.p.fdr = m.1$ivw.p.fdr, egger.i.p = m.1$egger.i.p,
                outcome.name = phenotype, exposure.name = "Celiac 2010")

  p <- plot_grid(ggdraw(p1), ggdraw(p2))
  ggsave(filename = paste0("~/Bioinformatics/plots/combined/2010/", phenotype, ".png"), plot = p, width = 17)
}





