#####Fam and Rand vQTL FINAL#####
library("qtl")
library("vqtl")
#read in data
random <-read.cross(file = url("https://raw.githubusercontent.com/tbillman/Stapleton-Lab/master/vQTL%20Random%20and%20Family/data/tidied/Random2.csv"))
random <- drop.nullmarkers(random)
#scan with variance
random <- calc.genoprob(random)
routv <- scanonevar(cross = random,
                    mean.formula = height.in. ~ mean.QTL.add + mean.QTL.dom,
                    var.formula = ~ var.QTL.add + var.QTL.dom)
#####Set up our own function to extract effect sizes from mean_var_plot function#####
library("dplyr")
effect.sizes = function (cross, phenotype.name, focal.groups = NULL, nuisance.groups = NULL, 
                         genotype.names = c("AA", "AB", "BB"), xlim = NULL, ylim = NULL, 
                         title = paste(phenotype.name, "by", paste(focal.groups, 
                                                                   collapse = ", ")), draw_ribbons = TRUE, se_line_size = 1, 
                         point_size = 1) 
{
  indiv.mean.estim <- indiv.mean.lb <- indiv.mean.ub <- "fake_global_for_CRAN"
  indiv.sd.estim <- indiv.sd.lb <- indiv.sd.ub <- "fake_global_for_CRAN"
  group.mean.estim <- group.mean.ub <- group.mean.lb <- "fake_global_for_CRAN"
  group.sd.estim <- group.sd.ub <- group.sd.lb <- "fake_global_for_CRAN"
  modeling.df <- dplyr::data_frame(placeholder = rep(NA, qtl::nind(cross)))
  modeling.df[[phenotype.name]] <- cross[["pheno"]][[phenotype.name]]
  marker.names <- c(focal.groups[focal.groups %in% colnames(qtl::pull.geno(cross = cross))], 
                    nuisance.groups[nuisance.groups %in% colnames(qtl::pull.geno(cross = cross))])
  phen.names <- c(focal.groups[focal.groups %in% colnames(qtl::pull.pheno(cross = cross))], 
                  nuisance.groups[nuisance.groups %in% colnames(qtl::pull.pheno(cross = cross))])
  for (marker.name in marker.names) {
    modeling.df[[marker.name]] <- factor(x = qtl::pull.geno(cross = cross)[, 
                                                                           marker.name], labels = genotype.names)
  }
  for (phen.name in phen.names) {
    modeling.df[[phen.name]] <- factor(qtl::pull.pheno(cross = cross)[[phen.name]])
  }
  modeling.df[["placeholder"]] <- NULL
  covar.form <- paste(focal.groups, collapse = "+")
  if (!is.null(nuisance.groups)) {
    covar.form <- paste(covar.form, "+", paste(nuisance.groups, 
                                               collapse = "+"))
  }
  mean.form <- paste(phenotype.name, "~", covar.form)
  var.form <- paste("~", covar.form)
  dglm.fit <- dglm::dglm(formula = stats::formula(mean.form), 
                         dformula = stats::formula(var.form), data = modeling.df)
  mean.pred <- stats::predict(dglm.fit, se.fit = TRUE)
  mean.estim <- mean.pred$fit
  mean.se <- mean.pred$se.fit
  sd.pred <- stats::predict(dglm.fit$dispersion.fit, se.fit = TRUE)
  sd.estim <- sd.pred$fit/sd.pred$residual.scale
  sd.se <- sd.pred$se.fit
  indiv.prediction.tbl <- dplyr::bind_cols(stats::na.omit(modeling.df), 
                                           dplyr::data_frame(indiv.mean.estim = mean.estim, indiv.mean.lb = mean.estim - 
                                                               mean.se, indiv.mean.ub = mean.estim + mean.se, indiv.sd.estim = exp(sd.estim), 
                                                             indiv.sd.lb = exp(sd.estim - sd.se), indiv.sd.ub = exp(sd.estim + 
                                                                                                                      sd.se)))
  group.prediction.tbl <- indiv.prediction.tbl %>% dplyr::group_by_(.dots = c(focal.groups)) %>% 
    dplyr::summarise(group.mean.estim = mean(indiv.mean.estim), 
                     group.mean.lb = mean(indiv.mean.lb), group.mean.ub = mean(indiv.mean.ub), 
                     group.sd.estim = mean(indiv.sd.estim), group.sd.lb = mean(indiv.sd.lb), 
                     group.sd.ub = mean(indiv.sd.ub))
  return(group.prediction.tbl)
}
#set up a vector to run the function on
y = 1:length(routv$result$loc.name)
#effect sizes can not be computed for these 3 SNPs so we remove them from the vector
y = y[-c(458,2482,2483)]
#populating a dataframe with effect size estimates
rsizedf = sapply(y, function(x){
  tempm =  effect.sizes(cross = random,
                        phenotype.name = "height.in.",
                        genotype.names = c("AA","BB"),
                        focal.groups = routv$result$loc.name[x])
  tempv = c(tempm[1,2:7],tempm[2,2:7])
  return(unlist(tempv))
})
#gathering data from the initial scan
routvdf<- data.frame(routv$result$loc.name,
                     routv$result$pos,
                     routv$result$mean.lod,
                     routv$result$mean.asymp.p,
                     routv$result$var.lod,
                     routv$result$var.asymp.p,
                     routv$result$joint.lod,
                     routv$result$joint.asymp.p)
#dropping the SNPs whose effect sizes could not be computed
routvdf = routvdf[-c(458,2482,2483),]
#combining both 
routvdf = cbind(routvdf,t(rsizedf))
colnames(routvdf) = c("SNP Name",
                      "Position (cM)",
                      "Mean LOD",
                      "Mean P Value",
                      "Variance LOD",
                      "Variance P Value",
                      "Joint LOD",
                      "Joint P Value",
                      "A Mean Est",
                      "A Mean Lower Bound",
                      "A Mean Upper Bound",
                      "A Standard Deviation Est",
                      "A Standard Deviation Lower Bound",
                      "A Standard Deviation Upper Bound",
                      "B Mean Est",
                      "B Mean Lower Bound",
                      "B Mean Upper Bound",
                      "B Standard Deviation Est",
                      "B Standard Deviation Lower Bound",
                      "B Standard Deviation Upper Bound")
write.csv(routvdf, file = "RandomvQTL_LOD,Pvals,EffectSizes1010.csv")
#now with fam
#read in data
fam <-read.cross(file = url("https://raw.githubusercontent.com/tbillman/Stapleton-Lab/master/vQTL%20Random%20and%20Family/data/tidied/Family.csv"))
fam <- drop.nullmarkers(fam)
#scan with variance
fam <- calc.genoprob(fam)
foutv <- scanonevar(cross = fam,
                    mean.formula = PlantHeight ~ mean.QTL.add + mean.QTL.dom,
                    var.formula = ~ var.QTL.add + var.QTL.dom)
#set up a vector to run the function on
y = 1:length(foutv$result$loc.name)
#effect sizes can not be computed for these 3 SNPs so we remove them from the vector
y = y[-826]
#populating a dataframe with effect size estimates
fsizedf = sapply(y, function(x){
  print(x)
  tempm =  effect.sizes(cross = fam,
                        phenotype.name = "PlantHeight",
                        genotype.names = c("AA","BB"),
                        focal.groups = foutv$result$loc.name[x])
  tempv = c(tempm[1,2:7],tempm[2,2:7])
  return(unlist(tempv))
})
#gathering data from the initial scan
foutvdf<- data.frame(foutv$result$loc.name,
                     foutv$result$pos,
                     foutv$result$mean.lod,
                     foutv$result$mean.asymp.p,
                     foutv$result$var.lod,
                     foutv$result$var.asymp.p,
                     foutv$result$joint.lod,
                     foutv$result$joint.asymp.p)
#dropping the SNPs whose effect sizes could not be computed
foutvdf = foutvdf[-826,]
#combining both 
foutvdf = cbind(foutvdf,t(fsizedf))
colnames(foutvdf) = c("SNP Name",
                      "Position (cM)",
                      "Mean LOD",
                      "Mean P Value",
                      "Variance LOD",
                      "Variance P Value",
                      "Joint LOD",
                      "Joint P Value",
                      "A Mean Est",
                      "A Mean Lower Bound",
                      "A Mean Upper Bound",
                      "A Standard Deviation Est",
                      "A Standard Deviation Lower Bound",
                      "A Standard Deviation Upper Bound",
                      "B Mean Est",
                      "B Mean Lower Bound",
                      "B Mean Upper Bound",
                      "B Standard Deviation Est",
                      "B Standard Deviation Lower Bound",
                      "B Standard Deviation Upper Bound")
write.csv(foutvdf, file = "famvQTL_LOD,Pvals,EffectSizes1010.csv")
