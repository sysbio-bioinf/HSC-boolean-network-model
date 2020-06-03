
#Required library
library(BoolNet)
par(xpd = T)
par(mar = c(5,10,3,2))

#LoadNetwork
HSC <- loadNetwork("./HSC_RULES_FOR_SUBMISSION")

#Grouping of genes in the final table
groupingHSC<-list(class=rev(c("", "", "", "")), 
                  index= rev(list(rev(c("External_quiescence", "External_cycling")),
                                  rev(c("PI3K", "TSC1_2", "mTORC1", "FOXO3A", "ATM", "ROS", "Mitochondria", "Autophagy")),
                                  rev(c("RAS","ETS", "MEF", "GSK3b", "CTNNB1", "cMYC", "BMI1", "MDM2", "TP53", "CDKN1C", "CDKN1A", "CDKN1B", "GFI1", "RB", "E2F", "CCND1", "CCNE1", "S_phase")), 
                                  rev(c("AKT", "CDKN2D", "CDKN2A", "Pro_apoptotic_proteins", "Anti_apoptotic_proteins", "CYCS", "Apoptosis", "Senescence")))))


######################Attractors for unperturbed Network
AllInAttr <- getAttractors(HSC, method="sat.exhaustive")
plotAttractors(AllInAttr, mode = "table",  grouping = groupingHSC, reverse= T, onColor = "#6fc381", offColor = "steelblue4", drawLegend = FALSE, title = "Main phenotypes for HSC", allInOnePlot = T)


#####################Progressions towards entering of cell cycle: from LT to ST to cycling HSC
ProgressionFromLToSTHSC <- plotSequence (HSC, startState = c(1,1,0,1,0,1,1,0,0,0,0,0,0,1,0,0,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0), onColor = "#6fc381", offColor = "steelblue4", title = "From LT-HSC to activated ST-HSC" , drawLegend = "FALSE", grouping = groupingHSC, reverse= T)
ProgressionFromSTHSCtoCycling <-  plotSequence (HSC, grouping = groupingHSC, reverse = T, startState = c(0,1,1,0,1,0,0,1,1,0,1,0,1,0,1,1,1,1,0,1,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0), onColor = "#6fc381", offColor = "steelblue4", title = "From ST-HSC to cycling HSC" , drawLegend = "FALSE")
ProgressionFromLTHSCtoCycling <- plotSequence(HSC, startState = c(0,1,0,1,0,1,1,0,0,0,0,0,0,1,0,0,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0), onColor = "#6fc381", offColor = "steelblue4", title = "From LT-HSC to cycling HSC" , drawLegend = "FALSE", grouping = groupingHSC, reverse= T)


##################################simulations of K.O. models for Main Manuscript

#FOXO3A knockout
par(mar=c(3,12,2,2))
HSC_FOXO3A_LOSS <- fixGenes(HSC, "FOXO3A", 0)
AttrFOXO3A <- getAttractors(HSC_FOXO3A_LOSS, method="sat.exhaustive")
plotAttractors(AttrFOXO3A, grouping = groupingHSC, reverse = T, onColor = "#6fc381", offColor = "steelblue4", drawLegend = FALSE, title = "FOXO3A K.O.")

#ATM knockout

par(mar=c(3,12,2,2))
HSC_ATM_LOSS <- fixGenes(HSC, "ATM", 0)
AttrATM <- getAttractors(HSC_ATM_LOSS, method="sat.exhaustive")
plotAttractors(AttrATM, grouping = groupingHSC, reverse = T, onColor = "#6fc381", offColor = "steelblue4", drawLegend = FALSE, title = "ATM K.O.")

#BMI knockout

par(mar=c(3,12,2,2))
HSC_BMI1_LOSS <- fixGenes(HSC, "BMI1", 0)
AttrBMI1 <- getAttractors(HSC_BMI1_LOSS, method="sat.exhaustive")
plotAttractors(AttrBMI1, grouping = groupingHSC, reverse = T, onColor = "#6fc381", offColor = "steelblue4", drawLegend = FALSE, title = "BMI-1 K.O.")

#BMI AND P53 knockout

par(mar=c(2,10,3,1))
HSC_P53_LOSS <- fixGenes(HSC, "TP53", 0)
HSC_BMI_P53_LOSS <- fixGenes(HSC_P53_LOSS, "BMI1", 0 )
AttrP53BMI <- getAttractors(HSC_BMI_P53_LOSS, method="sat.exhaustive")
plotAttractors(AttrP53BMI, grouping = groupingHSC, reverse = T, onColor = "#6fc381", offColor = "steelblue4", drawLegend = FALSE, title = "Resque of BMI-1 phenotype by loss of TP53")


#TP53 knockout

par(mar=c(2,10,3,1))
HSC_P53_LOSS <- fixGenes(HSC, "TP53", 0)
AttrP53 <- getAttractors(HSC_P53_LOSS, method = "sat.exhaustive")
plotAttractors(AttrP53, grouping = groupingHSC, reverse = T, onColor = "#6fc381", offColor = "steelblue4", drawLegend = FALSE, title = "p53 K.O.")

#MEF knockout

par(mar=c(2,10,3,1))
HSC_MEF_LOSS <- fixGenes(HSC, "MEF", 0)
Attr_MEF_LOSS <- getAttractors(HSC_MEF_LOSS, method= "sat.exhaustive")
plotAttractors(Attr_MEF_LOSS, grouping = groupingHSC, reverse = T, onColor = "#6fc381", offColor = "steelblue4", drawLegend = FALSE, title = "MEF K.O.", allInOnePlot = T)



#########################Robustness analysis
#Hamming Distance
set.seed(10000)
Hamming<-testNetworkProperties(HSC, numRandomNets = 1000, 
                               testFunction = "testTransitionRobustness",
                               testFunctionParams = list(numSamples=1000),
                               alternative="less")


#########################Cellcyclecheck   

#plot attractor of the Faure cell cycle model
par(mar=c(2,10,3,1))
data("cellcycle")
attractorsCellCycle<- getAttractors(cellcycle)
plotAttractors(attractorsCellCycle, allInOnePlot = F, drawLegend = F, onColor = "#6fc381", offColor = "steelblue4")

#sequence to attractor by using the as startstrates the ones of our cycling HSC attractor
#nodes state: CyclinD=1, E2F=1, RB=0, P27=0, CyclinE=0)

plotSequence(cellcycle, startState = c(1,0,1,1,0,0,0,0,0,0), onColor = "#6fc381", offColor = "steelblue4", drawLegend = F)


########################### SINGLE-PETURBATION EXPERIMENTS ####################
net <- HSC
perturbation <- c(0,1)
names(perturbation) <- c("knock-out", "overexpression")

groupingHSC<-list(class=rev(c("", "", "", "")), 
                  index= rev(list(rev(c("External_quiescence", "External_cycling")),
                                  rev(c("PI3K", "TSC1_2", "mTORC1", "FOXO3A", "ATM","ROS","Mitochondria","Autophagy")),
                                  rev(c("RAS","ETS", "MEF", "GSK3b", "CTNNB1", "cMYC", "BMI1", "MDM2", "TP53", "CDKN1C", "CDKN1A", "CDKN1B", "GFI1", "RB", "E2F", "CCND1", "CCNE1", "S_phase")), 
                                  rev(c("AKT", "CDKN2D", "CDKN2A", "Pro_apoptotic_proteins", "Anti_apoptotic_proteins", "CYCS", "Apoptosis", "Senescence")))))

for(g in seq_along(net$genes))
  for(p in perturbation)
  {
    tmp_net <- fixGenes(net, fixIndices = c(g), values=c(p))
    tmp_attr <- getAttractors(tmp_net,method = "sat.exhaustive")
    par(xpd = T)
    par(mar = c(5,10,3,2))
    par()
    plotAttractors(tmp_attr, 
                   title = paste(names(perturbation)[p + 1], " of ", net$genes[g], sep =""), 
                   grouping = groupingHSC, 
                   drawLegend = F,
                   onColor = "#6fc381", offColor = "steelblue4")
  }


##########################EXTRA REVIEWER 2: dependency of TP53 from external signals and FOXO3A###########################
########################## RAS HYPERACTIVATION DOWNREGULTES P53 IN ABSENCE OF EXTERNAL CYCLING #############################
RAS_KI <- fixGenes(HSC, "RAS", 1)
RAS_KI_CYCLING_KO <- fixGenes(RAS_KI, "External_cycling" , 0 )
AttrRAS <- getAttractors(RAS_KI, method="sat.exhaustive")
AttrRAS_EXT <- getAttractors(RAS_KI_CYCLING_KO, method="sat.exhaustive")

#attractor for only RAS/PI3K consititutively active (here all populations are depicted)
plotAttractors(AttrRAS, grouping = groupingHSC, reverse = T, onColor = "#6fc381", offColor = "steelblue4", drawLegend = FALSE, title = "RAS K.I-")

#attractor for RAS/PI3K constitutively active and no presence at all of external cycling signals
plotAttractors(AttrRAS_EXT, grouping = groupingHSC, reverse = T, onColor = "#6fc381", offColor = "steelblue4", drawLegend = FALSE, title = "RAS K.I. and EXTERNAL K.O.")


############################# IN ABSENCE OF CYCLING SIGNALS CONST FOXO3A and ATM ACTIVATE TP53###############################
#absence of external cycling signals fixed
EXT_CYCLING_KO <- fixGenes(HSC, "External_cycling" , 0 )
#constitutive FOXO3A/ATM
FOXO3AKI <- fixGenes(EXT_CYCLING_KO, "FOXO3A", 1)
attrEXT_CYCLING <- getAttractors(EXT_CYCLING_KO, method="sat.exhaustive")
attrEXT_CYCLINg_FOXO3 <- getAttractors(FOXO3AKI, method="sat.exhaustive")

#attractor for only absence of external cycling (here all populations are depicted)
plotAttractors(attrEXT_CYCLING, grouping = groupingHSC, reverse = T, onColor = "#6fc381", offColor = "steelblue4", drawLegend = FALSE, title = "External cycling absence")
#acctractor for absence of external cycling signals and consitutive FOXO3A/ATM (results are similar dependening on TP53 activation)
plotAttractors(attrEXT_CYCLINg_FOXO3, grouping = groupingHSC, reverse = T, onColor = "#6fc381", offColor = "steelblue4", drawLegend = FALSE, title = "External cycling absence and constitutive FOXO3A")

