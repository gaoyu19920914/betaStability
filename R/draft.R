#   #### testable data ####
#   # vegan: varechem + varespec
#   # # vegan: dune.env + dune
#   # # vegan: BCI.env + BCI
#   # # vegan: mite.env + mite
#   # rioja::Ponds, SWAP # only 1 variable in SWAP: pH
#   # HSAUR3::birds, gardenflowers, watervoles
#   #
#
#   # clean /check the input data frame
#   varechem <- na.omit(varechem)
#
#   #### 1. betaDiv ####
#   bray.dist <- betaDiv(varespec)
#
#   #### 2. envDist ####
#   vare.envnorm <- BBmisc::normalize(varechem, method = "range", margin = 2)
#   # TODO: 1. remove BBmisc dependency and add different normalize methods
#   # vare.envnorm <- scale(varechem)
#
#
#   # hypothesis1 : all environmental variables are of same importance and can be
#   # normalized to a range of 0~1.
#   vare.envdist <- dist(vare.envnorm, method = "euclidean")
#
#   varedf <- cbind(bray.dist, vare.envdist)
#   colnames(varedf) <- c("bray.dist", "vare.envdist")
#
#   #### 3. relPred ####
#
#   ##### 3.8 predict community then calculate diversity ####
#   # MicroEcoTools
#   # specificity
#   # microbiomeSeq
#
