plotAreaInROPE = function( mcmcChain , maxROPEradius , compVal=0.0 , 
                           HDImass=0.95 , ... ) {
  ropeRadVec = seq( 0 , maxROPEradius , length=201 ) # arbitrary comb
  areaInRope = rep( NA , length(ropeRadVec) )
  for ( rIdx in 1:length(ropeRadVec) ) {
    areaInRope[rIdx] = ( sum( mcmcChain > (compVal-ropeRadVec[rIdx]) 
                              & mcmcChain < (compVal+ropeRadVec[rIdx]) ) 
                         / length(mcmcChain) )
  }
  plot( ropeRadVec , areaInRope , 
        xlab=bquote("Radius of ROPE around "*.(compVal)) , 
        ylab="Posterior in ROPE" ,
        type="l" , lwd=4 , col="darkred" , cex.lab=1.5 , ... )
  # From http://www.indiana.edu/~kruschke/DoingBayesianDataAnalysis/Programs/ 
  # get HDIofMCMC.R. Put it in the working directory of this script, or specify
  # the path in the next line's source() command:
  # source("./HDIofMCMC.R") 
  HDIlim = HDIofMCMC( mcmcChain , credMass=HDImass )
  farHDIlim = HDIlim[which.max(abs(HDIlim-compVal))]
  ropeRadHDI = abs(compVal-farHDIlim)
  areaInFarHDIlim = ( sum( mcmcChain > (compVal-ropeRadHDI) 
                           & mcmcChain < (compVal+ropeRadHDI) ) 
                      / length(mcmcChain) )
  lines( c(ropeRadHDI,ropeRadHDI) , c(-0.5,areaInFarHDIlim) , 
         lty="dashed" , col="darkred" )
  text( ropeRadHDI , 0 , 
        bquote( atop( .(100*HDImass)*"% HDI limit" ,
                      "farthest from "*.(compVal) ) ) , adj=c(0.5,0) )
  lines( c(-0.5,ropeRadHDI) ,c(areaInFarHDIlim,areaInFarHDIlim) , 
         lty="dashed" , col="darkred" )
  text( 0 , areaInFarHDIlim , bquote(.(signif(areaInFarHDIlim,3))) , 
        adj=c(0,1.1) )
  return( list( ROPEradius=ropeRadVec , areaInROPE=areaInRope ) )
}