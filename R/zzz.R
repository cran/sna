######################################################################
#
# zzz.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 10/10/06
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# .First.lib is run when the package is loaded with library(sna)
#
######################################################################

.First.lib <- function(lib, pkg){
    library.dynam("sna", pkg, lib)
    ehelp <- help(package="sna")$info[[1]]
    cat(paste(substring(ehelp[4],first=16),"\n",
              "Version ",substring(ehelp[2],first=16),
              " created on ",
               substring(ehelp[3],first=16),".\n", sep=""))
    cat(paste("copyright (c) 2005, Carter T. Butts, University of California-Irvine\n",sep=""))
    cat('Type help(package="sna") to get started.\n')
}
