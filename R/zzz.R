######################################################################
#
# zzz.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 11/13/04
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# .First.lib is run when the package is loaded with library(sna)
#
######################################################################

.First.lib <- function(lib, pkg){
    library.dynam("sna", pkg, lib)
}
