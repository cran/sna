######################################################################
#
# fileio.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 11/25/04
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains routines relating to file I/O.
#
# Contents:
#   read.nos
#
######################################################################


#read.nos - Read an input file in Neo-OrgStat format.  At this time, only the 
#graph stack is read; any coloring information is ignored.
read.nos<-function(file){
   #Get the formatting information
   f<-sapply(readLines(file,n=2),strsplit," ")
   #Parse the formatting information
   m<-as.numeric((f[[1]])[1])
   n<-as.numeric((f[[2]])[1])
   o<-as.numeric((f[[2]])[2])
   #Read the input data
   dat<-scan(file,skip=3)
   #Unpack the data in the proper order
   gstack<-array(dim=c(m,n,o))
   for(i in 1:m)
      for(j in 1:n)
         for(k in 1:o)
            gstack[i,j,k]<-dat[(i-1)*n*o+(j-1)*o+k]   
   #Return the stack
   gstack
}
