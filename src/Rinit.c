/*
######################################################################
#
# Rinit.c
#
# copyright (c) 2019, Carter T. Butts <buttsc@uci.edu>
# Last Modified 8/01/24
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later.
#
# Part of the R/sna package
#
# This file contains registration routines for R-callable C functions.
#
######################################################################
*/


#include <R.h>            
#include <Rinternals.h>            
#include <R_ext/Rdynload.h>

#include "cohesion.h"
#include "components.h"
#include "geodist.h"
#include "gli.h"
#include "layout.h"
#include "likelihood.h"
#include "nli.h"
#include "paths.h"
#include "randomgraph.h"
#include "triads.h"
#include "utils.h"

#define CALLDEF(name, n) {#name,(DL_FUNC) &name, n}

/*Call and C method definitions*/

static R_CallMethodDef CallEntries[] = {
  CALLDEF(bicomponents_R,3),                   /*cohesion.h*/
  CALLDEF(cliques_R,6),
  CALLDEF(reachability_R,3),                   /*components.h*/
  CALLDEF(geodist_R,6),                        /*geodist.h*/
  CALLDEF(geodist_val_R,6),
  CALLDEF(betweenness_R,9),                    /*nli.h*/
  CALLDEF(rgbern_R,5),                         /*randomgraph.h*/
  {NULL,NULL,0}
};

static R_CMethodDef CEntries[] = {
  CALLDEF(cutpointsDir_R,4),                   /*cohesion.h*/
  CALLDEF(cutpointsUndir_R,4),
  CALLDEF(kcores_R,7),
  CALLDEF(component_dist_R,3),                 /*components.h*/
  CALLDEF(compsizes_R,4),
  CALLDEF(undirComponents_R,4),
  CALLDEF(geodist_adj_R,4),                    /*geodist.h*/
  CALLDEF(maxflow_EK_R,5),
  CALLDEF(brokerage_R,5),                      /*gli.h*/
  CALLDEF(connectedness_R,4),
  CALLDEF(lubness_con_R,4),
  CALLDEF(gplot_layout_target_R,14),           /*layout.h*/
  CALLDEF(gplot_layout_fruchtermanreingold_R,15),
  CALLDEF(gplot_layout_fruchtermanreingold_old_R,10),
  CALLDEF(gplot_layout_kamadakawai_R,9),
  CALLDEF(gplot3d_layout_fruchtermanreingold_R,11),
  CALLDEF(gplot3d_layout_kamadakawai_R,10),
  CALLDEF(bn_dyadstats_R,3),                   /*likelihood.h*/
  CALLDEF(bn_triadstats_R,3),
  CALLDEF(bn_lpl_dyad_R,7),
  CALLDEF(bn_lpl_triad_R,8),
  CALLDEF(bn_ptriad_R,5),
  CALLDEF(degree_R,6),                         /*nli.h*/
  CALLDEF(evcent_R,8),
  CALLDEF(stresscent_R,5),
  CALLDEF(gilschmidt_R,5),
  CALLDEF(cycleCensus_R,9),                    /*paths.h*/
  CALLDEF(pathCensus_R,11),
  CALLDEF(bn_cftp_R,8),                        /*randomgraph.h*/
  CALLDEF(bn_mcmc_R,13),
  CALLDEF(udrewire_R,4),
  CALLDEF(wsrewire_R,5),
  CALLDEF(transitivity_R,6),                   /*triads.h*/
  CALLDEF(triad_census_R,6),
  CALLDEF(triad_classify_R,3),
  CALLDEF(aggarray3d_R,5),                     /*utils.h*/
  CALLDEF(dyadcode_R,4),
  CALLDEF(logadd_R,3),
  CALLDEF(logsub_R,4),
  {NULL,NULL,0}
};

void R_init_sna(DllInfo *dll)
{
   R_registerRoutines(dll,CEntries,CallEntries, NULL, NULL);
   R_useDynamicSymbols(dll, FALSE);
}

