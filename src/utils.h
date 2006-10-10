/*
######################################################################
#
# utils.h
#
# copyright (c) 2006, Carter T. Butts <buttsc@uci.edu>
# Last Modified 8/24/06
# Licensed under the GNU General Public License version 2 (June, 1991)
# Portions taken from the NetStat library by Carter T. Butts (2002)
#  (self-licensed under GPL)
#
# Part of the R/sna package
#
# This file contains headers for utils.c.
#
######################################################################
*/
#ifndef UTILS_H
#define UTILS_H

/*DECLARATIONS/INCLUSIONS---------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>

/*The element datatype; contains a double value, an abstract pointer, and a
pointer to the next element. The purpose of the element is to serve in generic
stacks and queues, which are needed for a number of purposes (including the
BFS algorithms used here).*/

typedef struct elementtype{
   double val;
   void *dp;
   struct elementtype *next;
} element;

/*The snaNet datatype; contains incoming/outgoing edge lists for each vertex,
as well as network size information.  This is not intended to be a very fancy
structure (for that, use the network package), but is a useful, relatively
light-weight tool for backend processing of large, sparse graphs.*/

typedef struct snaNettype{
   int n, *outdeg,*indeg;
   element **oel,**iel;
} snaNet;


/*INTERNAL ROUTINES---------------------------------------------------------*/

/*snaNet ALLOCATION/MANIPULATION ROUTINES*/

snaNet *adjMatTosnaNet(double *mat, int *n);

int snaIsAdjacent(int i, int j, snaNet *g);


/*STACK/QUEUE/LIST ROUTINES*/

int isInList(element *head, double val);

element *listInsert(element *head, double val, void *dp);

element pop(element *head);

element *push(element *head, double val, void *dp);

element *clearstack(element *head);

long int stacklen(element *head);

char isinstack(element *head,double val);

element stackdel(element *head,double val);

element dequeue(element *head);

element *enqueue(element *head, double val, void *dp);

element *clearqueue(element *head);

long int queuelen(element *head);

char isinqueue(element *head,double val);

element queuedel(element *head,double val);


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

void aggarray3d_R(double *a, double *w, double *mat, int *m, int *n);

#endif
