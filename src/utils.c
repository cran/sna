/*
######################################################################
#
# utils.c
#
# copyright (c) 2006, Carter T. Butts <buttsc@uci.edu>
# Last Modified 8/24/06
# Licensed under the GNU General Public License version 2 (June, 1991)
# Portions taken from the NetStat library by Carter T. Butts (2002)
#  (self-licensed under GPL)
#
# Part of the R/sna package
#
# This file contains various utility routines to be called by other
# sna functions.
#
######################################################################
*/

#include "utils.h"

/*snaNet ALLOCATION/MANIPULATION ROUTINES---------------------------------*/

snaNet *adjMatTosnaNet(double *mat, int *n)
/*Given a vectorized adjacency matrix, create a snaNet object.  Note that any
non-zero/missing elements of mat are assumed to be valid adjacencies!  The
non-zero element values are stored within the snaNet object, in double form.*/
{
  snaNet *g;
  int i,j;
  double *dval;

  /*Allocate memory for the new object*/
  g=(snaNet *)R_alloc(1,sizeof(struct snaNettype));
  g->n=(*n);
  g->indeg=(int *)R_alloc(g->n,sizeof(int));
  g->outdeg=(int *)R_alloc(g->n,sizeof(int));
  g->iel=(element **)R_alloc(g->n,sizeof(element *));
  g->oel=(element **)R_alloc(g->n,sizeof(element *));

  /*Initialize the graph*/
  for(i=0;i<*n;i++){
    g->indeg[i]=0;
    g->outdeg[i]=0;
    g->iel[i]=NULL;
    g->oel[i]=NULL;
  }

  /*Add the edges*/
  for(i=0;i<*n;i++)
    for(j=0;j<*n;j++)
      if((!ISNAN(mat[i+j*(*n)]))&&(mat[i+j*(*n)]!=0.0)){
        dval=(double *)R_alloc(1,sizeof(double));   /*Create iel element*/
        dval[0]=mat[i+j*(*n)];
        g->iel[j]=listInsert(g->iel[j],(double)i,(void *)dval);
        g->indeg[j]++;
        dval=(double *)R_alloc(1,sizeof(double));   /*Create oel element*/
        dval[0]=mat[i+j*(*n)];
        g->oel[i]=listInsert(g->oel[i],(double)j,(void *)dval);
        g->outdeg[i]++;
      }

  /*Return the result*/
  return g;
}


int snaIsAdjacent(int i, int j, snaNet *g)
/*Determine whether i sends an edge to j.*/
{
  if(g->indeg[j]>g->outdeg[i])           /*Check shortest available list*/
    return isInList(g->oel[i],(double)j);
  else
    return isInList(g->iel[j],(double)i);
}


/*QUEUE/STACK/LIST ROUTINES-------------------------------------------------*/


int isInList(element *head, double val)
/*Is val in the sorted list pointed to by head?*/
{
  element *ep;
  
  for(ep=head;(ep!=NULL)&&(ep->val<val);ep=ep->next);
  if(ep==NULL)
    return 0;
  if(ep->val==val)
    return 1;
  return 0;
}


element *listInsert(element *head, double val, void *dp)
/*Add a new element to a sorted list, returning a pointer to the updated
list.*/
{
  element *elem,*ep;

  /*Initialize the element*/
  elem=(element *)R_alloc(1,sizeof(struct elementtype));
  elem->val=val;
  elem->dp=dp;
  elem->next=NULL;


  if(head==NULL){  /*If this is the only element, make it the head*/
    return elem;
  }else if(head->val>val){  /*If this is first, make it the head*/
    elem->next=head;
    return elem;
  }else{          /*Otherwise, traverse until we get to the right spot*/
    for(ep=head;(ep->next!=NULL)&&(ep->next->val<val);ep=ep->next);
    if(ep->next==NULL){   /*We ran out of list, apparently*/
      ep->next=elem;
      return head;
    }else{                /*We need to add elem after ep*/
      elem->next=ep->next;
      ep->next=elem;
      return head;
    }
  }
}


element pop(element *head)
/*Pop an element from the stack pointed to by head*/
{
element rval;

if(head==NULL){
   rval.val=-1.0;
   rval.dp=NULL;
   rval.next=head;
}else{
   if(head->next==NULL){
      rval.val=head->val;
      rval.dp=head->dp;
      head=NULL;
      rval.next=NULL;
   }else{
      rval.next=head->next;
      rval.val=head->val;
      rval.dp=head->dp;
      head=rval.next;
   }
}

return rval;
}


element *push(element *head, double val, void *dp)
/*Adds element with value val to the stack, returning the head 
pointer.*/
{
element *newnode;

/*Create the new entry*/
newnode=(element *)R_alloc(1,sizeof(struct elementtype));

newnode->val=val;   /*Set the element value*/
newnode->dp=dp;

/*Set the next pointer equal to the current first entry (if any)*/
newnode->next=head;

/*Place the new node at the head of the stack*/
head=newnode;

return head;
}


long int stacklen(element *head)
/*Returns the length of the stack pointed to by head*/
{
element *p;
int count=0;

for(p=head;p!=NULL;p=p->next)
   count++;

return count;
}


char isinstack(element *head,double val)
/*Returns a 1 if val is in the stack pointed to by head, otherwise 0*/
{
element *p;

for(p=head;p!=NULL;p=p->next)
   if(p->val==val)
      return 1;

return 0;
}


element stackdel(element *head,double val)
/*! Find the element val in the stack pointed to by head and delete it,
returning the deleted element.*/
{
element rval,*p;

if(head==NULL){
   rval.val=-1.0;
   rval.dp=NULL;
   rval.next=NULL;
}else if(head->val==val){
   rval.val=head->val;
   rval.dp=head->dp;
   rval.next=head->next;
   head=rval.next;
}else{
   for(p=head;(p->next!=NULL)&&(p->next->val!=val);p=p->next);
   if(p->next==NULL){
      rval.val=-1.0;
      rval.dp=NULL;
      rval.next=NULL;
   }else{
      rval.val=p->next->val;
      rval.dp=p->next->dp;
      rval.next=p->next->next;
      p->next=rval.next;
   }
}
      
return rval;
}


element dequeue(element *head)
/*Dequeue an element from the queue pointed to by head*/
{
element rval,*p;

if(head==NULL){
   rval.val=-1.0;
   rval.dp=NULL;
   rval.next=head;
}else{
   if(head->next==NULL){
      rval.val=head->val;
      rval.dp=head->dp; 
      head=NULL;
      rval.next=NULL;
   }else{
      for(p=head;p->next->next!=NULL;p=p->next);
      rval.next=NULL;
      rval.val=p->next->val;
      rval.dp=p->next->dp;
      p->next=NULL;
   }
}
      
return rval;
}


element *enqueue(element *head, double val, void *dp)
/*Adds element with value val to the queue, returning the head 
pointer.*/
{
element *newnode;

/*Create the new entry*/
newnode=(element *)R_alloc(1,sizeof(struct elementtype));

newnode->val=val;   /*Set the element value*/
newnode->dp=dp;

/*Set the next pointer equal to the current first entry (if any)*/
newnode->next=head;

/*Place the new node at the head of the queuek*/
head=newnode;

return head;
}

 
long int queuelen(element *head)
/*Returns the length of the queue pointed to by head*/
{
element *p; 
int count=0;

for(p=head;p!=NULL;p=p->next)
   count++;

return count;
}


char isinqueue(element *head,double val)
/*Returns a 1 if val is in the queue pointed to by head, otherwise 0*/
{
element *p;

for(p=head;p!=NULL;p=p->next)
   if(p->val==val)
      return 1;

return 0;
}


element queuedel(element *head,double val)
/*Find the element val in the queue pointed to by head and delete it,
returning the deleted element.*/
{
element rval,*p;

if(head==NULL){
   rval.val=-1.0;
   rval.dp=NULL;
   rval.next=NULL;
}else if(head->val==val){
   rval.val=head->val;
   rval.dp=head->dp;
   rval.next=head->next;
   head=rval.next;
}else{
   for(p=head;(p->next!=NULL)&&(p->next->val!=val);p=p->next);
   if(p->next==NULL){
      rval.val=-1.0;
      rval.dp=NULL;
      rval.next=NULL;
   }else{
      rval.val=p->next->val;
      rval.dp=p->next->dp;
      rval.next=p->next->next;
      p->next=rval.next;
   }
}
      
return rval;
}


