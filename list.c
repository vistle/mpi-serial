/*
 *  (C) 2000 UNIVERSITY OF CHICAGO
 *      See COPYRIGHT in top-level directory.
 */



#include <stdlib.h>
#include <stdio.h>
#include "listops.h"
#include "listP.h"

/*
 * list management code
 *
 * For storing singly-linked lists of pointers.
 *
 */


static int itemcount=0;
static int headcount=0;


/*
 * AP_listitem_malloc()
 *
 * malloc a new ilist item and return a pointer to it.
 *
 */

static pListitem AP_listitem_malloc(void)
{
  pListitem item;

  itemcount++;
  item=(pListitem)malloc( (unsigned) sizeof(Listitem) );

  if (!item)
    {
      perror("AP_listitem_malloc: malloc failure");
      abort();
    }

  return(item);
}



/*
 * AP_listitem_free(listitem)
 *
 * Free a listitem generated by AP_listitem_malloc()
 *
 */

static void AP_listitem_free(pListitem listitem)
{
  free(listitem);
  itemcount--;
}



/*
 * AP_listitem_verify(void)
 *
 * Checks to see if there are any outstanding listitems that have been
 * malloc'd.  Returns true if there are any.
 *
 */

int AP_listitem_verify(void)
{
  if (itemcount!=0)
    fprintf(stderr,"AP_list_verify: outstanding items, count=%d\n",
	    itemcount);

  if (headcount!=0)
    fprintf(stderr,"AP_list_verify: outstanding lists, count=%d\n",
	    headcount);

  return( (itemcount!=0) || (headcount!=0) );
}




pListitem AP_listitem_prev(pListitem listitem)
{
  return(listitem->prev);
}



pListitem AP_listitem_next(pListitem listitem)
{
  return(listitem->next);
}




void *AP_listitem_data(pListitem listitem)
{
  return(listitem->data );
}




/***************************************************************/



/*
 * AP_list_new(void)
 *
 * allocate an empty list return a pointer to it
 *
 */

pList AP_list_new(void)
{
  pList list;

  list=(pList)malloc(sizeof(List));

  if (!list)
    {
      perror("AP_list_new: malloc failure\n");
      abort();
    }

  list->head=NULL;
  list->tail=NULL;
  list->count=0;

  headcount++;
  return(list);
}





/*
 * AP_list_free(list)
 *
 * Free an entire list
 *
 */

void AP_list_free(pList list)
{
  pListitem next,cur;
  int count;

  count=0;
  cur=list->head;

  while(cur)
    {
      next=cur->next;

      AP_listitem_free(cur);
      count++;

      cur=next;
    }

  if (count!=list->count)
    {
      fprintf(stderr,"AP_list_free: count %d does not match actual length %d\n",
	      list->count,count);
      abort();
    }

  headcount--;
  free(list);
}



/*
 * AP_list_size(list)
 *
 * return the number of items in an ilist
 *
 */

int AP_list_size(pList list)
{
  return(list->count);
}



/*
 * AP_list_prepend(list,data)
 *
 * Prepend item to the front of list.
 *
 */

pListitem AP_list_prepend(pList list, void *data)
{
  pListitem new;

  new=AP_listitem_malloc();

  new->data=data;
  new->prev=NULL;
  new->next=list->head;

#ifdef CHECKS
  new->list=list;
#endif

  if (list->head)
    list->head->prev=new;

  list->head=new;
  if (!list->tail)
    list->tail=new;

  (list->count)++;

  return(new);
}



/*
 * AP_list_append(list,data)
 *
 * append item to end of list
 *
 */

pListitem AP_list_append(pList list, void *data)
{
  pListitem new;

  new=AP_listitem_malloc();
  new->data=data;
  new->prev=list->tail;
  new->next= NULL;

#ifdef CHECKS
  new->list= list;
#endif

  if (list->tail)
    list->tail->next=new;
  else
    list->head=new;

  list->tail=new;
  (list->count)++;

  return(new);
}





/*
 * AP_list_delete(list,data)
 *
 * delete item from list; return TRUE if successful
 *
 */

int AP_list_delete(pList list, void *data)
{
  pListitem item;

  if ((item=AP_list_search(list,data)))
    {
      AP_list_delete_item(list,item);
      return(1);
    }

  return(0);
}



void AP_list_delete_item(pList list, pListitem item)
{

#ifdef CHECKS
  if (item->list != list)
    {
      fprintf(stderr,"AP_list_delete_item: item is not in list\n");
      abort();
    }
#endif

  /* set pointer of prior listitem */

  if (item == list->head)
    list->head = item->next;
  else
    item->prev->next = item->next;

  /* set pointer of following listitem */

  if (item == list->tail)
    list->tail = item->prev;
  else
    item->next->prev = item->prev;

  AP_listitem_free(item);
  (list->count)--;
}




pListitem AP_list_head_item(pList list)
{
  return(list->head);
}



int AP_list_head(pList list, void **data)
{
  if (list->head)
    {
      *data=list->head->data;
      return(1);
    }
  else
    return(0);
}



int AP_list_tail(pList list, void **data)
{
  if (list->tail)
    {
      *data=list->tail->data;
      return(1);
    }
  else
    return(0);
}





/*
 * AP_list_print(str,list)
 *
 * Print out the message string followed by the
 * items in the list
 *
 */

void AP_list_print(char *str, pList list)
{
  pListitem cur;

  printf("%s (%d items): ",str,list->count);

  cur=list->head;
  while(cur)
    {
      printf("%ld ",(long int)cur->data);
      cur=cur->next;
    }

  printf("\n");
}




/*
 * AP_list_revprint(str,list)
 *
 * Print out the message string followed by the
 * items in the list
 *
 */

void AP_list_revprint(char *str, pList list)
{
  pListitem cur;

  printf("%s (%d items): ",str,list->count);

  cur=list->tail;
  while(cur)
    {
      printf("%ld ",(long int)cur->data);
      cur=cur->prev;
    }

  printf("\n");
}




/*
 * AP_list_search(list,data)
 *
 * Returns listitem if item appears in the list, otherwise NULL.
 *
 */


pListitem AP_list_search(pList list, void *data)
{
  pListitem cur;

  cur=list->head;

  while (cur)
    {
      if (cur->data == data)
        return(cur);

      cur=cur->next;
    }

  return(NULL);
}


/*
 * AP_list_search_func(list,func,data)
 *
 * Returns listitem if func(listitem->data,data) returns true
 *
 */


pListitem AP_list_search_func(pList list,
			      int (*func)(void *item_data, void *fixed_data),
			      void *fixed_data)
{
  pListitem cur;

  cur=list->head;

  while (cur)
    {
      if ( (*func)(cur->data,fixed_data) )
        return(cur);

      cur=cur->next;
    }

  return(NULL);
}



/*
 * AP_list_next(list,data,temp)
 *
 * like PList_next() except handles NULL pointers properly.
 *
 * initially, pass in (void **) NULL in 'temp'
 * returns next list item through 'item'
 * returns nonzero if there is a next item
 *
 */

int AP_list_next(pList list, void **data, void **temp)
{
  pListitem cur;

  if (*temp)                             /* temp is previous item */
    {
      cur=(pListitem)(*temp);
      cur=cur->next;
    }
  else                                  /* First item */
    cur=list->head;

  if (cur)
    {
      *temp=(void *)cur;
      *data=cur->data;
      return(1);
    }
  else
    return(0);
}


/*
 * Compatibility routine for scorec list traversal
 * Does not provide any way to differentiate
 * between NULL in the list, and the end of the list
 *
 */

void *AP_list_braindead_next(pList list, void **temp)
{
  void *item;

  if (AP_list_next(list,&item,temp))
    return(item);
  else
    return(NULL);
}



/*
 * AP_list_duplicate(list)
 *
 * return a copy of the list
 * (Note: caller is responsible for freeing this list)
 *
 */

pList AP_list_duplicate(pList list)
{
  pList newlist;
  pListitem cur,new,prev;

  newlist=AP_list_new();
  prev=NULL;

  cur=list->head;
  while(cur)
    {
      new=AP_listitem_malloc();
      new->data=cur->data;
      new->prev=prev;

      if (prev)
	prev->next=new;
      else
	newlist->head=new;

      prev=new;

      cur=cur->next;
    }

  if (prev)
    prev->next=NULL;

  newlist->tail=prev;
  newlist->count=list->count;
  return(newlist);
}



int AP_list_apply(pList list,
		  int (*func)(void *item_data, void *fixed_data),
		  void *fixed_data)
{
  pListitem cur;
  int total;

  total=0;
  cur=list->head;

  while (cur)
    {
      total += (*func)(cur->data,fixed_data);

      cur=cur->next;
    }

  return(total);
}




/*
 * main for debugging
 *
 */


#ifdef LISTMAIN

int main()
{
  pList mylist, list2;
  int i;
  void *temp,*item;
  pListitem next;

  mylist=AP_list_new();

  for (i=1; i<10; i++)
    {
      AP_list_prepend(mylist,(void *)i);
      AP_list_print("current",mylist);
      AP_list_revprint("    rev",mylist);
    }

  printf("Size %d\n",AP_list_size(mylist));

  for (i=10; i<15; i++)
    {
      AP_list_append(mylist,(void *)i);
      AP_list_print("new",mylist);
      AP_list_revprint("    rev",mylist);
    }

  AP_list_delete(mylist,(void *)5);
  AP_list_print("less 5",mylist);
  AP_list_revprint("  rev",mylist);

  AP_list_delete(mylist,(void *)9);
  AP_list_print("less 9",mylist);
  AP_list_revprint("  rev",mylist);

  AP_list_delete(mylist,(void *)14);
  AP_list_print("less 14",mylist);
  AP_list_revprint("  rev",mylist);

  AP_list_delete(mylist,(void *)2);
  AP_list_print("less 2",mylist);
  AP_list_revprint("  rev",mylist);

  if (!AP_list_delete(mylist,(void *)0))
    printf("(did not delete 0)\n");
  else
    printf("ERROR - found 0\n");
  AP_list_print("less 0",mylist);
  AP_list_revprint("  rev",mylist);

  if (AP_list_search(mylist,(void *)4))
    printf("Found 4\n");
  else
    printf("Did not find 4\n");

  if (AP_list_search(mylist,(void *)9))
    printf("Found 9\n");
  else
    printf("Did not find 9\n");

  printf("Traversal by AP_list_next()\n");
  temp=NULL;
  while (AP_list_next(mylist,&item,&temp))
    printf("  Got item %d\n",(int)item);

  printf("Traversal by AP_listitem_next()\n");
  for (item=AP_list_head_item(mylist); item; item=AP_listitem_next(item))
    printf("  Got item %d\n",(int)(AP_listitem_data(item)));


  list2=AP_list_duplicate(mylist);
  AP_list_print("Original list",mylist);
  AP_list_revprint("  rev",mylist);
  AP_list_print("Duplicate    ",list2);
  AP_list_revprint("  rev",list2);

  AP_list_append(list2,(void *)99);
  AP_list_print("Dup add 99   ",list2);
  AP_list_revprint("  rev",list2);


  printf("Traversal by AP_listitem_next(), deleting\n");
  i=0;
  for (item=AP_list_head_item(list2); item; )
    {
      printf("  Got item %d",(int)(AP_listitem_data(item)));

      next=AP_listitem_next(item);

      if (i%2)
	{
	  AP_list_delete_item(list2,item);
	  printf(" - deleted\n");
	}
      else
	printf("\n");

      item=next;
      i++;
    }

  AP_list_print("After delete-traversal",list2);

  AP_list_free(mylist);
  AP_list_print("After del    ",list2);
  AP_list_revprint("  rev",list2);

  AP_list_free(list2);

  AP_listitem_verify();

  return(0);
}
#endif
