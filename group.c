
#include "mpiP.h"


/*********/


int FC_FUNC( mpi_group_incl, MPI_GROUP_INCL )
     (int *group, int *n, int *ranks, int *newgroup, int *ierror)
{
  *ierror= MPI_Group_incl(*group, *n, ranks, newgroup);
  return MPI_SUCCESS;
}


int MPI_Group_incl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup)
{

  if (group==MPI_GROUP_NULL)
    {
      fprintf(stderr,"MPI_Group_incl: null group passed in\n");
      abort();
    }

  if (group==MPI_GROUP_EMPTY || n==0)
    *newgroup=MPI_GROUP_EMPTY;
  else
    if (n==1 && ranks[0]==0)
      *newgroup=MPI_GROUP_ONE;
    else
      {
	fprintf(stderr,"MPI_Group_incl: more than 1 proc in group\n");
	abort();
      }

  return(MPI_SUCCESS);
}

int MPI_Group_excl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup)
{
  if (group==MPI_GROUP_NULL)
    {
      fprintf(stderr,"MPI_Group_excl: null group passed in\n");
      abort();
    }

    if (n>1)
      {
	fprintf(stderr,"MPI_Group_excl: more than 1 proc in group\n");
	abort();
      }

  if (group==MPI_GROUP_EMPTY || n==1)
    *newgroup=MPI_GROUP_EMPTY;
  else
      *newgroup=MPI_GROUP_ONE;

  return(MPI_SUCCESS);
}

int MPI_Group_size(MPI_Group group, int *size)
{
    if (group == MPI_GROUP_NULL)
        return MPI_ERR_GROUP;
    if (group == MPI_GROUP_EMPTY)
        *size = 0;
    else
        *size = 1;
  return(MPI_SUCCESS);
}

int MPI_Group_rank(MPI_Group group, int *rank)
{
    if (group == MPI_GROUP_NULL)
        return MPI_ERR_GROUP;
    if (group == MPI_GROUP_EMPTY)
        *rank = MPI_UNDEFINED;
    else
        *rank = 0;
  return(MPI_SUCCESS);
}




/*********/


/* MPI_Group_range_incl
 * Include a strided range of ranks in a group.  For one processor, if
 * "0" is included in any of these ranges, it can only be the first rank.
 * Thus, if rank 0 is specified, include it, otherwise use GROUP_NULL
 */


int FC_FUNC( mpi_group_range_incl, MPI_GROUP_RANGE_INCL )
     (int *group, int *n, int ranges[][3], int *newgroup, int *ierror)
{
  *ierror= MPI_Group_range_incl(*group, *n, ranges, newgroup);
  return MPI_SUCCESS;
}


int MPI_Group_range_incl(MPI_Group group, int n, int ranges[][3],
			 MPI_Group *newgroup)
{

  if (group==MPI_GROUP_NULL)
    {
      fprintf(stderr,"MPI_Group_range_incl: null group passed in\n");
      abort();
    }

  if (group==MPI_GROUP_EMPTY || n==0)
    *newgroup=MPI_GROUP_EMPTY;
  else
    if (n==1 && ranges[0][0]==0 && ranges[0][1]==0)
      *newgroup=MPI_GROUP_ONE;
    else
      {
	fprintf(stderr,"MPI_Group_range_incl: more than 1 proc in group\n");
	abort();
      }

  return(MPI_SUCCESS);
}




/*********/



int FC_FUNC( mpi_group_union, MPI_GROUP_UNION )
     (int *group1, int *group2, int *newgroup, int *ierror)
{
  *ierror= MPI_Group_union(*group1,*group2,newgroup);
  return MPI_SUCCESS;
}



int MPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup)
{

  if (group1==MPI_GROUP_NULL || group2==MPI_GROUP_NULL)
    {
      fprintf(stderr,"MPI_Group_union: null group passed in\n");
      abort();
    }

  if (group1==MPI_GROUP_ONE || group2==MPI_GROUP_ONE)
    *newgroup=MPI_GROUP_ONE;
  else
    *newgroup=MPI_GROUP_EMPTY;


  return(MPI_SUCCESS);
}

/*********/



int FC_FUNC( mpi_group_intersection, MPI_GROUP_INTERSECTION )
     (int *group1, int *group2, int *newgroup, int *ierror)
{
  *ierror= MPI_Group_intersection(*group1,*group2,newgroup);
  return MPI_SUCCESS;
}



int MPI_Group_intersection(MPI_Group group1, MPI_Group group2,
			   MPI_Group *newgroup)
{

  if (group1==MPI_GROUP_NULL || group2==MPI_GROUP_NULL)
    {
      fprintf(stderr,"MPI_Group_intersection: null group passed in\n");
      abort();
    }

  if (group1==MPI_GROUP_ONE && group2==MPI_GROUP_ONE)
    *newgroup=MPI_GROUP_ONE;
  else
    *newgroup=MPI_GROUP_EMPTY;


  return(MPI_SUCCESS);
}


/*********/



int FC_FUNC( mpi_group_difference, MPI_GROUP_DIFFERENCE )
     (int *group1, int *group2, int *newgroup, int *ierror)
{
  *ierror= MPI_Group_difference(*group1,*group2,newgroup);
  return MPI_SUCCESS;
}



int MPI_Group_difference(MPI_Group group1, MPI_Group group2,
			 MPI_Group *newgroup)
{

  if (group1==MPI_GROUP_NULL || group2==MPI_GROUP_NULL)
    {
      fprintf(stderr,"MPI_Group_intersection: null group passed in\n");
      abort();
    }

  if (group1==MPI_GROUP_EMPTY || group2==MPI_GROUP_ONE)
    *newgroup=MPI_GROUP_EMPTY;
  else
    *newgroup=MPI_GROUP_ONE;

  return(MPI_SUCCESS);
}



/*********/


int FC_FUNC( mpi_group_free, MPI_GROUP_FREE )(int *group, int *ierror)
{
  *ierror= MPI_Group_free(group);
  return MPI_SUCCESS;
}


int MPI_Group_free(MPI_Group *group)
{
  *group= MPI_GROUP_NULL;

  return(MPI_SUCCESS);
}


/*********/



int FC_FUNC( mpi_group_translate_ranks, MPI_GROUP_TRANSLATE_RANKS )
     ( int *group1, int *n, int *ranks1,
       int *group2, int *ranks2, int *ierror)
{
  *ierror= MPI_Group_translate_ranks(*group1,*n,ranks1,*group2,ranks2);
  return MPI_SUCCESS;
}



int MPI_Group_translate_ranks(MPI_Group group1, int n, const int ranks1[],
			      MPI_Group group2, int *ranks2)
{
  int i;

  if (group1==MPI_GROUP_NULL || group2==MPI_GROUP_NULL)
    {
      fprintf(stderr,"MPI_Group_translate_ranks: null group passed in\n");
      abort();
    }

  if (n==0)
    return(MPI_SUCCESS);

  if (group1==MPI_GROUP_EMPTY)
    {
      fprintf(stderr,"MPI_Group_translate_ranks: empty input group\n");
      abort();
    }

  for (i=0; i<n; i++)
    {
      if (ranks1[i]!=0)
	{
	  fprintf(stderr,"MPI_Group_translate_ranks: bad input rank: %d\n",
		  ranks1[i]);
	  abort();
	}

      if (group2!=MPI_GROUP_EMPTY)
	ranks2[i]=ranks1[i];
      else
	ranks2[i]=MPI_UNDEFINED;
    }


  return(MPI_SUCCESS);

}

int MPI_Group_compare(MPI_Comm group1, MPI_Comm group2, int *result)
{
  if (group1==MPI_GROUP_NULL || group2==MPI_GROUP_NULL)
  {
      fprintf(stderr,"MPI_Group_compare: null input group\n");
      abort();
      return MPI_ERR_GROUP;
  }

  if (group1==MPI_GROUP_EMPTY && group2==MPI_GROUP_EMPTY)
      *result = MPI_IDENT;
  else if (group1 != MPI_GROUP_EMPTY && group2 != MPI_GROUP_EMPTY)
      *result = MPI_IDENT;
  else
      *result = MPI_UNEQUAL;
  return MPI_SUCCESS;
}


/*********/


MPI_Group MPI_Group_f2c(MPI_Fint group)
{
  return(group);
}


/*********/


MPI_Fint MPI_Group_c2f(MPI_Group group)
{
  return(group);
}

