/********************************************************************************

  $Workfile: binding.h $ 

  $Author: Hernani $ 
  $Modtime: 9-12-99 12:05 $ 
  $Revision: 5 $ 


********************************************************************************/

#ifndef __MPI_BINDINGS
#define __MPI_BINDINGS

#include "mpi.h"

__declspec( dllexport ) int MPI_Send(void*, int, MPI_Datatype, int, int, MPI_Comm);
__declspec( dllexport ) int MPI_Recv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *);
__declspec( dllexport ) int MPI_Get_count(MPI_Status *, MPI_Datatype, int *);
__declspec( dllexport ) int MPI_Bsend(void*, int, MPI_Datatype, int, int, MPI_Comm);
__declspec( dllexport ) int MPI_Ssend(void*, int, MPI_Datatype, int, int, MPI_Comm);
__declspec( dllexport ) int MPI_Rsend(void*, int, MPI_Datatype, int, int, MPI_Comm);
__declspec( dllexport ) int MPI_Buffer_attach( void*, int);
__declspec( dllexport ) int MPI_Buffer_detach( void*, int*);
__declspec( dllexport ) int MPI_Isend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
__declspec( dllexport ) int MPI_Ibsend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
__declspec( dllexport ) int MPI_Issend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
__declspec( dllexport ) int MPI_Irsend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
__declspec( dllexport ) int MPI_Irecv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
__declspec( dllexport ) int MPI_Wait(MPI_Request *, MPI_Status *);
__declspec( dllexport ) int MPI_Test(MPI_Request *, int *, MPI_Status *);
__declspec( dllexport ) int MPI_Request_free(MPI_Request *);
__declspec( dllexport ) int MPI_Waitany(int, MPI_Request *, int *, MPI_Status *);
__declspec( dllexport ) int MPI_Testany(int, MPI_Request *, int *, int *, MPI_Status *);
__declspec( dllexport ) int MPI_Waitall(int, MPI_Request *, MPI_Status *);
__declspec( dllexport ) int MPI_Testall(int, MPI_Request *, int *, MPI_Status *);
__declspec( dllexport ) int MPI_Waitsome(int, MPI_Request *, int *, int *, MPI_Status *);
__declspec( dllexport ) int MPI_Testsome(int, MPI_Request *, int *, int *, MPI_Status *);
__declspec( dllexport ) int MPI_Iprobe(int, int, MPI_Comm, int *flag, MPI_Status *);
__declspec( dllexport ) int MPI_Probe(int, int, MPI_Comm, MPI_Status *);
__declspec( dllexport ) int MPI_Cancel(MPI_Request *);
__declspec( dllexport ) int MPI_Test_cancelled(MPI_Status *, int *flag);
__declspec( dllexport ) int MPI_Send_init(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
__declspec( dllexport ) int MPI_Bsend_init(void*, int, MPI_Datatype, int,int, MPI_Comm, MPI_Request *);
__declspec( dllexport ) int MPI_Ssend_init(void*, int, MPI_Datatype, int,int, MPI_Comm, MPI_Request *);
__declspec( dllexport ) int MPI_Rsend_init(void*, int, MPI_Datatype, int,int, MPI_Comm, MPI_Request *);
__declspec( dllexport ) int MPI_Recv_init(void*, int, MPI_Datatype, int,int, MPI_Comm, MPI_Request *);
__declspec( dllexport ) int MPI_Start(MPI_Request *);
__declspec( dllexport ) int MPI_Startall(int, MPI_Request *);
__declspec( dllexport ) int MPI_Sendrecv(void *, int, MPI_Datatype,int, int, void *, int, 
		 MPI_Datatype, int, int, MPI_Comm, MPI_Status *);
__declspec( dllexport ) int MPI_Sendrecv_replace(void*, int, MPI_Datatype, 
			 int, int, int, int, MPI_Comm, MPI_Status *);
__declspec( dllexport ) int MPI_Type_contiguous(int, MPI_Datatype, MPI_Datatype *);
__declspec( dllexport ) int MPI_Type_vector(int, int, int, MPI_Datatype, MPI_Datatype *);
__declspec( dllexport ) int MPI_Type_hvector(int, int, MPI_Aint, MPI_Datatype, MPI_Datatype *);
__declspec( dllexport ) int MPI_Type_indexed(int, int *, int *, MPI_Datatype, MPI_Datatype *);
__declspec( dllexport ) int MPI_Type_hindexed(int, int *, MPI_Aint *, MPI_Datatype, MPI_Datatype *);
__declspec( dllexport ) int MPI_Type_struct(int, int *, MPI_Aint *, MPI_Datatype *, MPI_Datatype *);
__declspec( dllexport ) int MPI_Address(void*, MPI_Aint *);
__declspec( dllexport ) int MPI_Type_extent(MPI_Datatype, MPI_Aint *);

__declspec( dllexport ) int MPI_Type_size(MPI_Datatype, int *);
__declspec( dllexport ) int MPI_Type_count(MPI_Datatype, int *);
__declspec( dllexport ) int MPI_Type_lb(MPI_Datatype, MPI_Aint*);
__declspec( dllexport ) int MPI_Type_ub(MPI_Datatype, MPI_Aint*);
__declspec( dllexport ) int MPI_Type_commit(MPI_Datatype *);
__declspec( dllexport ) int MPI_Type_free(MPI_Datatype *);
__declspec( dllexport ) int MPI_Get_elements(MPI_Status *, MPI_Datatype, int *);
__declspec( dllexport ) int MPI_Pack(void* inbuf, int, MPI_Datatype, void *outbuf, 
	     int outsize, int *position,  MPI_Comm);
__declspec( dllexport ) int MPI_Unpack(void* inbuf, int insize, int *position, void *outbuf, 
	       int , MPI_Datatype, MPI_Comm);
__declspec( dllexport ) int MPI_Pack_size(int, MPI_Datatype, MPI_Comm, 
		  int *);
__declspec( dllexport ) int MPI_Barrier(MPI_Comm );

__declspec( dllexport ) int MPI_Bcast(void*fer, int, MPI_Datatype, int root, 
	      MPI_Comm );
__declspec( dllexport ) int MPI_Gather(void* , int, MPI_Datatype, 
	       void*, int, MPI_Datatype, 
	       int root, MPI_Comm); 
__declspec( dllexport ) int MPI_Gatherv(void* , int, MPI_Datatype, 
		void*, int *recvcounts, int *displs, 
		MPI_Datatype, int root, MPI_Comm); 
__declspec( dllexport ) int MPI_Scatter(void* , int, MPI_Datatype, 
		void*, int, MPI_Datatype, 
		int root, MPI_Comm);
__declspec( dllexport ) int MPI_Scatterv(void* , int *sendcounts, int *displs, 
		 MPI_Datatype, void*, int, 
		 MPI_Datatype, int root, MPI_Comm);
__declspec( dllexport ) int MPI_Allgather(void* , int, MPI_Datatype, 
		  void*, int, MPI_Datatype, 
		  MPI_Comm);
__declspec( dllexport ) int MPI_Allgatherv(void* , int, MPI_Datatype, 
		   void*, int *recvcounts, int *displs, 
		   MPI_Datatype, MPI_Comm);
__declspec( dllexport ) int MPI_Alltoall(void* , int, MPI_Datatype, 
		 void*, int, MPI_Datatype, 
		 MPI_Comm);
__declspec( dllexport ) int MPI_Alltoallv(void* , int *sendcounts, int *sdispls, 
		  MPI_Datatype, void*, int *recvcounts, 
		  int *rdispls, MPI_Datatype, MPI_Comm);
__declspec( dllexport ) int MPI_Reduce(void* , void*, int, 
	       MPI_Datatype, MPI_Op op, int root, MPI_Comm);
__declspec( dllexport ) int MPI_Op_create(MPI_User_function *, int, MPI_Op *);
__declspec( dllexport ) int MPI_Op_free( MPI_Op *);
__declspec( dllexport ) int MPI_Allreduce(void* , void*, int, 
		  MPI_Datatype, MPI_Op op, MPI_Comm);
__declspec( dllexport ) int MPI_Reduce_scatter(void* , void*, int *recvcounts, 
		       MPI_Datatype, MPI_Op op, MPI_Comm);
__declspec( dllexport ) int MPI_Scan(void* , void*, int, MPI_Datatype, 
	     MPI_Op op, MPI_Comm );
__declspec( dllexport ) int MPI_Group_size(MPI_Group group, int *);
__declspec( dllexport ) int MPI_Group_rank(MPI_Group group, int *rank);
__declspec( dllexport ) int MPI_Group_translate_ranks (MPI_Group group1, int n, int *ranks1, 
			       MPI_Group group2, int *ranks2);
__declspec( dllexport ) int MPI_Group_compare(MPI_Group group1, MPI_Group group2, int *result);
__declspec( dllexport ) int MPI_Comm_group(MPI_Comm, MPI_Group *);
__declspec( dllexport ) int MPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup);
__declspec( dllexport ) int MPI_Group_intersection(MPI_Group group1, MPI_Group group2, 
			   MPI_Group *newgroup);
__declspec( dllexport ) int MPI_Group_difference(MPI_Group group1, MPI_Group group2, 
			 MPI_Group *newgroup);
__declspec( dllexport ) int MPI_Group_incl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup);
__declspec( dllexport ) int MPI_Group_excl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup);
__declspec( dllexport ) int MPI_Group_range_incl(MPI_Group group, int n, int ranges[][3], 
			 MPI_Group *newgroup);
__declspec( dllexport ) int MPI_Group_range_excl(MPI_Group group, int n, int ranges[][3], 
			 MPI_Group *newgroup);
__declspec( dllexport ) int MPI_Group_free(MPI_Group *);
__declspec( dllexport ) int MPI_Comm_size(MPI_Comm, int *);
__declspec( dllexport ) int MPI_Comm_rank(MPI_Comm, int *rank);
__declspec( dllexport ) int MPI_Comm_compare(MPI_Comm, MPI_Comm, int *result);
__declspec( dllexport ) int MPI_Comm_dup(MPI_Comm, MPI_Comm *newcomm);
__declspec( dllexport ) int MPI_Comm_create(MPI_Comm, MPI_Group group, MPI_Comm *newcomm);
__declspec( dllexport ) int MPI_Comm_split(MPI_Comm, int color, int key, MPI_Comm *newcomm);
__declspec( dllexport ) int MPI_Comm_free(MPI_Comm *comm);
__declspec( dllexport ) int MPI_Comm_test_inter(MPI_Comm, int *flag);
__declspec( dllexport ) int MPI_Comm_remote_size(MPI_Comm, int *);
__declspec( dllexport ) int MPI_Comm_remote_group(MPI_Comm, MPI_Group *);
__declspec( dllexport ) int MPI_Intercomm_create(MPI_Comm local_comm, int local_leader, 
			 MPI_Comm peer_comm, int remote_leader, 
			 int, MPI_Comm *newintercomm);
__declspec( dllexport ) int MPI_Intercomm_merge(MPI_Comm intercomm, int high, MPI_Comm *newintracomm);
__declspec( dllexport ) int MPI_Keyval_create(MPI_Copy_function *copy_fn, 
		      MPI_Delete_function *delete_fn, 
		      int *keyval, void* extra_state);
__declspec( dllexport ) int MPI_Keyval_free(int *keyval);
__declspec( dllexport ) int MPI_Attr_put(MPI_Comm, int keyval, void* attribute_val);
__declspec( dllexport ) int MPI_Attr_get(MPI_Comm, int keyval, void *attribute_val, int *flag);
__declspec( dllexport ) int MPI_Attr_delete(MPI_Comm, int keyval);
__declspec( dllexport ) int MPI_Topo_test(MPI_Comm, int *);
__declspec( dllexport ) int MPI_Cart_create(MPI_Comm comm_old, int ndims, int *dims, int *periods,
		    int reorder, MPI_Comm *comm_cart);
__declspec( dllexport ) int MPI_Dims_create(int nnodes, int ndims, int *dims);
__declspec( dllexport ) int MPI_Graph_create(MPI_Comm, int, int *, int *, int, MPI_Comm *);
__declspec( dllexport ) int MPI_Graphdims_get(MPI_Comm, int *nnodes, int *nedges);
__declspec( dllexport ) int MPI_Graph_get(MPI_Comm, int, int, int *, int *);
__declspec( dllexport ) int MPI_Cartdim_get(MPI_Comm, int *ndims);
__declspec( dllexport ) int MPI_Cart_get(MPI_Comm, int maxdims, int *dims, int *periods,
		 int *coords);
__declspec( dllexport ) int MPI_Cart_rank(MPI_Comm, int *coords, int *rank);
__declspec( dllexport ) int MPI_Cart_coords(MPI_Comm, int rank, int maxdims, int *coords);
__declspec( dllexport ) int MPI_Graph_neighbors_count(MPI_Comm, int rank, int *nneighbors);
__declspec( dllexport ) int MPI_Graph_neighbors(MPI_Comm, int rank, int maxneighbors,
			int *neighbors);
__declspec( dllexport ) int MPI_Cart_shift(MPI_Comm, int direction, int disp, 
		   int *rank_source, int *rank_dest);
__declspec( dllexport ) int MPI_Cart_sub(MPI_Comm, int *remain_dims, MPI_Comm *newcomm);
__declspec( dllexport ) int MPI_Cart_map(MPI_Comm, int ndims, int *dims, int *periods, 
		 int *newrank);
__declspec( dllexport ) int MPI_Graph_map(MPI_Comm, int, int *, int *, int *);
__declspec( dllexport ) int MPI_Get_processor_name(char *name, int *result_len);
__declspec( dllexport ) int MPI_Errhandler_create(MPI_Handler_function *function, 
			  MPI_Errhandler *errhandler);
__declspec( dllexport ) int MPI_Errhandler_set(MPI_Comm, MPI_Errhandler errhandler);
__declspec( dllexport ) int MPI_Errhandler_get(MPI_Comm, MPI_Errhandler *errhandler);
__declspec( dllexport ) int MPI_Errhandler_free(MPI_Errhandler *errhandler);
__declspec( dllexport ) int MPI_Error_string(int errorcode, char *string, int *result_len);
__declspec( dllexport ) int MPI_Error_class(int errorcode, int *errorclass);
__declspec( dllexport ) double MPI_Wtime(void);
__declspec( dllexport ) double MPI_Wtick(void);
__declspec( dllexport ) int MPI_Init(int *argc, char ***argv);
__declspec( dllexport ) int MPI_Finalize(void);
__declspec( dllexport ) int MPI_Initialized(int *flag);
__declspec( dllexport ) int MPI_Abort(MPI_Comm, int errorcode);

/* MPI-2 communicator naming functions */
__declspec( dllexport ) int MPI_Comm_set_name(MPI_Comm, char *);
__declspec( dllexport ) int MPI_Comm_get_name(MPI_Comm, char **);


/* ---------------------- MPI-2 MPI and Threads --------------------------------- */
_declspec (dllexport) int MPI_Init_thread
	(
    int *,
    char ***,
    int,
    int *
	);

_declspec (dllexport) int MPI_Query_thread
	(
    int *
	);

_declspec (dllexport) int MPI_Is_thread_main
	(
    int *
	);
/* ------------------------------------------------------------------------------ */

__declspec( dllexport ) int MPI_NULL_COPY_FN ( MPI_Comm oldcomm, int keyval, void *extra_state, 
		       void *attr_in, void *attr_out, int *flag );
__declspec( dllexport ) int MPI_NULL_DELETE_FN ( MPI_Comm, int keyval, void *attr, 
			 void *extra_state );
__declspec( dllexport ) int MPI_DUP_FN ( MPI_Comm, int keyval, void *extra_state, void *attr_in,
		 void *attr_out, int *flag );

#endif
