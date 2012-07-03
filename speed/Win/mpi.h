/********************************************************************************

  $Workfile: mpi.h $ 

  $Author: Hernani $ 
  $Modtime: 9-12-99 12:05 $ 
  $Revision: 5 $ 


********************************************************************************/

/* user include file for MPI programs */

#ifndef MPI_H_
#define MPI_H_


/* Keep C++ compilers from getting confused */
#if defined(__cplusplus)
extern "C" {
#endif


 /* Results of the compare operations */
/* These should stay ordered */
#define MPI_IDENT     0
#define MPI_CONGRUENT 1
#define MPI_SIMILAR   2
#define MPI_UNEQUAL   3

/* Data types
 * A more aggressive yet homogeneous implementation might want to 
 * make the values here the number of bytes in the basic type, with
 * a simple test against a max limit (e.g., 16 for long double), and
 * non-contiguous structures with indices greater than that.
 */
typedef int MPI_Datatype;
#define MPI_CHAR           ((MPI_Datatype)1)
#define MPI_UNSIGNED_CHAR  ((MPI_Datatype)2)
#define MPI_BYTE           ((MPI_Datatype)3)
#define MPI_SHORT          ((MPI_Datatype)4)
#define MPI_UNSIGNED_SHORT ((MPI_Datatype)5)
#define MPI_INT            ((MPI_Datatype)6)
#define MPI_UNSIGNED       ((MPI_Datatype)7)
#define MPI_LONG           ((MPI_Datatype)8)
#define MPI_UNSIGNED_LONG  ((MPI_Datatype)9)
#define MPI_FLOAT          ((MPI_Datatype)10)
#define MPI_DOUBLE         ((MPI_Datatype)11)
#define MPI_LONG_DOUBLE    ((MPI_Datatype)12)
#define MPI_LONG_LONG_INT  ((MPI_Datatype)13)

#define MPI_PACKED         ((MPI_Datatype)14)
#define MPI_LB             ((MPI_Datatype)15)
#define MPI_UB             ((MPI_Datatype)16)

/* 
   The layouts for the types MPI_DOUBLE_INT etc are simply
   struct { 
       double var;
       int    loc;
   }
   This is documented in the man pages on the various datatypes.   
 */
#define MPI_FLOAT_INT      ((MPI_Datatype)17)
#define MPI_DOUBLE_INT     ((MPI_Datatype)18)
#define MPI_LONG_INT       ((MPI_Datatype)19)
#define MPI_SHORT_INT      ((MPI_Datatype)20)
#define MPI_2INT           ((MPI_Datatype)21)
#define MPI_LONG_DOUBLE_INT ((MPI_Datatype)22)

/* Fortran types */
#define MPI_COMPLEX        ((MPI_Datatype)23)
#define MPI_DOUBLE_COMPLEX ((MPI_Datatype)24)
#define MPI_LOGICAL        ((MPI_Datatype)25)
#define MPI_REAL           ((MPI_Datatype)26)
#define MPI_DOUBLE_PRECISION ((MPI_Datatype)27)
#define MPI_INTEGER        ((MPI_Datatype)28)
#define MPI_2INTEGER       ((MPI_Datatype)29)
#define MPI_2COMPLEX       ((MPI_Datatype)30)
#define MPI_2DOUBLE_COMPLEX   ((MPI_Datatype)31)
#define MPI_2REAL             ((MPI_Datatype)32)
#define MPI_2DOUBLE_PRECISION ((MPI_Datatype)33)
#define MPI_CHARACTER         ((MPI_Datatype)1)

/* Communicators */
typedef int MPI_Comm;
#define MPI_COMM_WORLD 91
#define MPI_COMM_SELF 92

/* Groups */
typedef int MPI_Group;
#define MPI_GROUP_EMPTY 90

/* Collective operations */
typedef int MPI_Op;

#define MPI_MAX    (MPI_Op)(100)
#define MPI_MIN    (MPI_Op)(101)
#define MPI_SUM    (MPI_Op)(102)
#define MPI_PROD   (MPI_Op)(103)
#define MPI_LAND   (MPI_Op)(104)
#define MPI_BAND   (MPI_Op)(105)
#define MPI_LOR    (MPI_Op)(106)
#define MPI_BOR    (MPI_Op)(107)
#define MPI_LXOR   (MPI_Op)(108)
#define MPI_BXOR   (MPI_Op)(109)
#define MPI_MINLOC (MPI_Op)(110)
#define MPI_MAXLOC (MPI_Op)(111)

/* Permanent key values */
/* C Versions (return pointer to value) */
#define MPI_TAG_UB 81
#define MPI_HOST 83
#define MPI_IO 85
#define MPI_WTIME_IS_GLOBAL 87
/* Fortran Versions (return integer value) */
#define WMLF_TAG_UB 80
#define WML_HOST 82
#define WML_IO 84
#define WML_WTIME_IS_GLOBAL 86

/* Define some null objects */
#define MPI_COMM_NULL      ((MPI_Comm)0)
#define MPI_OP_NULL        ((MPI_Op)0)
#define MPI_GROUP_NULL     ((MPI_Group)0)
#define MPI_DATATYPE_NULL  ((MPI_Datatype)0)
#define MPI_REQUEST_NULL   ((MPI_Request)0)
#define MPI_ERRHANDLER_NULL ((MPI_Errhandler )0)

/* These are only guesses*/
#define MPI_MAX_PROCESSOR_NAME 256
#define MPI_MAX_ERROR_STRING   512
#define MPI_MAX_NAME_STRING     63		

/* Pre-defined constants */
#define MPI_UNDEFINED      (-32766)
#define MPI_UNDEFINED_RANK MPI_UNDEFINED
#define MPI_KEYVAL_INVALID 0

/* Upper bound on the overhead in bsend for each message buffer */
#define MPI_BSEND_OVERHEAD 512

/* Topology types */
#define MPI_GRAPH  1
#define MPI_CART   2

#define MPI_BOTTOM      (void *)0

#define MPI_PROC_NULL   (-1)
#define MPI_ANY_SOURCE 	(-2)
#define MPI_ANY_TAG	(-1)

/* 
   Status object.  It is the only user-visible MPI data-structure 
   The "count" field is PRIVATE; use MPI_Get_count to access it. 
 */
typedef struct { 
    int count;
    int MPI_SOURCE;
    int MPI_TAG;
    int MPI_ERROR;
} MPI_Status;

typedef long MPI_Aint;


/* MPI Error handlers.  Systems that don't support stdargs can't use
   this definition
 */
typedef void (MPI_Handler_function) ();


#define MPI_ERRORS_ARE_FATAL ((MPI_Errhandler)119)
#define MPI_ERRORS_RETURN    ((MPI_Errhandler)120)
#define WML_ERRORS_WARN     ((MPI_Errhandler)121)
typedef int MPI_Errhandler;
/* Make the C names for the null functions all upper-case.  Note that 
   this is required for systems that use all uppercase names for Fortran 
   externals.  */
#define MPI_NULL_COPY_FN   WML_null_copy_fn
#define MPI_NULL_DELETE_FN WML_null_delete_fn
#define MPI_DUP_FN         WML_dup_fn

/* MPI request opjects */
typedef union WML_Handle_ *MPI_Request;

/* User combination function */
typedef void (MPI_User_function) ( void *, void *, int *, MPI_Datatype *); 

/* MPI Attribute copy and delete functions */
typedef int (MPI_Copy_function) ( MPI_Comm, int, void *, void *, void *, int *);
typedef int (MPI_Delete_function) ( MPI_Comm, int, void *, void * );

#define MPI_VERSION 1
#define MPI_SUBVERSION 1


/* for info */
typedef struct WML_Info *MPI_Info;
# define MPI_INFO_NULL         ((MPI_Info) 0)
# define MPI_MAX_INFO_KEY       255
# define MPI_MAX_INFO_VAL      1024

#define MPI_THREAD_SINGLE     0
#define MPI_THREAD_FUNNELED   1
#define MPI_THREAD_SERIALIZED 2
#define MPI_THREAD_MULTIPLE   3

/********************** MPI-2 FEATURES END HERE ***************************/

/* MPI's error classes */
#include "mpi_errno.h"

/* Bindings of the MPI routines */
#include "binding.h"


#if defined(__cplusplus)
}
#endif


#endif MPI_H_
