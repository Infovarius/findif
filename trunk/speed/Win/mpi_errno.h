/********************************************************************************

  $Workfile: mpi_errno.h $ 

  $Author: Hernani $ 
  $Modtime: 7-07-99 16:59 $ 
  $Revision: 1 $ 


********************************************************************************/

/* error codes for MPI programs*/

#ifndef MPI_ERRNO_H_
#define MPI_ERRNO_H_

/* error return classes */
#define MPI_SUCCESS          0      /* Successful return code */
#define MPI_ERR_BUFFER       1      /* Invalid buffer pointer */
#define MPI_ERR_COUNT        2      /* Invalid count argument */
#define MPI_ERR_TYPE         3      /* Invalid datatype argument */
#define MPI_ERR_TAG          4      /* Invalid tag argument */
#define MPI_ERR_COMM         5      /* Invalid communicator */
#define MPI_ERR_RANK         6      /* Invalid rank */
#define MPI_ERR_ROOT         7      /* Invalid root */
#define MPI_ERR_GROUP        8      /* Null group passed to function */
#define MPI_ERR_OP           9      /* Invalid operation */
#define MPI_ERR_TOPOLOGY    10      /* Invalid topology */
#define MPI_ERR_DIMS        11      /* Illegal dimension argument */
#define MPI_ERR_ARG         12      /* Invalid argument */
#define MPI_ERR_UNKNOWN     13      /* Unknown error */
#define MPI_ERR_TRUNCATE    14      /* message truncated on receive */
#define MPI_ERR_OTHER       15      /* Other error; use Error_string */
#define MPI_ERR_INTERN      16      /* internal error code    */
#define MPI_ERR_IN_STATUS   17      /* Look in status for error value */
#define MPI_ERR_PENDING     18      /* Pending request */
#define MPI_ERR_REQUEST     19      /* illegal mpi_request handle */
#define MPI_ERR_LASTCODE    (256*16+18)      /* Last error code*/


#endif //MPI_ERRNO_H_




