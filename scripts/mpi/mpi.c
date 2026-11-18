#include <mpi.h>
#include <stddef.h>

int mpi_init(void)
{
    int initialized;
    MPI_Initialized(&initialized);

    if (!initialized) { return MPI_Init(NULL, NULL); }

    return MPI_SUCCESS;
}

int mpi_barrier(void)
{ return MPI_Barrier(MPI_COMM_WORLD); }

int mpi_rank(int *rank)
{ return MPI_Comm_rank(MPI_COMM_WORLD, rank); }

int mpi_size(int *size)
{ return MPI_Comm_size(MPI_COMM_WORLD, size); }

int mpi_allgather_int( int sendvalue, int *recvvalues)
{ return MPI_Allgather( &sendvalue, 1, MPI_INT, recvvalues, 1, MPI_INT, MPI_COMM_WORLD); }

int mpi_bcast_int(int *value, int root)
{ return MPI_Bcast( value, 1, MPI_INT, root, MPI_COMM_WORLD); }

/* Gather variable-size byte buffers */
int mpi_gatherv_bytes( const char *sendbuf, int sendbytes, char *recvbuf, const int *recvcounts, const int *displs, int root)
{ return MPI_Gatherv( sendbuf, sendbytes, MPI_BYTE, recvbuf, recvcounts, displs, MPI_BYTE, root, MPI_COMM_WORLD); }

/* Allgather variable-length byte buffers */
int mpi_allgatherv_bytes( const char *sendbuf, int sendbytes, char *recvbuf, const int *recvcounts, const int *displs)
{ return MPI_Allgatherv( sendbuf, sendbytes, MPI_BYTE, recvbuf, recvcounts, displs, MPI_BYTE, MPI_COMM_WORLD); }

/* Broadcast byte buffer */
int mpi_bcast_bytes( char *buffer, int nbytes, int root)
{ return MPI_Bcast( buffer, nbytes, MPI_BYTE, root, MPI_COMM_WORLD); }

/* scatterv  buffer */
int mpi_scatter_int(int *sendbuf, int *recvbuf, int root)
{ return MPI_Scatter( sendbuf, 1, MPI_INT, recvbuf, 1, MPI_INT, root, MPI_COMM_WORLD); }

int mpi_scatterv_bytes( char *sendbuf, int *sendcounts, int *displs, char *recvbuf, int recvcount, int root)
{ return MPI_Scatterv( sendbuf, sendcounts, displs, MPI_BYTE, recvbuf, recvcount, MPI_BYTE, root, MPI_COMM_WORLD); }

/* Reduce double arrays */
int mpi_reduce_double( const double *sendbuf, double *recvbuf, int count, int operation, int root)
{
    MPI_Op op;
    switch (operation) {
        case 0: op = MPI_SUM; break;
        case 1: op = MPI_MAX; break;
        case 2: op = MPI_MIN; break;
        case 3: op = MPI_PROD; break;
        default: return MPI_ERR_OP;
    }
    return MPI_Reduce( sendbuf, recvbuf, count, MPI_DOUBLE, op, root, MPI_COMM_WORLD);
}

/* Allreduce integer message sizes */
int mpi_allreduce_double( const double *sendbuf, double *recvbuf, int count, int operation)
{ 
    MPI_Op op;
    switch (operation) {
        case 0: op = MPI_SUM; break;
        case 1: op = MPI_MAX; break;
        case 2: op = MPI_MIN; break;
        case 3: op = MPI_PROD; break;
        default: return MPI_ERR_OP;
    }
    return MPI_Allreduce( sendbuf, recvbuf, count, MPI_DOUBLE, op, MPI_COMM_WORLD);
}

int mpi_finalize(void)
{
    int initialized, finalized;
    MPI_Initialized(&initialized);
    MPI_Finalized(&finalized);
    if (initialized && !finalized) { return MPI_Finalize(); }
    return MPI_SUCCESS;
}
