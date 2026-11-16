#include <mpi.h>
#include <stddef.h>

int mpi_init(void)
{
    int initialized;
    MPI_Initialized(&initialized);

    if (!initialized) {
        return MPI_Init(NULL, NULL);
    }

    return MPI_SUCCESS;
}

int mpi_barrier(void)
{
    return MPI_Barrier(MPI_COMM_WORLD);
}

int mpi_rank(int *rank)
{
    return MPI_Comm_rank(MPI_COMM_WORLD, rank);
}

int mpi_size(int *size)
{
    return MPI_Comm_size(MPI_COMM_WORLD, size);
}

int mpi_allgather_int(
    int sendvalue,
    int *recvvalues
)
{
    return MPI_Allgather(
        &sendvalue, 1, MPI_INT,
        recvvalues, 1, MPI_INT,
        MPI_COMM_WORLD
    );
}

/* Gather arbitrary byte buffers */
int mpi_gather_bytes(
    const char *sendbuf,
    int sendbytes,
    char *recvbuf,
    int root
)
{
    return MPI_Gather(
        sendbuf,
        sendbytes,
        MPI_BYTE,
        recvbuf,
        sendbytes,
        MPI_BYTE,
        root,
        MPI_COMM_WORLD
    );
}

/* Gather variable-size byte buffers */
int mpi_gatherv_bytes(
    const char *sendbuf,
    int sendbytes,
    char *recvbuf,
    const int *recvcounts,
    const int *displs,
    int root
)
{
    return MPI_Gatherv(
        sendbuf,
        sendbytes,
        MPI_BYTE,
        recvbuf,
        recvcounts,
        displs,
        MPI_BYTE,
        root,
        MPI_COMM_WORLD
    );
}

int mpi_finalize(void)
{
    int initialized, finalized;

    MPI_Initialized(&initialized);
    MPI_Finalized(&finalized);

    if (initialized && !finalized) {
        return MPI_Finalize();
    }

    return MPI_SUCCESS;
}
