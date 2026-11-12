#include <mpi.h>

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
