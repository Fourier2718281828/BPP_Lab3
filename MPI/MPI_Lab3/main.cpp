#include <mpi.h>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <array>
#include <vector>
#include <stdlib.h>

#define PRINTER if (is_main_process) std::cout 

#define PRINT_MATRIX(chunk, M, N)                \
{                                                \
    for (std::size_t i = 0u; i < M; ++i)         \
    {                                            \
        for (std::size_t j = 0u; j < N; ++j)     \
        {                                        \
            PRINTER << chunk[i * N + j] << ' ';  \
        }                                        \
        PRINTER << '\n';                         \
    }                                            \
}

using std::cout;

using value_type = float;
using matrix_type = std::vector<value_type>;
using row_type = std::vector<value_type>;
using fixed_values_type = std::vector<std::tuple<std::size_t, std::size_t, value_type>>;

std::size_t M = 0u;
std::size_t N = 0u;


matrix_type get_zero_matrix(const std::size_t length)
{
    return matrix_type(length);
}

fixed_values_type get_default_fixed(const std::size_t m, const std::size_t n)
{
    fixed_values_type fixed_values;

    for (std::size_t i = 0; i < m; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            if (i == 0 || i == m - 1 || j == 0 || j == n - 1)
            {
                fixed_values.emplace_back(i, j, 1.f);
            }
        }
    }

    fixed_values.emplace_back(m >> 1, n >> 1, 1.f);
    return fixed_values;
}

void set_fixed(matrix_type& matrix, const fixed_values_type& fixed_values)
{
    for (const auto& [x, y, value] : fixed_values)
    {
        matrix[x * N + y] = value;
    }
}

inline bool is_index_within_rank(const std::size_t i, const std::size_t j, const int rank, const std::size_t chunk_size)
{
    const std::size_t urank = static_cast<std::size_t>(rank);
    const std::size_t first_index = (urank) * chunk_size;
    const std::size_t last_index  = (urank + 1) * chunk_size;
    const std::size_t linear_index = i * N + j;
    
    return linear_index >= first_index && linear_index < last_index;
}

inline std::size_t to_chunk_index(const std::size_t index, const int rank, const std::size_t chunk_size)
{
    return index - rank * chunk_size;
}

void set_fixed_according_to_rank
(
    matrix_type& chunk, 
    const fixed_values_type& fixed_values, 
    const int rank, 
    const std::size_t chunk_size
)
{
    for (const auto& [x, y, value] : fixed_values)
    {
        if (is_index_within_rank(x, y, rank, chunk_size))
        {
            const std::size_t actual_index = x * N + y;
            const std::size_t index_in_chunk = to_chunk_index(actual_index, rank, chunk_size);

            chunk[index_in_chunk] = value;
        }
    }
}

inline value_type smooth_neighborhood
(
    const value_type i_j,
    const value_type li_j,
    const value_type i_lj,
    const value_type ri_j,
    const value_type i_rj
)
{
    return (i_j + li_j + i_lj + ri_j + i_rj) / 5;
}

void smooth_chunk
(
    matrix_type& chunk,
    const fixed_values_type& fixed_values,
    const std::size_t iteration_count
)
{
    const std::size_t m_chunk = chunk.size() / N;
    const std::size_t n = N;

    int number_of_processors;
    int rank;

    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const bool has_upper_neighbor  = rank != 0;
    const bool has_lower_neighbor = rank != number_of_processors - 1;

    const int upper_neighbor_rank = rank - 1;
    const int lower_neighbor_rank  = rank + 1;

    matrix_type next_chunk(chunk.size());

    for (std::size_t iteration = 0; iteration < iteration_count; ++iteration)
    {
        row_type upper_row_received(n);
        row_type lower_row_received(n);

        if (has_upper_neighbor)
        {
            const value_type* upper_row_ptr = chunk.data();

            MPI_Send(upper_row_ptr, n, MPI_FLOAT, upper_neighbor_rank, 0, MPI_COMM_WORLD);
            MPI_Recv(upper_row_received.data(), n, MPI_FLOAT, upper_neighbor_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (has_lower_neighbor)
        {
            const value_type* lower_row_ptr = chunk.data() + chunk.size() - n;

            MPI_Send(lower_row_ptr, n, MPI_FLOAT, lower_neighbor_rank, 0, MPI_COMM_WORLD);
            MPI_Recv(lower_row_received.data(), n, MPI_FLOAT, lower_neighbor_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (has_upper_neighbor)
        {
            for (std::size_t j = 1u; j < n - 1; ++j)
            {
                next_chunk[0u * n + j] = smooth_neighborhood
                (
                    chunk[0u * n +        j],
                    upper_row_received[j],
                    chunk[0u * n + (j - 1u)],
                    chunk[1u * n +        j],
                    chunk[0u * n + (j + 1u)]
                );
            }
        }
        
        if (has_lower_neighbor)
        {
            for (std::size_t j = 1u; j < n - 1; ++j)
            {
                next_chunk[(m_chunk - 1) * n + j] = smooth_neighborhood
                (
                    chunk[(m_chunk - 1u) * n +        j],
                    chunk[(m_chunk - 2u) * n +        j],
                    chunk[(m_chunk - 1u) * n + (j - 1u)],
                    upper_row_received[j],
                    chunk[(m_chunk - 1u) * n + (j + 1u)]
                );
            }
        }

        for (std::size_t i = 1u; i < m_chunk - 1; ++i)
        {
            for (std::size_t j = 1u; j < n; ++j)
            {
                next_chunk[i * n + j] = smooth_neighborhood
                (
                    chunk[i * n + j],
                    chunk[(i - 1u) * n + j],
                    chunk[i * n + (j - 1u)],
                    chunk[(i + 1u) * n + j],
                    chunk[i * n + (j + 1u)]
                );
            }
        }

        set_fixed_according_to_rank(next_chunk, fixed_values, rank, chunk.size());

        std::swap(next_chunk, chunk);
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

int main(int argc, char* argv[])
{
    std::cout << std::fixed << std::setprecision(3);

    if (auto err = MPI_Init(&argc, &argv))
        return err;

    int number_of_processors;
    int rank;

    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const bool is_main_process = rank == 0;

    if (argc != 4)
    {
        PRINTER << "Usage: " << argv[0] << " <iterations> <m> <n>" << '\n';
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    const std::size_t iterations = std::atoi(argv[1]);
    const std::size_t m = M = std::atoi(argv[2]);
    const std::size_t n = N = std::atoi(argv[3]);

    if (iterations <= 0 || m <= 0 || n <= 0)
    {
        PRINTER << "All arguments must be positive integers.\n";
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    if (number_of_processors >= m)
    {
        PRINTER << "Number of processors must be less then row count (m).";
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    PRINTER << "Number of processors: " << number_of_processors << '\n';
    PRINTER << "Iterations: " << iterations << '\n';
    PRINTER << "Matrix dimensions (m x n): " << m << " x " << n << '\n';

    const std::size_t total_rows_count = m % number_of_processors == 0
        ? m
        : m + (number_of_processors - m % number_of_processors);

    const std::size_t total_length = total_rows_count * n;
    const std::size_t chunk_size   = total_length / number_of_processors;

    PRINTER << "Chunk size: " << chunk_size << '\n';

    matrix_type matrix = get_zero_matrix(total_length);
    fixed_values_type fixed_values = get_default_fixed(m, n);

    set_fixed(matrix, fixed_values);

    if (m <= 20 && n <= 20)
    {
        PRINTER << "\nInitial matrix:\n";
        PRINT_MATRIX(matrix, m, n);
    }

    PRINTER << '\n';

    const double start_time = MPI_Wtime();

    matrix_type chunk(chunk_size);
    MPI_Scatter(matrix.data(), chunk_size, MPI_FLOAT, chunk.data(), chunk_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    smooth_chunk(chunk, fixed_values, iterations);

    matrix_type result(matrix.size());
    MPI_Gather(chunk.data(), chunk_size, MPI_FLOAT, result.data(), chunk_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    const double end_time = MPI_Wtime();

    if (m <= 20 && n <= 20)
    {
        PRINTER << "Result:\n";
        PRINT_MATRIX(result, m, n);
    }

    double elapsed_time = end_time - start_time;
    std::cout << std::setprecision(15);
    PRINTER << "Time taken: " << elapsed_time << "s." << std::endl;

    MPI_Finalize();
    return EXIT_SUCCESS;
}