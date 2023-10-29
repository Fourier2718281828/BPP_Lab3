#include <mpi.h>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <array>
#include <vector>
#include <stdlib.h>
#include <chrono>

#define OUTPUT_RESULT

using std::cout;

constexpr static std::size_t fixed_count = 3ull;
constexpr static std::size_t chunk_size = 3ull; 

using value_type = double;
using fixed_values = std::array<std::pair<std::size_t, value_type>, fixed_count>;
using vector_type = std::vector<value_type>;

void set_fixed(vector_type& vector, const fixed_values& fixed_values)
{
    for (std::size_t i = 0ull; i < fixed_values.size(); ++i)
    {
        const auto& [index, value] = fixed_values[i];
        vector[index] = value;
    }
}

int main2(int argc, char* argv[])
{
    int number_of_processors;
    int rank;

    if (auto err = MPI_Init(&argc, &argv))
        return EXIT_FAILURE;

    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const bool is_main_process = rank == 0;

    std::size_t m = 1ull;
    std::size_t n = 10ull;

    try
    {
        if (argc < 2 || argc > 3)
        {
            std::cerr << "Usage: MPI_Lab <iterations> [processors] [size]" << std::endl;
            return EXIT_FAILURE;
        }

        m = std::stoi(argv[1]);

        if (argc >= 3)
            n = std::stoi(argv[2]);
    }
    catch (const std::invalid_argument& e)
    {
        std::cerr << "Invalid argument: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    catch (const std::out_of_range& e)
    {
        std::cerr << "Out of range argument: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }


    if (rank == 0)
    {
        std::cout << "n = " << n << '\n';
        std::cout << "m = " << m << '\n';
    }

    fixed_values fixed_values{ { {0, 1.}, {n >> 1, 1.}, {n - 1, 1.} } };

    if (is_main_process)
    {
        std::cout << "Fixed elements (index, value): ";
        for (const auto& [index, value] : fixed_values)
        {
            std::cout << '{' << index << ',' << ' ' << value << '}' << ' ';
        }
        std::cout << '\n';
    }

    const int total_length = n % number_of_processors == 0
        ? n
        : (number_of_processors - n % number_of_processors) + n;
    
    const int chunk_size = total_length / number_of_processors;

    vector_type vector(total_length, 0.0);
    set_fixed(vector, fixed_values);

#ifdef OUTPUT_RESULT
    if (is_main_process)
    {
        std::cout << "Initial vector: ";
        for (std::size_t i = 0ull; i < n; ++i)
        {
            std::cout << vector[i] << ' ';
        }
        std::cout << '\n';
    }
#endif

    vector_type chunk(chunk_size, 0.0);
    vector_type next_chunk(chunk_size, 0.0);

    const bool is_left_edge = rank == 0;
    const bool is_right_edge = rank == number_of_processors - 1;
    const bool has_left_neighbor = !is_left_edge;
    const bool has_right_neighbor = !is_right_edge;

    const int left_neighbor_rank = rank - 1;
    const int right_neighbor_rank = rank + 1;

    const double start_time = MPI_Wtime();

    MPI_Scatter(vector.data(), chunk_size, MPI_DOUBLE, chunk.data(), chunk_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (std::size_t i = 0; i < m; ++i)
    {
        value_type left_edge_to_receive = 0.0;
        value_type right_edge_to_receive = 0.0;

        if (is_left_edge)
        {
            next_chunk[0] = 1;
        }

        if (is_right_edge)
        {
            next_chunk[chunk_size - 1ull] = 1;
        }

        if (has_left_neighbor)
        {
            MPI_Send(&chunk[0], 1, MPI_DOUBLE, left_neighbor_rank, 0, MPI_COMM_WORLD);
            MPI_Recv(&left_edge_to_receive, 1, MPI_DOUBLE, left_neighbor_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (has_right_neighbor)
        {
            MPI_Send(&chunk[chunk_size - 1ull], 1, MPI_DOUBLE, right_neighbor_rank, 0, MPI_COMM_WORLD);
            MPI_Recv(&right_edge_to_receive, 1, MPI_DOUBLE, right_neighbor_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (has_left_neighbor)
        {
            next_chunk[0] = chunk[0] == 1
                ? 1
                : (left_edge_to_receive + chunk[1]) / 2;
        }

        if (has_right_neighbor)
        {
            next_chunk[chunk_size - 1ull] = chunk[chunk_size - 1ull] == 1
                ? 1
                : (right_edge_to_receive + chunk[chunk_size - 2ull]) / 2;
        }

        for (int j = 1; j < chunk_size - 1; ++j)
        {
            if (chunk[j] == 1)
                next_chunk[j] = 1;
            else
                next_chunk[j] = (chunk[j - 1ull] + chunk[j + 1ull]) / 2;
        }

        std::swap(next_chunk, chunk);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    vector_type resulting_vector(vector.size(), 0.0);
    MPI_Gather(chunk.data(), chunk_size, MPI_DOUBLE, resulting_vector.data(), chunk_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    const double end_time = MPI_Wtime();

#ifdef OUTPUT_RESULT
    if (is_main_process)
    {
        std::cout << "Result: ";
        for (std::size_t i = 0ull; i < n; ++i)
        {
            std::cout << resulting_vector[i] << ' ';
        }
        std::cout << '\n';
    }
#endif

    if (is_main_process)
    {
        double elapsed_time = end_time - start_time;
        std::cout << "Time taken: " << elapsed_time << "s." << std::endl;
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}

