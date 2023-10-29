#include <mpi.h>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <array>
#include <vector>
#include <stdlib.h>

using std::cout;

#define PRINTER if (is_main_process) std::cout 

template<typename T>
using row_type = std::vector<T>;

using value_type = float;
using matrix_type = row_type<row_type<value_type>>;
using fixed_values_type = std::vector<std::tuple<std::size_t, std::size_t, value_type>>;

matrix_type get_zero_matrix(const std::size_t m, const std::size_t n)
{
    return matrix_type(m, row_type<value_type>(n));
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

void print_matrix(const matrix_type& matrix)
{
    for (const auto& row : matrix)
    {
        for (auto&& element : row)
        {
            std::cout << element << ' ';
        }
        std::cout << '\n';
    }
}

void set_fixed(matrix_type& matrix, const fixed_values_type& fixed_values)
{
    for (const auto& [x, y, value] : fixed_values)
    {
        matrix[x][y] = value;
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

    int iterations = std::atoi(argv[1]);
    int m = std::atoi(argv[2]);
    int n = std::atoi(argv[3]);

    if (iterations <= 0 || m <= 0 || n <= 0)
    {
        PRINTER << "All arguments must be positive integers.\n";
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    PRINTER << "Number of processors: " << number_of_processors << '\n';
    PRINTER << "Iterations: " << iterations << '\n';
    PRINTER << "Matrix dimensions (m x n): " << m << " x " << n << '\n';

    MPI_Finalize();
    return EXIT_SUCCESS;
}