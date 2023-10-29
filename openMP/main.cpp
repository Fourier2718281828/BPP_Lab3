#include <iostream>
#include <iomanip>
#include <algorithm>
#include <array>
#include <vector>
#include <stdlib.h>
#include <chrono>
#include <omp.h>

template<typename T>
using row_type = std::vector<T>;

using value_type = float;
using matrix_type = row_type<row_type<value_type>>;
using fixed_values_type = std::vector<std::tuple<std::size_t, std::size_t, value_type>>;

matrix_type get_zero_matrix(const std::size_t m, const std::size_t n)
{
    return matrix_type (m, row_type<value_type>(n));
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

matrix_type smooth_matrix
(
    const matrix_type& matrix, 
    const fixed_values_type& fixed_values, 
    const std::size_t iterations_count
)
{
    matrix_type prev_iteration { matrix.cbegin(), matrix.cend() };
    matrix_type curr_iteration { matrix.cbegin(), matrix.cend() };
    const std::size_t m = matrix.size();
    const std::size_t n = matrix[0].size();

    for (std::size_t iteration = 0; iteration < iterations_count; ++iteration)
    {
        std::size_t i;

        #pragma omp parallel for
        for (i = 1ull; i < m - 1; ++i)
        {
            std::size_t j;

            #pragma omp parallel for
            for (j = 1ull; j < n - 1; ++j)
            {
                curr_iteration[i][j] = (
                    prev_iteration[i][j]     +
                    prev_iteration[i - 1][j] +
                    prev_iteration[i][j - 1] +
                    prev_iteration[i + 1][j] +
                    prev_iteration[i][j + 1]
                ) / 5;
            }
        }

        set_fixed(curr_iteration, fixed_values);
        std::swap(prev_iteration, curr_iteration);
    }

    return prev_iteration;
}

int main(int argc, const char* argv[])
{
    std::cout << std::fixed << std::setprecision(3);

    if (argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " <num_of_processors> <iterations> <m> <n>" << std::endl;
        return EXIT_FAILURE;
    }

    int num_of_processors = std::atoi(argv[1]);
    int iterations = std::atoi(argv[2]);
    int m = std::atoi(argv[3]);
    int n = std::atoi(argv[4]);

    if (num_of_processors <= 0 || iterations <= 0 || m <= 0 || n <= 0)
    {
        std::cerr << "All arguments must be positive integers.\n";
        return EXIT_FAILURE;
    }

    std::cout << "Number of processors: " << num_of_processors << '\n';
    std::cout << "Iterations: " << iterations << '\n';
    std::cout << "Matrix dimensions (m x n): " << m << " x " << n << '\n';

    omp_set_num_threads(num_of_processors);

    matrix_type matrix = get_zero_matrix(m, n);
    fixed_values_type fixed_values = get_default_fixed(m, n);

    set_fixed(matrix, fixed_values);

    if (m <= 20 && n <= 20)
    {
        std::cout << "\nInitial matrix:\n";
        print_matrix(matrix);
    }

    std::cout << '\n';

    
    auto start_time = std::chrono::high_resolution_clock::now();
    matrix_type result = smooth_matrix(matrix, fixed_values, iterations);
    auto end_time = std::chrono::high_resolution_clock::now();

    if (m <= 20 && n <= 20)
    {
        std::cout << "Result:\n";
        print_matrix(result);
    }

    std::cout << std::setprecision(15);

    std::chrono::duration<double> elapsed_time = end_time - start_time;
    std::cout << "\nTime taken: " << elapsed_time.count() << "s." << std::endl;

    return EXIT_SUCCESS;
}
