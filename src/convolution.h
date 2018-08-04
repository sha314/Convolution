//
// Created by shahnoor on 2/1/2018.
//

#ifndef CONVOLUTION_CONVOLUTION_H
#define CONVOLUTION_CONVOLUTION_H

#include <vector>
#include <cstddef>

std::vector<double> convolve_v1(const std::vector<double>& data);
std::vector<double> convolve_v2(const std::vector<double>& data);
std::vector<double> convolve_v3(const std::vector<double>& data);

std::vector<std::vector<double>> convolve_multi_v1(const std::vector<std::vector<double>>& data);
std::vector<std::vector<double>> convolve_multi_v2(const std::vector<std::vector<double>>& data);
// 7 July, 2018
std::vector<std::vector<double>> convolve_multi_threaded_v1(const std::vector<std::vector<double>> &data_in);
// 8 July, 2018
std::vector<std::vector<double>> convolve_multi_threaded_v2(const std::vector<std::vector<double>> &data_in);

// 7 July, 2018
void calculate_convolution_for_rows(
        const size_t start,
        const size_t end,
        const std::vector<std::vector<double>> &data,
        std::vector<std::vector<double>> &out_data
);

// 8 july, 2018
void calculate_convolution_for_rows_v2(
        const size_t start,
        const size_t end,
        const std::vector<std::vector<double>> &data,
        std::vector<std::vector<double>> &out_data,
        const size_t id,
        std::vector<size_t> &iterations
);


/**
 * A Class to make using convolution user friendly
 */
class Convolution{
    std::vector<double> _forward_factor;
    std::vector<double> _backward_factor;
    size_t _number_of_data{};
    bool _initialized{false};
    double _count{}; // keep record of progress percentage. very useful for multithreaded environment ???
    double _time_elapsed_initialization{};
    double _time_elapsed_convolution{};
public:
    ~Convolution() = default;
    Convolution() = default;
    explicit Convolution(size_t n);

    // array of values
    // single column version
    std::vector<double> run(std::vector<double>& data_in); // single column and single threaded version
    std::vector<double> run_omp(std::vector<double>& data_in);
    std::vector<double> run_acc(std::vector<double>& data_in);
    std::vector<double> run_pthread(std::vector<double> &data_in);

    // array of columns of values
    // multiple column version
    std::vector<std::vector<double>> run_multi(std::vector<std::vector<double>>& data_in);
    std::vector<std::vector<double>> run_multi_omp(std::vector<std::vector<double>>& data_in);
    std::vector<std::vector<double>> run_multi_pthread(std::vector<std::vector<double>>& data_in);

    void timeElapsed() const {
        std::cout << "Initialization time " << _time_elapsed_initialization << " sec" << std::endl;
        std::cout << "Convolution time " << _time_elapsed_convolution << " sec" << std::endl;
    }

    double progress() {return _count/_number_of_data;}
private:
    void initialize(size_t n) ;

    void convolution_single_range(long row_start, long row_stop, const std::vector<double> &data_in,
                                  std::vector<double> &data_out);

    void convolution_multi_range(long row_start, long row_stop, const std::vector<std::vector<double>>  &data_in,
                                 std::vector<std::vector<double>> &data_out);
};


#endif //CONVOLUTION_CONVOLUTION_H
