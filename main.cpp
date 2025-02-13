#include <inttypes.h>
#include <locale.h>
#include <stdio.h>

#include <chrono>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "gmp.h"
#include "gmpxx.h"

#include "ParallelExecutorWithMaxPendingTasks.hpp"

using ClockType = std::chrono::system_clock;
using int_type = mpz_class;

struct ModClass {
    int_type linear_coefficient;
    int_type constant;
};

bool is_odd(const int_type & k) {
    return mpz_odd_p(k.get_mpz_t());
}

bool is_even(const int_type & k) {
    return mpz_even_p(k.get_mpz_t());
}

int_type pow2(int_type exp) {
    int_type res;
    mpz_mul_2exp(res.get_mpz_t(), (1_mpz).get_mpz_t(), exp.get_ui());
    return res;
}

size_t log2(int_type k) {
    return mpz_sizeinbase(k.get_mpz_t(), 2) - 1;
}

int_type u64tompz(uint64_t u64)
{
    int_type z;
    mpz_import(z.get_mpz_t(), 1, -1, sizeof u64, 0, 0, &u64);
    return z;
}

struct CollatzClass {
    uint16_t c; // c is the number of increases
    int_type f_k_b;
};

namespace {
const constexpr uint64_t max_k = 36;
const int_type max_k_mpz = u64tompz(max_k);
const int_type max_sieve_size = pow2(max_k_mpz);

class CollatzClasses {
public:
    CollatzClasses(int k) : k_(k) {
//        std::filesystem::path base("/Volumes/tank/collatz");
        std::filesystem::path base("/Volumes/980Pro2TB");

        if (!is_directory(base)) {
            throw std::runtime_error("Base is not a directory!");
        }

        std::stringstream sieve_filename_ss;
        sieve_filename_ss << "sieve_2^" << k << ".txt";
        auto sieve_file_path = base / sieve_filename_ss.str();
        sieve_fptr_ = fopen(sieve_file_path.c_str(), "r");
        if (sieve_fptr_ == NULL) {
            throw std::runtime_error("Unable to open sieve file");
        }
    }

    CollatzClasses(const CollatzClasses &) = delete;
    CollatzClasses(CollatzClasses && other) {
        operator=(std::move(other));
    }

    CollatzClasses & operator=(const CollatzClasses &) = delete;
    CollatzClasses & operator=(CollatzClasses && other) {
        std::swap(k_, other.k_);
        std::swap(sieve_fptr_, other.sieve_fptr_);

        return *this;
    };

    ~CollatzClasses() {
        if (sieve_fptr_ != NULL) {
            fclose(sieve_fptr_);
        }
    }

    CollatzClass operator[](size_t idx) const {
        const off_t stride = 5 + 1 + 20 + 1; // c (5 bytes) + <space> + f_k_b (16 bytes) + "\n";
        if (fseek(sieve_fptr_, stride * idx, SEEK_SET) != 0) {
            throw std::runtime_error(strerror(errno));
        }
        uint16_t c = 0;
        char f_k_b[17];
        if (fscanf(sieve_fptr_, "%5" PRIu16 " %16s", &c, &f_k_b) == EOF) {
            throw std::runtime_error("fscanf returned EOF");
        }

        return CollatzClass{c, int_type(f_k_b, 10)};
    }

private:
    FILE * sieve_fptr_ = NULL;
    int k_;
};

class CollatzClassesByK {
public:
    CollatzClassesByK() {
        for (int k = 1; k <= max_k; ++k) {
            collatz_classes_by_k_.emplace_back(k);
        }
    }

    const CollatzClasses & operator[](size_t idx) const {
        return collatz_classes_by_k_[idx - 1];
    }


private:
    std::vector<CollatzClasses> collatz_classes_by_k_;
};
} // anonymous namespace


// returns the collatz number after applying shortcut collatz function step number of times
int_type collatz(int_type test_num, int_type steps, const CollatzClassesByK & collatz_classes_by_k) {
    if (steps == 0) {
        return test_num;
    }

    while (is_even(test_num)) {
        test_num /= 2;
        --steps;

        if (steps == 0) {
            return test_num;
        }
    }

    if (test_num == 1) {
        return collatz(4, steps - 1, collatz_classes_by_k);
    }

    if (steps > max_k_mpz) {
        int_type res = collatz(test_num, max_k_mpz, collatz_classes_by_k);
        return collatz(res, steps - max_k_mpz, collatz_classes_by_k);
    }

    //size_t num_size_in_bits = mpz_sizeinbase(test_num.get_mpz_t(), 2);
    int_type A, B;

    if (steps <= max_k_mpz) {
        const int_type sieve_size = pow2(steps);
        A = test_num / sieve_size;
        B = test_num % sieve_size;
    } else { // test_num is greater than our max sieve
        A = test_num / max_sieve_size;
        B = test_num - max_sieve_size;
    }

    const CollatzClass & collatz_class = collatz_classes_by_k[steps.get_ui()][int_type(B / 2).get_ui()];
    uint16_t c = collatz_class.c;
    int_type f_k_b = collatz_class.f_k_b;

    int_type res;
    mpz_ui_pow_ui(res.get_mpz_t(), 3, c);
    res = A * res;
    res += f_k_b;
//    std::cout << test_num << " " << sieve_k << " " << steps << " " << A << " " << B << " " << c << " " << f_k_b << " " << res << std::endl;

    return res;
}

int main() {
    auto start = ClockType::now();

    struct NumPeakResult {
        int_type num;
        int_type peak;
        int_type steps;
    };

    ParallelExecutorWithMaxPendingTasks<std::vector<NumPeakResult>> tp(std::thread::hardware_concurrency()-1, 32);

//    NumPeakResult global_peak_result{2358909599867980429759_mpz,494114680693632162318569565806744468898092_mpz,0};
    NumPeakResult global_peak_result{3_mpz,0_mpz,0};
    uint64_t global_peak_result_number = 2;
    std::thread check_peak_result_thread([&tp, &global_peak_result, &global_peak_result_number](){
        for (auto peak_results_opt = tp.get(); peak_results_opt != std::nullopt; peak_results_opt = tp.get()) {
            const auto & peak_results = *peak_results_opt;
            for (const auto & peak_result : peak_results) {
                if (peak_result.peak > global_peak_result.peak) {
                    if (peak_result.num > global_peak_result.num) {
                        ++global_peak_result_number;
                    }
                    std::cout << global_peak_result_number << ": collatz^" << peak_result.steps << "(" << peak_result.num << ") = " << peak_result.peak << std::endl;
                    global_peak_result = peak_result;
                }
            }
        }
    });

    const int_type batch_size = 32768_mpz;
    const int_type max_test_num = 23589095998679804297590_mpz;
    for (int_type test_num = global_peak_result.num; test_num < max_test_num; test_num += (2 * batch_size)) {

        auto test_func = [test_num, max_test_num, batch_size]() -> std::vector<NumPeakResult> {
            CollatzClassesByK collatz_classes_by_k;
            std::vector<NumPeakResult> peak_results;
            for (int_type batch_test_num = test_num; (batch_test_num < (test_num + 2 * batch_size)) && (batch_test_num < max_test_num); batch_test_num+=2) {
                NumPeakResult peak_result{batch_test_num,0,0};
                int_type res = batch_test_num;

                if (batch_test_num % 3 == 2) {
                    // Don't need check
                } else if (batch_test_num % 9 == 4) {
                    // Don't need check
                } else {
                    for (int_type steps = 1; res >= test_num; ++steps) {
                        res = collatz(batch_test_num, steps, collatz_classes_by_k);

                        if (res > peak_result.peak) {
                            peak_result.peak = res;
                            peak_result.steps = steps;
                            peak_results.emplace_back(peak_result);
                        }
                    }
                }
            }

            return peak_results;
        };
        tp.post(test_func);
    }
    tp.finished_posting();

    check_peak_result_thread.join();

    auto end = ClockType::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << " ms" << std::endl;
}
