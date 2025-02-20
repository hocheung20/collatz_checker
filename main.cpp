#include <inttypes.h>
#include <locale.h>
#include <stdio.h>

#include <chrono>
#include <cstring>
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
    static const int_type one = 1_mpz;
    int_type res;
    mpz_mul_2exp(res.get_mpz_t(), one.get_mpz_t(), exp.get_ui());
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
    uint64_t f_k_b;
};

namespace {
const constexpr uint64_t max_k = 39;
const int_type max_k_mpz = u64tompz(max_k);
const int_type max_sieve_size = pow2(max_k_mpz);

class CollatzClasses {
public:
    CollatzClasses(int k) :
    k_(k),
    f_k_b_len_(k_ < 40 ? 20 : 22)
    {
//        std::filesystem::path base("/Volumes/tank/collatz");
        std::filesystem::path base("/mnt/tank/collatz");

        if (!is_directory(base)) {
            throw std::runtime_error("Base is not a directory!");
        }

        std::stringstream sieve_filename_ss;
        sieve_filename_ss << "sieve_2^" << k << ".txt";
        auto sieve_file_path = base / sieve_filename_ss.str();

        sieve_fptr_ = fopen(sieve_file_path.c_str(), "r");
        if (sieve_fptr_ == NULL) {
            throw std::runtime_error("Unable to open sieve file" + sieve_file_path.string());
        }
    }

    CollatzClasses(const CollatzClasses &) = delete;
    CollatzClasses& operator=(const CollatzClasses &) = delete;

    CollatzClasses(CollatzClasses && other) {
        operator=(std::move(other));
    }

    CollatzClasses& operator=(CollatzClasses && other) {
        std::swap(sieve_fptr_, other.sieve_fptr_);
        std::swap(k_, other.k_);
        std::swap(f_k_b_len_, other.f_k_b_len_);
        return *this;
    }

    ~CollatzClasses() {
        if (sieve_fptr_ != NULL) {
            fclose(sieve_fptr_);
        }

        sieve_fptr_ = NULL;
    }

    CollatzClass operator[](size_t idx) const {
        const off_t stride = 5 + 1 + f_k_b_len_ + 1;
        if (0 != fseek(sieve_fptr_, idx * stride, SEEK_SET)) {
            throw std::runtime_error("fseek failed");
        }

        char buf[stride + 1]; memset(buf, '\0', sizeof(buf));
        uint16_t c;
        uint64_t f_k_b;
        if (1 != fread(buf, stride, 1, sieve_fptr_)) {
            throw std::runtime_error("fread faild");
        }

        if (k_ < 40) {
            if (sscanf(buf, legacy_printf_specifier_, &c, &f_k_b) == EOF) {
                throw std::runtime_error("legacy sscanf failed");
            }
        } else {
            if (scanf(buf, printf_specifier_, &c, &f_k_b) == EOF) {
                throw std::runtime_error("sscanf failed");
            }
        }

        return CollatzClass{c, f_k_b};
    }

private:
    FILE * sieve_fptr_ = NULL;

    int k_;
    const char * legacy_printf_specifier_ = "%5hu %20llu";
    const char * printf_specifier_ = "%5hu %22llu";
    size_t f_k_b_len_;
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
        const int_type res = collatz(test_num, max_k_mpz, collatz_classes_by_k);
        return collatz(res, steps - max_k_mpz, collatz_classes_by_k);
    }

    const int_type sieve_size = pow2(steps);
    const int_type A = test_num / sieve_size;
    const int_type B = test_num % sieve_size;

    const CollatzClass & collatz_class = collatz_classes_by_k[steps.get_ui()][int_type(B / 2).get_ui()];
    uint16_t c = collatz_class.c;
    const int_type f_k_b = u64tompz(collatz_class.f_k_b);

    int_type res;
    mpz_ui_pow_ui(res.get_mpz_t(), 3, c);
    res = A * res;
    res += f_k_b;

    return res;
}

int main() {
    auto start = ClockType::now();

    struct NumPeakResult {
        int_type num;
        int_type peak;
        int_type steps;
    };

    ParallelExecutorWithMaxPendingTasks<std::vector<NumPeakResult>> tp(15, 32);
//    ParallelExecutorWithMaxPendingTasks<std::vector<NumPeakResult>> tp(1, 1);


//    NumPeakResult global_peak_result{2358909599867980429759_mpz,0_mpz,0};
//    uint64_t global_peak_result_number = 97;
    NumPeakResult global_peak_result{3_mpz,0_mpz,0};
    uint64_t global_peak_result_number = 2;
    std::thread check_peak_result_thread([&tp, &global_peak_result, &global_peak_result_number](){
        for (auto peak_results_opt = tp.get(); peak_results_opt != std::nullopt; peak_results_opt = tp.get()) {
            const auto & peak_results = *peak_results_opt;
            for (const auto & peak_result : peak_results) {
                if (peak_result.peak > global_peak_result.peak) {
                    std::cout << global_peak_result_number << ": collatz^" << peak_result.steps << "(" << peak_result.num << ") = " << peak_result.peak << std::endl;
                    global_peak_result = peak_result;
                    ++global_peak_result_number;
                }
            }
        }
    });

    const int_type batch_size = 1024_mpz;
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
                        }
                    }
                }

                peak_results.emplace_back(peak_result);
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
