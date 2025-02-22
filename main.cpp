#include <inttypes.h>
#include <locale.h>
#include <stdio.h>

#include <sys/mman.h>

#include <chrono>
#include <cstring>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "gmp.h"
#include "gmpxx.h"

#include "boost/iostreams/device/mapped_file.hpp"

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

struct CollatzAccelClass {
    uint16_t c; // c is the number of increases
    uint64_t f_k_b;
};

namespace {
const constexpr uint64_t MAX_K = 32;
const constexpr size_t ADVISE_RANDOM_ABOVE_K = 33;
const constexpr size_t MAX_K_FOR_CACHE = 32;

class CollatzSieve {
public:
    bool need_to_test(int_type test_num) {
        for (auto & [A, B] : sieve_list_) {
            if (((test_num - B) % A) == 0) {
                return false;
            }
        }

        return true;
    }

private:
    std::vector<std::pair<size_t, size_t>> sieve_list_ = {{9, 2}, {9, 4}, {9, 5}, {9, 8}, {8, 5}};
};

class CollatzAccelClasses {
public:
    CollatzAccelClasses(int k) :
    k_(k),
    f_k_b_len_(k_ < 40 ? 20 : 22)
    {
        std::filesystem::path base("/Volumes/tank/collatz");

        if (!is_directory(base)) {
            throw std::runtime_error("Base is not a directory!");
        }

        std::stringstream sieve_filename_ss;
        sieve_filename_ss << "sieve_2^" << k << ".txt";
        auto sieve_file_path = base / sieve_filename_ss.str();

        mfs_.open(sieve_file_path);
        if (k > ADVISE_RANDOM_ABOVE_K) {
            posix_madvise((void *)mfs_.data(), mfs_.size(), POSIX_MADV_RANDOM);
        } else if (k_ <= MAX_K_FOR_CACHE){
            posix_madvise((void *)mfs_.data(), mfs_.size(), POSIX_MADV_DONTNEED);
            const size_t cache_size = pow2(u64tompz(k_-1)).get_ui();
            cache_.reserve(cache_size);

            for (size_t i = 0; i < cache_size; ++i) {
                cache_.emplace_back(read_collatz_accel_class_from_file(i));
            }
            mfs_.close();
        } else {
            posix_madvise((void *)mfs_.data(), mfs_.size(), POSIX_MADV_WILLNEED);
        }
    }

    CollatzAccelClass operator[](size_t idx) const {
        if (k_ <= MAX_K_FOR_CACHE) {
            return cache_[idx];
        }

        return read_collatz_accel_class_from_file(idx);
    }

private:
    CollatzAccelClass read_collatz_accel_class_from_file(size_t idx) const {
        const off_t stride = 5 + 1 + f_k_b_len_ + 1;

        char buf[stride + 1]; memset(buf, '\0', sizeof(buf));
        uint16_t c;
        uint64_t f_k_b;
        memcpy(buf, mfs_.data() + (idx * stride), stride);

        if (k_ < 40) {
            if (sscanf(buf, legacy_printf_specifier_, &c, &f_k_b) == EOF) {
                throw std::runtime_error("legacy sscanf failed");
            }
        } else {
            if (scanf(buf, printf_specifier_, &c, &f_k_b) == EOF) {
                throw std::runtime_error("sscanf failed");
            }
        }

        return CollatzAccelClass{c, f_k_b};
    }

private:
    boost::iostreams::mapped_file_source mfs_;

    mutable std::vector<CollatzAccelClass> cache_;

    int k_;
    const char * legacy_printf_specifier_ = "%5hu %20llu";
    const char * printf_specifier_ = "%5hu %22llu";
    size_t f_k_b_len_;
};

class CollatzAccelClassesByK {
public:
    CollatzAccelClassesByK() {
        for (int k = 1; k <= MAX_K; ++k) {
            collatz_accel_classes_by_k_.emplace_back(k);
        }
    }

    const CollatzAccelClasses & operator[](size_t idx) const {
        return collatz_accel_classes_by_k_[idx - 1];
    }


private:
    std::vector<CollatzAccelClasses> collatz_accel_classes_by_k_;
};
} // anonymous namespace


// returns the collatz number after applying shortcut collatz function step number of times
int_type collatz(int_type test_num, int_type original_test_num, uint16_t steps, const CollatzAccelClassesByK & collatz_accel_classes_by_k) {
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
        return collatz(4, original_test_num, steps - 1, collatz_accel_classes_by_k);
    }

    if (steps > MAX_K) {
        const int_type res = collatz(test_num, original_test_num, MAX_K, collatz_accel_classes_by_k);
        return collatz(res, original_test_num, steps - MAX_K, collatz_accel_classes_by_k);
    }

    const int_type sieve_size = pow2(steps);
    const int_type A = test_num / sieve_size;
    const int_type B = test_num % sieve_size;

    const CollatzAccelClass & collatz_accel_class = collatz_accel_classes_by_k[steps][int_type(B / 2).get_ui()];
    uint16_t c = collatz_accel_class.c;
    const int_type f_k_b = u64tompz(collatz_accel_class.f_k_b);

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

    ParallelExecutorWithMaxPendingTasks<std::vector<NumPeakResult>> tp(48, 96);
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

    CollatzAccelClassesByK collatz_classes_by_k;
    const int_type batch_size = 1048576_mpz;
    const int_type max_test_num = 23589095998679804297590_mpz;
    for (int_type test_num = global_peak_result.num; test_num < max_test_num; test_num += (2 * batch_size)) {
        auto test_func = [&collatz_sieve, &collatz_classes_by_k, test_num, max_test_num, batch_size]() -> std::vector<NumPeakResult> {
            std::vector<NumPeakResult> peak_results;
            for (int_type batch_test_num = test_num; (batch_test_num < (test_num + 2 * batch_size)) && (batch_test_num < max_test_num); batch_test_num+=2) {
                NumPeakResult peak_result{batch_test_num,0,0};
                int_type res = batch_test_num;

                if (!collatz_sieve.need_to_test(batch_test_num)) {
                    continue;
                }

                for (uint16_t steps = 1; res >= test_num; ++steps) {
                    if (steps == std::numeric_limits<decltype(steps)>::max()) {
                        throw std::runtime_error("Steps would not fit in current representation");
                    }

                    res = collatz(batch_test_num, batch_test_num, steps, collatz_classes_by_k);

                    if (res > peak_result.peak) {
                        peak_result.peak = res;
                        peak_result.steps = steps;
                    }
                }

                peak_results.emplace_back(std::move(peak_result));
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
