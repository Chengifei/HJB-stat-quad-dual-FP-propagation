// Copyright (C) 2024 Y. Zheng
// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <cstddef>
#include <optional>
#include <cstdio>
#include <cmath>

struct iterative_solver_exc {
    char str[32];
    std::size_t i;
    double e;
};

struct fixed_point_no_check {
    void operator()(double) {}
};

class fixed_point_check {
    const char* str;
    std::optional<double> prev;
public:
    fixed_point_check(const char* str = "fxpt") : str(str) {}
    void operator()(std::size_t i, double diff) {
        if (!prev || (diff < *prev + 1)) [[likely]]
            prev = diff;
        else [[unlikely]] {
            iterative_solver_exc exc{ "", i, diff };
            std::snprintf(exc.str, sizeof exc.str, "%s %g", str, *prev);
            throw exc;
        }
    }
};

template <typename T, typename F, typename C = fixed_point_check>
std::pair<T, double> fixed_point(F&& f, T x, std::size_t fp_max, double rtol, double abstol, C&& c = {}) {
    T x_new = f(x);
    double norm_diff = (x - x_new).norm();
    const double norm_ref = std::max(rtol * norm_diff, abstol);
    for (std::size_t i = 0; i != fp_max && norm_diff > norm_ref; ++i) {
        c(i, norm_diff);
        x = std::move(x_new);
        x_new = f(x);
        norm_diff = (x - x_new).norm();
    }
    return std::make_pair(std::move(x_new), norm_diff);
}

template <typename T, typename F, typename J>
std::pair<T, double> newton(F&& f, J&& j, T x, std::size_t iter_max, double rtol, double abstol) {
    T res = f(x);
    double norm = res.norm();
    const double norm_ref = std::max(rtol * norm, abstol);
    for (std::size_t i = 0; i != iter_max && norm > norm_ref; ++i) {
        if (!std::isfinite(norm))
            throw iterative_solver_exc{ "newton", i, norm };
        x -= j(x).colPivHouseholderQr().solve(res);
        res = f(x);
        norm = res.norm();
    }
    return std::make_pair(std::move(x), norm);
}
