#pragma once
#include "pde_info.h"
#include <cassert>
#include <vector>
#include <functional>
#include <Eigen/LU>
#include <unsupported/Eigen/MatrixFunctions>

template <typename T>
struct RS_base {
    auto R(double t) const {
        return static_cast<const T*>(this)->R(t);
    }
    auto S(double t) const {
        return static_cast<const T*>(this)->S(t);
    }
    auto S_begin() const {
        return static_cast<const T*>(this)->S_begin();
    }
    auto R_end() const {
        return static_cast<const T*>(this)->R_end();
    }
    template <typename V>
    Eigen::Matrix<double, dim, V::ColsAtCompileTime> S_inv(double t, V vec) const {
        return static_cast<const T*>(this)->S_inv(t, vec);
    }
};

class RS_nearest : public RS_base<RS_nearest> {
    const double init;
    const double len;
    const double step_inv;
    std::vector<Eigen::Matrix<double, dim, dim>> _R;
    std::vector<Eigen::Matrix<double, dim, dim>> _S;
    std::vector<Eigen::PartialPivLU<Eigen::Matrix<double, dim, dim>>> _S_inv;
public:
    __attribute__((pure)) std::ptrdiff_t time_idx(double t) const noexcept {
        assert(t >= init && t - init <= len);
        return std::lround((t - init) * step_inv);
    }
    RS_nearest(const PDE_info& info, double init, double end, std::size_t n) noexcept : init(init), len(end - init), step_inv(n / len) {
        const double step = len / n;
        _R.reserve(n + 1);
        _S.reserve(n + 1);
        _S_inv.reserve(n + 1);
        Eigen::Matrix<double, dim, dim> S_end = (len * info.A_bar).exp().bottomRightCorner<dim, dim>();
        Eigen::Matrix<double, dim, dim> corr = S_end.colPivHouseholderQr().inverse();
        for (std::size_t i = 0; i <= n; ++i) {
            Eigen::Matrix<double, 2 * dim, dim> exp = (info.A_bar * (i * step)).exp().rightCols<dim>() * corr;
            _R.emplace_back(exp.topRows<dim>());
            _S.emplace_back(exp.bottomRows<dim>());
            _S_inv.emplace_back(_S.back());
        }
    }
    __attribute__((pure)) const Eigen::Matrix<double, dim, dim>& R(double t) const noexcept {
        return _R[time_idx(t)];
    }
    __attribute__((pure)) const Eigen::Matrix<double, dim, dim>& S(double t) const noexcept {
        return _S[time_idx(t)];
    }
    const Eigen::Matrix<double, dim, dim>& S_begin() const noexcept {
        return _S[0];
    }
    const Eigen::Matrix<double, dim, dim>& R_end() const noexcept {
        return _R.back();
    }
/*
    Eigen::Matrix<double, dim, 1> S_inv(double t, const Eigen::Matrix<double, dim, 1>& vec) const {
        return _S_inv[time_idx(t)].solve(vec);
    }
*/
    Eigen::Matrix<double, dim, 1> S_inv(std::size_t t, const Eigen::Matrix<double, dim, 1>& vec) const {
        return _S_inv[t].solve(vec);
    }
    auto operator[](std::ptrdiff_t n) const noexcept {
        return std::make_pair(std::cref(_R[n]), std::cref(_S[n]));
    }
};
