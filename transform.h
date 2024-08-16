#pragma once
#include "pde_info.h"
#include <cassert>
#include <vector>
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
    template <typename V>
    Eigen::Matrix<double, dim, V::ColsAtCompileTime> S(double t, V vec) const {
        return static_cast<const T*>(this)->S(t, vec);
    }
    auto S_begin() const {
        return static_cast<const T*>(this)->S_begin();
    }
    template <typename V>
    Eigen::Matrix<double, dim, V::ColsAtCompileTime> S_inv(double t, V vec) const {
        return static_cast<const T*>(this)->S_inv(t, vec);
    }
    template <typename V>
    Eigen::Matrix<double, dim, V::ColsAtCompileTime> S_inv_nearest(double t, V vec) const {
        return static_cast<const T*>(this)->S_inv_nearest(t, vec);
    }
};

class RS_interp : public RS_base<RS_interp> {
    const double init;
    const double len;
    std::size_t n;
    std::vector<Eigen::Matrix<double, dim, dim>> _R;
    std::vector<Eigen::Matrix<double, dim, dim>> _S;
    std::vector<Eigen::FullPivLU<Eigen::Matrix<double, dim, dim>>> _S_inv;
public:
    __attribute__((pure)) std::pair<std::ptrdiff_t, double> time_idx(double t) const noexcept {
        double i;
        double frac = std::modf((t - init) / len * n, &i);
        if (frac < 1e-14)
            frac = 0;
        return std::make_pair(static_cast<std::ptrdiff_t>(i), frac);
    }
    __attribute__((pure)) std::ptrdiff_t time_nearest(double t) const {
        return std::lround((t - init) / len * n);
    }
    RS_interp(const PDE_info& info, double init, double end, std::size_t n) noexcept : init(init), len(end - init), n(n) {
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
    __attribute__((pure)) Eigen::Matrix<double, dim, dim> R(double t) const noexcept {
        assert(t >= init && t - init <= len);
        auto [idx, frac] = time_idx(t);
        if (frac) [[unlikely]]
            return (1 - frac) * _R[idx] + frac * _R[idx + 1];
        else
            return _R[idx];
    }
    __attribute__((pure)) Eigen::Matrix<double, dim, dim> S(double t) const noexcept {
        assert(t >= init && t - init <= len);
        auto [idx, frac] = time_idx(t);
        if (frac) [[unlikely]]
            return (1 - frac) * _S[idx] + frac * _S[idx + 1];
        else
            return _S[idx];
    }
    const Eigen::Matrix<double, dim, dim>& S_begin() const {
        return _S[0];
    }
    const Eigen::Matrix<double, dim, dim>& R_end() const {
        return _R.back();
    }
    template <typename V>
    Eigen::Matrix<double, dim, V::ColsAtCompileTime> S(double t, V vec) const {
        return S(t) * vec;
    }
    template <typename V>
    Eigen::Matrix<double, dim, V::ColsAtCompileTime> S_inv(double t, V vec) const {
        assert(t >= init && t - init <= len);
        auto [idx, frac] = time_idx(t);
        return (1 - frac) * _S_inv[idx].solve(vec) + frac * _S_inv[idx + 1].solve(vec);
    }
    template <typename V>
    Eigen::Matrix<double, dim, V::ColsAtCompileTime> S_inv_nearest(double t, V vec) const {
        assert(t >= init && t - init <= len);
        return _S_inv[time_nearest(t)].solve(vec);
    }
    auto operator[](std::ptrdiff_t n) const {
        return std::make_pair(_R[n], _S[n]);
    }
};
