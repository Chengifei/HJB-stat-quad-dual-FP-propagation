// Copyright (C) 2024 Y. Zheng
// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <vector>
#include <array>
#include <memory>
#include <iterator>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/LU>
#include <span>
#include "solvers.h"

constexpr int dim = 25;
constexpr int M_dim = 1;
constexpr int L_dim = 1;
constexpr int mu_dim = M_dim + L_dim;

namespace nonlinearity {
    template <typename T>
    concept alpha_sized = T::RowsAtCompileTime == M_dim;

    template <alpha_sized T>
    __attribute__((pure)) inline auto l(const Eigen::MatrixBase<T>& M0x) noexcept {
        return M0x.array().cos().matrix();
    }
    template <alpha_sized T>
    __attribute__((pure)) auto l_der(const Eigen::MatrixBase<T>& M0x) noexcept {
        return -M0x.array().sin().matrix();
    }
    __attribute__((pure)) inline Eigen::Matrix<double, M_dim, M_dim>
    l_dder(const Eigen::Matrix<double, M_dim, 1>& M0x) noexcept {
        Eigen::Matrix<double, M_dim, M_dim> ret;
        ret(0, 0) = -std::cos(M0x[0]);
        return ret;
    }

    template <alpha_sized T>
    __attribute__((pure)) auto f(const Eigen::MatrixBase<T>& M0x) noexcept {
        return M0x.array().atan().matrix();
    }
    template <alpha_sized T>
    __attribute__((pure)) auto f_der(const Eigen::MatrixBase<T>& M0x) noexcept {
        auto sqr = M0x.array().abs2();
        return (1 / (sqr + 1)).matrix();
    }
    __attribute__((pure)) inline std::array<Eigen::Matrix<double, M_dim, M_dim>, L_dim>
    f_dder(const Eigen::Matrix<double, M_dim, 1>& M0x) noexcept {
        std::array<Eigen::Matrix<double, M_dim, M_dim>, L_dim> ret;
        const double tmp = M0x[0] * M0x[0] + 1;
        ret[0] = -2 * M0x / (tmp * tmp);
        return ret;
    }
}

class Ntilde {
    Eigen::Matrix<double, mu_dim, mu_dim> C_hat;
    Eigen::FullPivLU<Eigen::Matrix<double, mu_dim, mu_dim>> C_inv;
public:
    Ntilde() = default;
    Ntilde(const Eigen::Matrix<double, mu_dim, mu_dim>& C_hat) : C_hat(C_hat), C_inv(C_hat) {}
    template <typename T>
    static Eigen::Matrix<double, 1, T::ColsAtCompileTime> eval(const Eigen::PlainObjectBase<T>& args);
    Eigen::RowVectorXd operator()(const Eigen::Matrix<double, mu_dim, Eigen::Dynamic>& args) const {
        return eval(args);
    }
    template <typename T>
    static Eigen::Matrix<double, mu_dim, T::ColsAtCompileTime> gradient(const Eigen::PlainObjectBase<T>& args);
    __attribute__((pure)) static Eigen::Matrix<double, mu_dim, mu_dim> dder(const Eigen::Matrix<double, mu_dim, 1>& args);
    template <typename M>
    std::pair<Eigen::Matrix<double, mu_dim, M::ColsAtCompileTime>, double> dual_argstat(const Eigen::MatrixBase<M>& y) const {
        // find argstat_x \{ operator()(x) - (x - y)' (this->C_hat/2) (x - y) \}
        // given y, solve gradient(x) = this->C_hat (x - y)
        // CAUTION: this->C_hat is initialized to -info.C_hat
        typedef Eigen::Matrix<double, mu_dim, M::ColsAtCompileTime> mu_matrix;
        mu_matrix ret(y.rows(), y.cols());
        Eigen::Matrix<double, 1, 1> C1_inv = C_inv.inverse().topLeftCorner(1, 1);
        Eigen::Matrix<double, L_dim, L_dim> C2_inv = C_inv.inverse().bottomRightCorner<L_dim, L_dim>();
        ret.template topRows<M_dim>() = fixed_point<Eigen::Matrix<double, M_dim, M::ColsAtCompileTime>>(
            [&](const Eigen::Matrix<double, M_dim, M::ColsAtCompileTime>& x) -> Eigen::Matrix<double, M_dim, M::ColsAtCompileTime> {
            mu_matrix xp;
            xp.template topRows<M_dim>() = x;
            xp.template bottomRows<L_dim>() = y.row(1) + C2_inv * nonlinearity::f(x.row(0));
            return C1_inv * gradient(x).template topRows<M_dim>() + y.template topRows<M_dim>();
            }, y.template topRows<M_dim>(), 200, 1e-4, 1e-4).first;
        ret.template bottomRows<L_dim>() = y.row(1) + C2_inv * nonlinearity::f(ret.template topRows<M_dim>());
        Eigen::Index col;
        double res = 0;
        for (Eigen::Index i = 0; i != y.cols(); ++i) {
            auto [sol, _res] = newton([&](const Eigen::Matrix<double, mu_dim, 1>& x) -> Eigen::Matrix<double, mu_dim, 1> {
                return C_hat * (x - y.col(i)) - gradient(x);
                },
                [&](const Eigen::Matrix<double, mu_dim, 1>& mu) -> Eigen::Matrix<double, mu_dim, mu_dim> {
                    return C_hat - dder(mu);
                }, ret.col(i).eval(), 25, 1e-10, 1e-8 * y.rows());
            ret.col(i) = sol;
            if (_res > res) {
                res = _res;
                col = i;
            }
        }
        return std::make_pair(std::move(ret), res);
    }
    __attribute__((pure)) Eigen::Matrix<double, mu_dim, mu_dim> dual_argstat_der(const Eigen::Matrix<double, mu_dim, 1>& mu) const {
        Eigen::Matrix<double, mu_dim, mu_dim> ret;
        ret.setIdentity();
        ret -= C_inv.solve(dder(dual_argstat(mu).first));
        return ret.fullPivLu().inverse();
    }
};

class uniform_linear_interpolator {
    double del;
    double t0;
    double t1;
public:
    template <typename Iter>
    uniform_linear_interpolator(Iter begin, Iter end) noexcept
        : t0(*begin), del(*(begin + 1) - *begin), t1(*(end - 1)) {}
    template <typename T, std::input_iterator It>
    T operator()(double a, It begin) const noexcept {
        assert(a >= t0 && a <= t1);
        double idx = (a - t0) / del;
        double i;
        double r = std::modf(idx, &i);
        begin += static_cast<It::difference_type>(i);
        if (r > 0x1p-20)
            return (1 - r) * (*begin) + r * (*(begin + 1));
        else
            return *begin;
    }
};

class linear_interpolator {
    std::span<const double> t;
public:
    template <typename Iter>
    linear_interpolator(Iter begin, Iter end) : t(begin, end) {
        assert(end != begin);
        for (auto it = begin + 1; it != end; ++it)
            assert(*(it - 1) < *it);
    }
    template <typename T, std::input_iterator It>
    T operator()(double a, It begin) const {
        auto it = std::upper_bound(t.begin(), t.end(), a);
        if (it == t.end())
            return *(begin + (t.size() - 1));
        else if (it == t.begin())
            return *begin;
        else {
            const double len = *it - *(it - 1);
            std::advance(begin, it - t.begin());
            if (a - *(it - 1) < 0x1p-20 * len)
                return *(begin - 1);
            return ((a - *(it - 1)) * (*begin) + (*it - a) * (*(begin - 1))) / len; // overflow?
        }
    }
};

struct PDE_info {
    Eigen::Matrix<double, dim, dim> A;
    Eigen::SparseMatrix<double> A_sparse;
    Eigen::Matrix<double, dim, dim> C;
    Eigen::Matrix<double, Eigen::Dynamic, dim, Eigen::AutoAlign, dim, dim> M_0;
    Eigen::Matrix<double, dim, Eigen::Dynamic, Eigen::AutoAlign, dim, dim> L_0;
    Eigen::Matrix<double, mu_dim, mu_dim> C_hat;
    Eigen::FullPivLU<Eigen::Matrix<double, mu_dim, mu_dim>> C_inv;
    Ntilde _;
    Eigen::Matrix<double, dim, dim> Gamma;
    std::vector<double> t;
    Eigen::Matrix<double, 2 * dim, 2 * dim> A_bar;
    template <typename T, std::input_iterator It>
    T interp(double a, It begin) const {
        linear_interpolator li{t.cbegin(), t.cend()};
        return li.template operator()<T, It>(a, begin);
    }
    template <typename T> requires std::derived_from<T, Eigen::PlainObjectBase<T>>
    Eigen::Matrix<double, mu_dim, T::ColsAtCompileTime> eta_inv(const Eigen::DenseBase<T>& eta) const {
        return eta.derived() + C_inv.solve(_.gradient(eta.derived()));
    }
    template <typename T> requires std::derived_from<T, Eigen::PlainObjectBase<T>>
    Eigen::Matrix<double, mu_dim, T::ColsAtCompileTime> C_hat_eta_inv(const Eigen::DenseBase<T>& eta) const {
        return C_hat * eta.derived() + _.gradient(eta.derived());
    }
    template <typename T>
    Eigen::Matrix<double, mu_dim, T::ColsAtCompileTime>
    eta_inv(const Eigen::DenseBase<T>& eta) const {
        auto _eta = eta.derived().matrix().eval();
        return eta_inv(_eta);
    }
    Eigen::Matrix<double, mu_dim, mu_dim> C_hat_jacobian_eta_inv(const Eigen::Matrix<double, mu_dim, 1>& eta) const {
        // C_hat * _.dual_argstat_der(mu).inverse()
        return C_hat + _.dder(eta);
    }
};
