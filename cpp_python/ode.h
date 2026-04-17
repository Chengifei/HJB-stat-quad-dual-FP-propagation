// Copyright (C) 2024 Y. Zheng
// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra_dispatcher.hpp>
#include "pde_info.h"
#include "transform.h"
typedef Eigen::Matrix<double, 2 * dim, 1> state_t;

namespace boost::numeric::odeint::detail {
    template<int Rows, int Cols>
    struct is_range<Eigen::Matrix<double, Rows, Cols>> :
    std::integral_constant<bool, (Rows < 0 || Cols < 0)>
    {};
}

struct no_control {
    const PDE_info& info;
    no_control(const PDE_info& info) noexcept : info(info) {}
    template <typename T>
    void operator()(const T& in, T& out, double) const {
        out.noalias() = info.A_bar * in;
    }
};

struct forward {
    const PDE_info& info;
    const Eigen::Matrix<double, mu_dim, 1> mu;
    forward(const PDE_info& info, Eigen::Matrix<double, mu_dim, 1> C_mu) noexcept
        : info(info), mu(C_mu) {}
    void operator()(const state_t& in, state_t& out, double) const {
        const Eigen::Index _dim = info.A.rows();
        out.noalias() = info.A_bar * in;
        out.template topRows<dim>(_dim) += info.L_0 * mu.bottomRows(info.L_0.cols());
        out.template bottomRows<dim>(_dim) -= info.M_0.transpose() * mu.topRows(info.M_0.rows());
    }
};

struct closed_loop {
    const PDE_info& info;
    closed_loop(const PDE_info& info) noexcept : info(info) {}
    template <typename T, typename U>
    void operator()(const T& in, U&& out, double) const {
        const Eigen::Index _dim = info.A.rows();
        Eigen::Matrix<double, mu_dim, 1> eta;
        eta.topRows(info.M_0.rows()) = info.M_0 * in.template topRows<dim>(_dim);
        eta.bottomRows(info.L_0.cols()) = info.L_0.transpose() * in.template bottomRows<dim>(_dim);
        const Eigen::Matrix<double, mu_dim, 1> C_mu = info.C_hat * info.eta_inv(eta);
        out.noalias() = info.A_bar * in;
        out.template topRows<dim>(_dim) += info.L_0 * C_mu.bottomRows(info.L_0.cols());
        out.template bottomRows<dim>(_dim) -= info.M_0.transpose() * C_mu.topRows(info.M_0.rows());
    }
};

struct closed_loop_with_r {
    const PDE_info& info;
    closed_loop_with_r(const PDE_info& info) noexcept : info(info) {}
    typedef Eigen::Matrix<double, 2 * dim + 1, 1> state_t;
    void operator()(const state_t& in, state_t& out, double) const {
        const Eigen::Index _dim = info.A.rows();
        auto z = in.topRows<dim>(_dim);
        auto lambda = in.middleRows<dim>(_dim, _dim);
        Eigen::Matrix<double, mu_dim, 1> eta;
        eta.topRows(info.M_0.rows()) = info.M_0 * z;
        eta.bottomRows(info.L_0.cols()) = info.L_0.transpose() * lambda;
        double N = info._.eval(eta)[0];
        const Eigen::Matrix<double, mu_dim, 1> C_mu = info.C_hat * info.eta_inv(eta);
        out.template topRows<dim + dim>(_dim + _dim).noalias() = info.A_bar * in.topRows<dim + dim>(_dim + _dim);
        out.template topRows<dim>(_dim) += info.L_0 * C_mu.bottomRows(info.L_0.cols());
        out.template middleRows<dim>(_dim, _dim) -= info.M_0.transpose() * C_mu.topRows(info.M_0.rows());
        eta.topRows(info.M_0.rows()).setZero();
        out[_dim + _dim] = 0.5 * (z.dot(info.C * z) + lambda.dot(info.Gamma * lambda) + eta.dot(info.C_hat * eta)) + N - C_mu.dot(eta);
    }
};

class closed_loop_with_der : closed_loop {
    Eigen::Matrix<double, 2 * dim, mu_dim> mult;
    Eigen::Matrix<double, mu_dim, 2 * dim> mult2;
public:
    closed_loop_with_der(const PDE_info& info) noexcept : closed_loop(info) {
        mult.setZero();
        mult.topRightCorner<dim, Eigen::Dynamic>(info.L_0.rows(), info.L_0.cols()) = info.L_0;
        mult.bottomLeftCorner<dim, Eigen::Dynamic>(info.M_0.cols(), info.M_0.rows()) = -info.M_0.transpose();
        mult2.setZero();
        mult2.topLeftCorner(info.M_0.rows(), dim) = info.M_0;
        mult2.bottomRightCorner(info.L_0.cols(), dim) = info.L_0.transpose();
    }
    template <int N>
    void operator()(const Eigen::Matrix<double, 2 * dim, N>& in, Eigen::Matrix<double, 2 * dim, N>& out, double t) const {
        const Eigen::Index _dim = info.A.rows();
        closed_loop::operator()(in.col(0), out.col(0), t);
        out.template rightCols<N - 1>(N - 1).noalias() = (info.A_bar + mult * info.C_hat_jacobian_eta_inv(mult2 * in.col(0)) * mult2) * in.template rightCols<N - 1>(N - 1);
    }
};

struct backward {
    const PDE_info& info;
    const Eigen::Matrix<double, mu_dim, Eigen::Dynamic>& mu_array;
    backward(const PDE_info& info, const Eigen::Matrix<double, mu_dim, Eigen::Dynamic>& mu_array) noexcept
        : info(info), mu_array(mu_array) {}
    void operator()(const state_t& in, state_t& out, double t) const {
        const Eigen::Index _dim = info.A.rows();
        Eigen::Matrix<double, mu_dim, 1> C_mu = info.C_hat * info.interp<Eigen::Matrix<double, mu_dim, 1>>(t, mu_array.colwise().cbegin()); // CAUTION already scaled!
        out.noalias() = info.A_bar * in;
        out.topRows<dim>(_dim) += info.L_0 * C_mu.bottomRows(info.L_0.cols());
        out.bottomRows<dim>(_dim) -= info.M_0.transpose() * C_mu.topRows(info.M_0.rows());
    }
};

struct RS_backward {
    const PDE_info& info;
    const RS_nearest& RS;
    const Eigen::Matrix<double, mu_dim, Eigen::Dynamic>& mu_array;
    Eigen::Matrix<double, dim, dim> Sc1;
    RS_backward(const PDE_info& info, const RS_nearest& RS, const Eigen::Matrix<double, mu_dim, Eigen::Dynamic>& mu_array) noexcept
        : info(info), RS(RS), mu_array(mu_array) {
        Sc1 = info.C - info.M_0.transpose() * info.C_hat.topLeftCorner<M_dim, M_dim>() * info.M_0;
    }
    void operator()(const state_t& in, state_t& out, double t) const noexcept {
        const Eigen::Index _dim = info.A.rows();
        Eigen::Matrix<double, mu_dim, 1> C_mu = info.C_hat * mu_array.col(RS.time_idx(t));
        auto p = in.topRows<dim>(_dim);
        Eigen::Matrix<double, dim, 1> tmp;
        tmp.noalias() = Sc1 * p + (info.M_0.transpose() * C_mu.topRows<M_dim>(info.M_0.rows()));
        out.bottomRows<dim>(_dim) = -RS.S_inv(t, tmp);
        out.topRows<dim>(_dim).noalias() = info.A * p - RS.R(t) * out.bottomRows<dim>(_dim) + info.L_0 * C_mu.bottomRows<L_dim>(info.L_0.cols());
    }
    state_t operator()(std::size_t t, const Eigen::Ref<const state_t>& in) const noexcept {
        const Eigen::Index _dim = info.A.rows();
        Eigen::Matrix<double, mu_dim, 1> C_mu = info.C_hat * mu_array.col(t);
        auto p = in.topRows<dim>(_dim);
        Eigen::Matrix<double, dim, 1> tmp;
        tmp = Sc1 * p + (info.M_0.transpose() * C_mu.topRows<M_dim>(info.M_0.rows()));
        state_t out;
        out.bottomRows<dim>(_dim) = -RS.S_inv(t, tmp);
        out.topRows<dim>(_dim) = info.A * p - RS[t].first * out.bottomRows<dim>(_dim) + (info.L_0 * C_mu.bottomRows<L_dim>(info.L_0.cols()));
        return out;
    }
};

struct RS_closed_loop {
    const PDE_info& info;
    const RS_nearest& RS;
    Eigen::Matrix<double, dim, dim> Sc1;
    RS_closed_loop(const PDE_info& info, const RS_nearest& RS) noexcept
        : info(info), RS(RS) {
        Sc1 = info.C - info.M_0.transpose() * info.C_hat.topLeftCorner<M_dim, M_dim>() * info.M_0;
    }
    void operator()(const state_t& in, state_t& out, double t) const noexcept {
        const Eigen::Index _dim = info.A.rows();
        auto p = in.topRows<dim>(_dim);
        auto q = in.bottomRows<dim>(_dim);
        Eigen::Matrix<double, mu_dim, 1> eta;
        eta.topRows<M_dim>() = info.M_0 * (p + RS.R(t) * q);
        eta.bottomRows<L_dim>() = info.L_0.transpose() * (RS.S(t) * q);
        Eigen::Matrix<double, mu_dim, 1> C_mu = info.C_hat_eta_inv(eta);
        Eigen::Matrix<double, dim, 1> tmp;
        tmp.noalias() = Sc1 * p + (info.M_0.transpose() * C_mu.topRows<M_dim>(info.M_0.rows()));
        out.bottomRows<dim>(_dim) = -RS.S_inv(t, tmp);
        out.topRows<dim>(_dim).noalias() = info.A * p - RS.R(t) * out.bottomRows<dim>(_dim) + info.L_0 * C_mu.bottomRows<L_dim>(info.L_0.cols());
    }
    state_t operator()(std::size_t t, const Eigen::Ref<const state_t>& in) const noexcept {
        const Eigen::Index _dim = info.A.rows();
        auto p = in.topRows<dim>(_dim);
        auto q = in.bottomRows<dim>(_dim);
        Eigen::Matrix<double, mu_dim, 1> eta;
        eta.topRows<M_dim>() = info.M_0 * (p + RS[t].first * q);
        eta.bottomRows<L_dim>() = info.L_0.transpose() * (RS[t].second * q);
        Eigen::Matrix<double, mu_dim, 1> C_mu = info.C_hat_eta_inv(eta);
        Eigen::Matrix<double, dim, 1> tmp;
        tmp = Sc1 * p + (info.M_0.transpose() * C_mu.topRows<M_dim>(info.M_0.rows()));
        state_t out;
        out.bottomRows<dim>(_dim) = -RS.S_inv(t, tmp);
        out.topRows<dim>(_dim) = info.A * p - RS[t].first * out.bottomRows<dim>(_dim) + (info.L_0 * C_mu.bottomRows<L_dim>(info.L_0.cols()));
        return out;
    }
};

template <typename _state_type = state_t>
struct pre_discretized_rk4_stepper {
    typedef _state_type state_type;
    typedef state_type deriv_type;
    typedef double time_type;
    typedef double value_type;
    typedef std::uint_least8_t order_type;
    typedef boost::numeric::odeint::stepper_tag stepper_category;
    static constexpr order_type order() {
        return 4;
    }
private:
    Eigen::Matrix<double, state_type::RowsAtCompileTime, 4, Eigen::RowMajor> k;
    alignas(64) static constexpr double _coeffs[4] = { 1. / 3., 2. / 3., 2. / 3., 1. / 3. };
    const double half_h;
    std::size_t counter = 0;
public:
    pre_discretized_rk4_stepper(double h) noexcept : half_h(0.5 * h) {}
    void reset() noexcept {
        counter = 0;
    }
    template <typename T>
    void do_step(const T& sys, Eigen::Ref<state_type> x, time_type, time_type) noexcept {
        k.col(0) = sys(counter, x);
        ++counter;
        k.col(1) = sys(counter, x + half_h * k.col(0));
        k.col(2) = sys(counter, x + half_h * k.col(1));
        ++counter;
        k.col(3) = sys(counter, x + (2 * half_h) * k.col(2));
        Eigen::Map<const Eigen::Matrix<double, 4, 1>, Eigen::Aligned64> rk_coeffs(_coeffs);
        x += k * (half_h * rk_coeffs);
    }
};
