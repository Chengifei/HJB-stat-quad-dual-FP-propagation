#pragma once
#include "pde_info.h"
#include "ode.h"
#include <Eigen/Core>
#include <Eigen/QR>
#include <cstdio>
#include <optional>

class solver {
protected:
    PDE_info info;
    Eigen::Matrix<double, mu_dim, dim + dim> mult;
    auto beta_idx(Eigen::Index t_idx = 0) const {
        return Eigen::seqN(t_idx * info.C_hat.rows() + info.M_0.rows(), info.L_0.cols());
    }
    auto idx(Eigen::Index t_idx) const {
        return Eigen::seqN(t_idx * info.C_hat.rows(), info.C_hat.rows());
    }
public:
    solver() = default;
    template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
    solver(T1&& A, T2&& C, T3&& M_0, T4&& L_0, T5&& C_hat, T6&& Gamma);
    std::tuple<double, Eigen::Matrix<double, dim, 1>, Eigen::Matrix<double, dim, 1>>
    solve(const Eigen::Matrix<double, dim, 1>& x, std::span<const double> t, std::span<const std::size_t> N, std::span<const double> err) const;
};

