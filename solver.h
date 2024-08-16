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
    cache_t cache;
    Eigen::Matrix<double, mu_dim, dim + dim> mult;
    auto beta_idx(Eigen::Index t_idx = 0) const {
        return Eigen::seqN(t_idx * info.C_hat.rows() + info.M_0.rows(), info.L_0.cols());
    }
    auto idx(Eigen::Index t_idx) const {
        return Eigen::seqN(t_idx * info.C_hat.rows(), info.C_hat.rows());
    }
public:
    solver() = default;
    template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
    solver(T1&& A, T2&& C, T3&& M_0, T4&& L_0, T5&& C_hat, T6&& Gamma, T7&& t);

    typedef Eigen::Matrix<double, mu_dim, Eigen::Dynamic> mu_matrix;

    std::tuple<double, double, Eigen::Matrix<double, dim, 1>>
    W(const Eigen::Matrix<double, dim, 1>& x, std::size_t T) const;

    template <std::random_access_iterator Iter>
    Eigen::Matrix<double, 2 * dim, Eigen::Dynamic>
    qint_LTI(const Eigen::Matrix<double, dim, 1>& x, const mu_matrix& mu_array, Iter T, Iter begin) const;
    template <std::random_access_iterator Iter>
    Eigen::Matrix<double, 2 * dim, Eigen::Dynamic>
    qint_LTI(const Eigen::Matrix<double, dim, 1>& x, const mu_matrix& mu_array, Iter T, const Iter begin, const Iter end) const;

    std::pair<Eigen::Matrix<double, 2 * dim, Eigen::Dynamic>, double> ode_init(const Eigen::Matrix<double, dim, 1>& x, std::size_t len) const;
    Eigen::Matrix<double, 2 * dim, 1>
    ode_rhs(const Eigen::Matrix<double, dim, 1>& x, Eigen::Matrix<double, 2 * dim, 1> state, double T) const;
    Eigen::Matrix<double, 2 * dim, Eigen::Dynamic>
    ode_next_euler(const Eigen::Matrix<double, dim, 1>& x, const Eigen::Matrix<double, 2 * dim, Eigen::Dynamic>& z_Pzq, Eigen::Index T) const;
    std::tuple<double, Eigen::Matrix<double, dim, 1>, Eigen::Matrix<double, dim, 1>>
    solve(const Eigen::Matrix<double, dim, 1>& x, std::span<const double> t, std::span<const std::size_t> N, std::span<const double> err) const;
    Eigen::Matrix<double, 2 * dim, Eigen::Dynamic>
    solve_times(const Eigen::Matrix<double, dim, 1>& x, std::size_t T) const;
};

