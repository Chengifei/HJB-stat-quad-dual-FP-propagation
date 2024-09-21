#include "solver.h"
#include "transform.h"

struct exception {
    std::string t;
};

void verify(bool b, std::string t = {}) {
    if (!b)
        throw exception{ std::move(t) };
}

__attribute__((const)) int approx_sqrt(std::size_t i) noexcept {
   return 1 << ((std::numeric_limits<std::size_t>::digits - __builtin_clzg(i)) / 2);
}

typedef Eigen::Matrix<double, mu_dim, Eigen::Dynamic> mu_matrix;
double write_mu(const PDE_info& info, mu_matrix& mu, const Eigen::Matrix<double, dim + dim, Eigen::Dynamic>& pq, const RS_nearest& RS) {
    const Eigen::Index N = pq.cols() - 1;
    mu_matrix old_mu = mu;
    for (Eigen::Index j = 0; j != N; ++j) {
        auto p0 = pq.col(j).topRows<dim>();
        auto q0 = pq.col(j).bottomRows<dim>();
        const auto& [R0, S0] = RS[2 * j];
        const auto& [R1, S1] = RS[2 * j + 1];
        Eigen::Matrix<double, mu_dim, 1> y;
        y.topRows<M_dim>().noalias() = info.M_0 * (p0 + R0 * q0);
        y.bottomRows<L_dim>().noalias() = info.L_0.transpose() * (S0 * q0);
        mu.col(2 * j) = info.C_hat_eta_inv(y);
        Eigen::Matrix<double, dim, 1> p1 = 0.5 * (p0 + pq.col(j + 1).topRows<dim>());
        Eigen::Matrix<double, dim, 1> q1 = 0.5 * (q0 + pq.col(j + 1).bottomRows<dim>());
        y.topRows<M_dim>().noalias() = info.M_0 * (p1 + R1 * q1);
        y.bottomRows<L_dim>().noalias() = info.L_0.transpose() * (S1 * q1);
        mu.col(2 * j + 1) = info.C_hat_eta_inv(y);
    }
    auto p0 = pq.col(N).topRows<dim>();
    auto q0 = pq.col(N).bottomRows<dim>();
    const auto& [R0, S0] = RS[2 * N];
    Eigen::Matrix<double, mu_dim, 1> y;
    y.topRows<M_dim>().noalias() = info.M_0 * (p0 + R0 * q0);
    y.bottomRows<L_dim>().noalias() = info.L_0.transpose() * (S0 * q0);
    mu.col(2 * N) = info.C_hat_eta_inv(y);
    mu = info.C_inv.solve(mu);
    return (mu - old_mu).norm();
}
Eigen::Matrix<double, dim, Eigen::Dynamic> solve_P(const PDE_info& info, const Eigen::Matrix<double, dim, 1>& p, const Eigen::Matrix<double, dim, 1>& q, std::span<const double> t) {
    typedef Eigen::Matrix<double, 2 * dim, 1> state_t;
    state_t tmp;
    tmp.topRows<dim>() = p;
    tmp.bottomRows<dim>() = q;
    Eigen::Matrix<double, dim, Eigen::Dynamic> P(dim, t.size());
    int i = 0;
    using namespace boost::numeric::odeint;
    auto stepper = make_controlled<runge_kutta_dopri5<state_t>>(1e-10, 1e-10);
    integrate_times(stepper, closed_loop{info}, tmp, t.begin(), t.end(), 1e-5, [&](const state_t& state, double) {
        P.col(i++) = state.topRows<dim>();
    });
    return P;
}

Eigen::Matrix<double, 2 * dim, Eigen::Dynamic> solve_PQ(const PDE_info& info, const Eigen::Matrix<double, dim, 1>& p, const Eigen::Matrix<double, dim, 1>& q, std::span<const double> t) {
    typedef Eigen::Matrix<double, 2 * dim, 1> state_t;
    state_t tmp;
    tmp.topRows<dim>() = p;
    tmp.bottomRows<dim>() = q;
    Eigen::Matrix<double, 2 * dim, Eigen::Dynamic> P(2 * dim, t.size());
    int i = 0;
    using namespace boost::numeric::odeint;
    auto stepper = make_controlled<runge_kutta_dopri5<state_t>>(1e-10, 1e-10);
    integrate_times(stepper, closed_loop{info}, tmp, t.begin(), t.end(), 1e-5, [&](const state_t& state, double) {
        P.col(i++) = state;
    });
    return P;
}

template <typename T>
inline auto split(const Eigen::DenseBase<T>& pq) {
    return std::make_pair(pq.template topRows<dim>(), pq.template middleRows<dim>(dim));
}

std::tuple<double, Eigen::Matrix<double, dim, 1>, Eigen::Matrix<double, dim, 1>>
solver::solve(const Eigen::Matrix<double, dim, 1>& x, std::span<const double> t, std::span<const std::size_t> N, std::span<const double> errs) const {
    std::vector<RS_nearest> RS;
    RS.reserve(t.size() - 1);
    for (std::size_t i = 0; i != t.size() - 1; ++i)
        RS.emplace_back(info, t[i], t[i + 1], 2 * N[i]);
    using namespace boost::numeric::odeint;
    typedef Eigen::Matrix<double, 2 * dim, 1> state_t;
    Eigen::Matrix<double, dim + dim, Eigen::Dynamic> PQ(dim + dim, t.size());
    auto P = PQ.topRows<dim>();
    auto Q = PQ.bottomRows<dim>();
    Q.setZero(); // FIXME
    P = solve_P(info, x, RS[0].S_begin() * Q.col(0), t);
    for (double err : errs) {
        bool down = false;
        int k = t.size() - 1;
        for (std::size_t counter = 0; counter < 500 * t.size() * t.size(); ++counter) {
            if (down) {
                down = false;
                ++k;
            }
            else if (k > 0)
                --k;
            else
                break;
            if (k == t.size() - 1) {
                Q.col(k + 1).setZero(); // FIXME
                continue;
            }
            pre_discretized_rk4_stepper stepper((t[k + 1] - t[k]) / N[k]);
            RS_closed_loop pq_ode{info, RS[k]};
            Eigen::Matrix<double, dim, 1> lambda = (k + 1 < t.size() - 1) ? (RS[k + 1].S_begin() * Q.col(k + 1)).eval() : Q.col(k + 1);
            for (uint_fast8_t i = 0; i != 255; ++i) {
                stepper.reset();
                state_t pq = PQ.col(k);
                for (std::size_t j = 0; j != N[k]; ++j)
                    stepper.do_step(pq_ode, pq, 0, 0);
                auto [p, q] = split(pq);
                Eigen::Matrix<double, dim, 1> z = p + RS[k].R_end() * q;
                if ((z - P.col(k + 1)).norm() > err) {
                    P.col(k + 1) = z;
                    down = true;
                    break;
                }
                else {
                    Eigen::Matrix<double, dim, 1> corr = lambda - pq.bottomRows<dim>();
                    Q.col(k) += corr;
                    if (corr.squaredNorm() < errs.back() * errs.back())
                        break;
                }
            }
        }
    }
    Eigen::Matrix<double, 2 * dim + 1, 1> state;
    state.topRows<dim>() = x;
    state.middleRows<dim>(dim) = RS[0].S_begin() * Q.col(0);
    state[2 * dim] = 0;
    {
    auto stepper = make_controlled<runge_kutta_dopri5<decltype(state)>>(1e-10, 1e-10);
    integrate_adaptive(stepper, closed_loop_with_r{ info },
        state, t.front(), t.back(), 1e-5);
    }
    return std::make_tuple(state[2 * dim], RS[0].S_begin() * Q.col(0), state.middleRows<dim>(dim).eval());
}
