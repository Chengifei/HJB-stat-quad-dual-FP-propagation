function P_dot = eqn_6_10(~, P, A, Gammatilde, Sc1)
    dim = size(Sc1, 1);
    P = reshape(P, dim, dim);
    P_dot = -Sc1 - A' * P - P * A - P * Gammatilde * P;
    P_dot = P_dot(:);
end