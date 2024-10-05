function P_Psi_dot = eqn_P_Psi(~, P_Psi, A, Gammatilde, Sc1)
    dim = size(Sc1, 1);
    P_Psi = reshape(P_Psi, 2 * dim, dim);
    P = P_Psi(1:dim, :);
    Psi = P_Psi(dim+1:end, :);
    P_dot = -Sc1 - A' * P - P * A + P * Gammatilde * P;
    Psi_dot = (P * Gammatilde - A') * Psi;
    P_Psi_dot = [P_dot; Psi_dot];
    P_Psi_dot = P_Psi_dot(:);
end