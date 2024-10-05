function Psi_dot = eqn_A_bar(t, Psi, P, A, Gammatilde)
    dim = size(A, 1);
    Psi = reshape(Psi, dim, dim);
    Psi_dot = (squeeze(P(t)) * Gammatilde - A') * Psi;
    Psi_dot = Psi_dot(:);
end