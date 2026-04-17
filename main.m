function [Warray, relbacksubserrpropmat, rel_error] = main(A, ...
    C, c1, M0, ynoncompv, Gammatilde, tfullv, X1, X2, ycenter, icomp, inoncompv, funcparams, bc)
assert(all(size(X1) == size(X2)));
t_rev = tfullv(end:-1:1); % reversed t vector
nt1=30; % number of time steps used in fixed-point iterations

odeopts = odeset('RelTol', 1e-7, 'AbsTol', 1e-12);

dim = size(A, 1); % dimension of state space

%% Solve the Riccati equation for P and cache the results
P_Psi_end = [zeros(dim); eye(dim)];
Sc1 = C-c1*M0'*M0;
[~, P_Psi] = ode45(@(t, P) eqn_P_Psi(t, P, A, Gammatilde, Sc1), t_rev, P_Psi_end(:), odeopts);
P_Psi = P_Psi(end:-1:1, :);
P_Psi = reshape(P_Psi, [], 2 * dim, dim);
P_Psi = permute(P_Psi, [2 3 1]);
Pv = P_Psi(1:dim, :, :);
Psiv = P_Psi(dim+1:end, :, :);

%% Compute the associated state transition matrix and cache the results
Abarv = pagemtimes(Pv, Gammatilde) - A';
assert(size(Psiv,3) == length(tfullv));
Phibarst=zeros(dim,dim,length(tfullv),length(tfullv));
for it0=1:length(tfullv)
    for it1=1:length(tfullv)
        Phibarst(:,:,it1,it0)=Psiv(:,:,it1) / Psiv(:,:,it0);
    end
end

%% Compute \breve{\mathcal{D}}_{\cdot, \cdot} (M^0)'
DMz = squeeze(pagemtimes(D_sigma_s(tfullv, Phibarst), M0'));

%% Main loop for computation
Warray=zeros([length(tfullv), size(X1')]);
relbacksubserrpropmat=zeros(size(X1'));

for ix1=1:size(X1, 2)
    for ix2=1:size(X2, 1)
        y = ycenter(:);
        y(icomp) = y(icomp) + [X1(ix2, ix1); X2(ix2, ix1)];
        y(inoncompv) = ynoncompv;
        [Warray(:,ix1,ix2), relbacksubserrpropmat(ix1, ix2)]=W(y,M0,Pv,Abarv,Phibarst,DMz,c1,tfullv,nt1,Gammatilde,funcparams);
    end % loop over ix2
end % loop over ix1

%% Back substitution error checks
if bc
    % Time derivative (finite-difference) using adjacent times
    Wt = diff(Warray, 1, 1) / (tfullv(2) - tfullv(1));
    dx = [X1(1, 2) - X1(1, 1); X2(2, 1) - X2(1, 1)];
    backsubchkW = zeros(size(Warray));
    scale = zeros(size(Warray));
    it = 1;
    % Gradient (finite difference) using adjacent grid points in plane
    % and additional evaluations for out-of-plane derivatives
    for ix1=2:size(X1, 2)-1
        for ix2=2:size(X2, 1)-1
            y = ycenter(:);
            y(icomp) = y(icomp) + [X1(ix2, ix1); X2(ix2, ix1)];
            y(inoncompv) = ynoncompv;
            xqterm=0.5*y'*C*y;
            gradv = zeros(dim, 1);
            gradv(icomp) = ([Warray(1, ix1 + 1, ix2); Warray(1, ix1, ix2 + 1)] - [Warray(1, ix1 - 1, ix2); Warray(1, ix1, ix2 - 1)]) .* (0.5 ./ dx);
            for jj = inoncompv
                y1 = y;
                y2 = y;
                y1(jj) = y1(jj) + 0.01;
                y2(jj) = y2(jj) - 0.01;
                Wv1 = W(y1,M0,Pv,Abarv,Phibarst,DMz,c1,tfullv,nt1,Gammatilde,funcparams);
                Wv2 = W(y2,M0,Pv,Abarv,Phibarst,DMz,c1,tfullv,nt1,Gammatilde,funcparams);
                gradv(jj) = (Wv1(1) - Wv2(1)) * (0.5 * 100);
            end
            bqterm=gradv'*A*y;
            pqterm=0.5*gradv'*Gammatilde*gradv;
            Ntildeatx=Ntildeviafxdpt(M0 * y, c1, funcparams);
            terms = [Wt(it, ix1, ix2), xqterm, bqterm, Ntildeatx, -pqterm];
            backsubchkW(it,ix1,ix2) = sum(terms);
            scale(it,ix1,ix2) = norm(terms, 1);
        end % loop over ix2
    end % loop over ix1
    rel_error = abs(backsubchkW) ./ scale;
end
end

function [Wv, relbsalpproperr] = W(y,M0,Pv,Abarv,Phibarst,DMz,c1,tfullv,nt1,Gammatilde,funcparams)
% Evaluates the value function at a point specified by y
% using the precomputed data Pv, Abarv, Phibarst, DMz, Gammatilde
%% Compute alpha process
    alpha0v = M0 * squeeze(pagemtimes(squeeze(Phibarst(:, :, 1, :)), "transpose", y, "none"));
    [allalpham,t1vfinal,relbsalpproperr]=...
        propalpnlpsiTbndim(M0,Abarv,Phibarst,...
        DMz,alpha0v,c1,tfullv,nt1,funcparams);
    alphafinalv=allalpham(end, :);
    if relbsalpproperr>0.1
        error('Max rel backsub err from propalp...Tbar >0.45. Pausing in lxl...')
    end
%% Evaluate the value function from alpha
    [~,~,qv,rv]=qrintegrals(t1vfinal,alphafinalv,Phibarst,...
        M0,c1,Gammatilde,funcparams);
    Wv = 0.5*squeeze(pagemtimes(y', pagemtimes(Pv, y)))' + y' * qv + 0.5 * rv;
end
