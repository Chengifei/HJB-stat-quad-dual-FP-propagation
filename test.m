% The driver script

% PROBLEM DEFINITION
t0=0;       % initial time (generally zero). This is "t" in the paper.
Tbar=0.5;   % Terminal time
tfullv=linspace(t0,Tbar,51);
% A matrix in dynamics
Amat=[0 1.5 0 1 0; -1 0 0 -0.5 0.5; 1.5 0 0 -1 0; 0 -1 0 -1 0; 1 0 0 1 0];
% extract state space dimension
dim=size(Amat,1);
% M^0, sets the domain of nonlinearity.
M0=[0 1 0 0 0];
% \tilde\Gamma in paper
Gammatilde=eye(dim,dim);
% C matrix in cost
Cmat=0.5*eye(dim,dim);
% coefficient related to stat-quad duality
c1=-3;

% REGION OF INTEREST (currently only handling planes parallel to coordinate axes)
% affine offset of the plane where the value will be evaluated
ycenter=zeros(dim,1); % not fully implemented --- see call signature of main
% width of the two parameters into the affine plane
widthx1=3;
widthx2=4;

[X1, X2] = meshgrid(linspace(-widthx1,widthx1,3), linspace(-widthx2,widthx2,2));
% checkingstat = 0;
icomp = [1 2]; % select indices of two components where x1 and x2 varies
funcparams = { % Nonlinearity
    @(x) 0.5 * sin(2 * x); % \Theta
    @(x) cos(2 * x); % \Theta'
    @(x) -2 * sin(2 * x); % \Theta''
    };
%%
[Warray, relbacksubserrpropmat, rel_error] = ...
    main(Amat, Cmat, c1, M0, [-3 1 -2], Gammatilde, ...
        tfullv, X1, X2, ycenter, icomp, [3 4 5], funcparams, true);

%% Plotting
surf(X1(1, :), X2(:, 1), squeeze(Warray(1,:,:))', 'EdgeAlpha', 0.35);
xlabel('$$y_1$$', 'interpreter', 'latex')
ylabel('$$y_2$$', 'interpreter', 'latex')
set(gcf, 'Position', [100 100 290 225])
view(75, 15)
daspect([1 1 1.2])
% title('Value function')
exportgraphics(gca, 'Value.pdf', 'ContentType', 'vector')
surf(X1, X2, squeeze(rel_error(1, :, :))', 'EdgeAlpha', 0.35);
% title('Relative back-substitution error (PDE)')
xlabel('$$y_1$$', 'interpreter', 'latex')
ylabel('$$y_2$$', 'interpreter', 'latex')
set(gcf, 'Position', [100 100 300 225])
view(75, 15)
daspect([1 1 0.002])
exportgraphics(gca, 'rel_error.pdf', 'ContentType', 'vector')
surf(X1(1, :), X2(:, 1),squeeze(relbacksubserrpropmat(:,:))' * 1e3, 'EdgeAlpha', 0.35);
% title('Relative back-substitution error into (9.1)')
xlabel('$$y_1$$', 'interpreter', 'latex')
ylabel('$$y_2$$', 'interpreter', 'latex')
zlabel('$$\times 10^{-3}$$', 'interpreter', 'latex')
set(gcf, 'Position', [100 100 300 225])
view(75, 15)
daspect([1 1 1.1])
exportgraphics(gca, 'Relativebacksubsplot_c1_0p5.pdf', 'ContentType', 'vector')
