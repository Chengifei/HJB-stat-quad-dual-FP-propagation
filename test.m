Amat=[0 1.5 0 1 0; -1 0 0 -0.5 0.5; 1.5 0 0 -1 0; 0 -1 0 -1 0; 1 0 0 1 0];
dim=size(Amat,1);
M0=[0 1 0 0 0]; % Sets y_1 dim to be domain/range of  nonlinearity.
Gammatilde=eye(dim,dim);
Cmat=0.5*eye(dim,dim);
c1=-3; % coefficient related to stat-quad duality
ycenter=zeros(dim,1);
widthx1=3;
widthx2=4;
t0=0;       % initial time (generally zero). This is "t" in the paper.
Tbar=0.5;
tfullv=linspace(t0,Tbar,51);
[X1, X2] = meshgrid(linspace(-widthx1,widthx1,75), linspace(-widthx2,widthx2,50));
checkingstat = 0;
icomp = [1 2];
funcparams = {
    @(x) 0.5 * sin(2 * x);
    @(x) cos(2 * x);
    @(x) -2 * sin(2 * x);
    };
%%
[Warray, relbacksubserrpropmat, rel_error] = main(Amat, Cmat, c1, M0, [-3 1 -2], Gammatilde, tfullv, X1, X2, ycenter, icomp, [3 4 5], funcparams, true);
%%
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