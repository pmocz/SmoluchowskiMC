clear all
close all
clc

%% Solve a Smoluchowski problem with monte carlo approach
% Philip Mocz (2019)
% Princeton University

% for some references, see:
% [1] finite volume approach of http://arxiv.org/pdf/1312.7240v1.pdf
% [2] Lecot (Quasi-Monte Carlo simulation of coagulation-fragmentation)  https://pdfs.semanticscholar.org/a302/b31dbaf28f31cbebcbf6f7b2692d391cad9e.pdf
% [3] http://waset.org/publications/pdfs/P005912011330.pdf

% method

% draw particles form initial distribution
% loop over all particles
% calculate reaction probabilities with particles above it (timestep such that reaction probability < alpha < 1)
% apply reactions / update masses

% The equation is:
% df/dt = 1/2 * int_0^x[K_A(x-y,y)f(x-y,t)f(y,t)]dy - int_0^inf[K_A(x,y)f(x,t)f(y,t)]dy
% where
% K_A(x,y) is the aggregation kernel
% Note that the integral of g := x*f is a conserved quantity



%% Consider a problem with known analytic solution

% [Example 1.] constant kernel: K = 1


% problem parameters
a = 0;
b = 100;
f0 = @(x) 4e3*exp(-x);
kernel = @(x,y) 0*x*y + 1;
source = 'none';
t_out = [0 0.5 1 3];
seed = 42;


% solve the problem
tic;
sol = smolMCsolve(a, b, f0, kernel, source, t_out, seed);
toc;

N = numel(sol{1});

% plot solution -- first the analytic solution
fh = figure;
x = linspace(a, b, 1000);
Nt = numel(t_out);
cc = 1;
for t = t_out
    f = N*(2/(2+t)).^2 .* exp(-2/(2+t).*x);
    semilogy(x, f, 'color', [(Nt-cc)/Nt 0 cc/Nt], 'linewidth',2);
    hold on
    cc = cc + 1;
end


% compare with numerical solution -- make histogram of sample
xbin = linspace(a,b,100);
dxbin = xbin(2) - xbin(1);
xcenter = 0.5 * (xbin(1:end-1) + xbin(2:end));
cc = 1;
for t = t_out
    h = histcounts(sol{cc},xbin,'normalization','countdensity');
    semilogy(xcenter, h, 'o', 'color', [(Nt-cc)/Nt 0 cc/Nt], 'linewidth',2);
    cc = cc + 1;
end


% figure final touches
xlabel('$x$','interpreter','latex','fontsize',14);
ylabel('$f(t,x)$','interpreter','latex','fontsize',14);
axis([0 12 N*1e-6 N*1e1])




%% [Example 2.] alternate case


% problem parameters
a = 5;
b = 100;
f0 = @(x) 5e4*x.^-2.35.*(x>5).*(x<50);
kernel = @(x,y) 0*x*y + 1;
source = 'none';
t_out = [0 10 20 40];
seed = 42;


% solve the problem
tic;
sol = smolMCsolve(a, b, f0, kernel, source, t_out, seed);
toc;

N = numel(sol{1});

% plot solution
fh = figure;
xbin = linspace(a,b+100,40);
dxbin = xbin(2) - xbin(1);
xcenter = 0.5 * (xbin(1:end-1) + xbin(2:end));
Nt = numel(t_out);
cc = 1;
for t = t_out
    h = histcounts(sol{cc},xbin,'normalization','countdensity');
    loglog(xcenter, h, 'o-', 'color', [(Nt-cc)/Nt 0 cc/Nt], 'linewidth',2); hold on
    cc = cc + 1;
end


% figure final touches
xlabel('$x$','interpreter','latex','fontsize',14);
ylabel('$f(t,x)$','interpreter','latex','fontsize',14);
axis([5 200 N*1e-4 N*1e0])
