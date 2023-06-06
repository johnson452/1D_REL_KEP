function [rho,u,grid] = make_grid()

%%% Initialize memory %%%
%[DEFAULT] Setup Grid: (Boundary Grid):
grid.Nx = 100; % Only specified here
grid.xmin = 0;
grid.xmax = 1.0;
grid.time = 0;
grid.t_max = 0.35;
grid.Output_interval = 100;
Nx = grid.Nx;

%[DEFAULT] Constants, updated in IC.m
grid.iter = 1;

%[DEFAULT] Grids, updated in IC.m
grid.dx = (grid.xmax - grid.xmin)/grid.Nx;
grid.time = 0;
grid.cfl = 0.98; %clf = udt/dx <= C_max
grid.dt = (1/4)*0.98*grid.dx/50;
grid.NT = ceil(grid.t_max/grid.dt);
grid.L = (grid.xmax - grid.xmin);

%Grid
grid.x = linspace(grid.xmin,grid.xmax,Nx);

%Quantities
rho = 1.5 + sin((2*grid.x*pi)*((Nx-1)/Nx)/(grid.L));
u = 1.5 + sin((2*grid.x*pi)*((Nx-1)/Nx)/(grid.L));


%Right and left
grid.R = mod( linspace(1,Nx,Nx), Nx) + 1; %Good
grid.L = mod( linspace(-1,Nx-2,Nx), Nx) + 1; %Good

% KE vs t
grid.E_vs_t = zeros(1,grid.NT);
grid.time_vec = linspace(0,grid.t_max,grid.NT);

gamma = sqrt(1+u.^2);
KE = (gamma - 1).*rho;
grid.E0 = sum(KE).*grid.dx;

end