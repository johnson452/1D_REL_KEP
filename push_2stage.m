function [rho,u] = push(rho,u,grid)

% Push with KEP scheme (rho*u*u constraint)
% Forward Euler Timestep
% CD Space derivatives (Unstable?)
% Should be horendously diffusive but conserve KE

%SSP-RK3 (4 stage)
rho0 = rho;
rho0_u0 = rho.*u;
rho_u = rho.*u;
%Stage 1:
[rho, rho_u] = stage1(rho, rho_u, grid);
%Stage 2:
[rho, rho_u] = stage2(rho, rho_u, rho0, rho0_u0, grid);

%Retrieve u:
u = rho_u./rho;
end


%SSP-RK3 (stage 2)
function [rho, rho_u] = stage2(rho, rho_u, rho0, rho0_u0, grid)

%First stage calculation
[rho_star, rho_u_star] = euler(rho, rho_u, grid);
rho = 0.5*rho0 + 0.5*rho_star;
rho_u = 0.5*rho0_u0 + 0.5*rho_u_star;

end

%SSP-RK3 (stage 1)
function [rho, rho_u] = stage1(rho, rho_u, grid)

%First stage calculation
[rho_star, rho_u_star] = euler(rho, rho_u, grid);
rho = rho_star;
rho_u = rho_u_star;

end

function [rho, rho_u] = euler(rho, rho_u, grid)

%Simplify the equations;
c = grid.dt/grid.dx;
R = grid.R;
L = grid.L;
u = rho_u./rho;
rho_n = rho;
rho_u_n = rho_u;
u_j = u;
u_j_plus1 = u(R);
u_j_minus1 = u(L);

%Central averages to get ...
pu_j_plus_half = 0.5*(rho(R).*u(R) + rho.*u);
pu_j_minus_half = 0.5*(rho(L).*u(L) + rho.*u);

%Update rho, u
% Euler Update:
rho_n_plus1 = rho_n - c*(pu_j_plus_half - pu_j_minus_half);
rho_u_n_plus1 = rho_u_n - (c*0.5*(u_j + u_j_plus1).*pu_j_plus_half ...
    - c*0.5*(u_j_minus1 + u_j).*pu_j_minus_half);

%Extract rho and u
rho = rho_n_plus1;
rho_u = rho_u_n_plus1;

end