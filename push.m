function [rho,u] = push(rho,u,grid)

% Push with KEP scheme (rho*u*u constraint)
% Forward Euler Timestep
% CD Space derivatives (Unstable?)
% Should be horendously diffusive but conserve KE

%SSP-RK3 (3 stage)
rho0 = rho;
rho0_u0 = rho.*u;
rho_u = rho.*u;
%Stage 1:
[rho, rho_u] = stage1(rho, rho_u, grid);
%Stage 2:
[rho, rho_u] = stage2(rho, rho_u, rho0, rho0_u0, grid);
%Stage 3:
[rho, rho_u] = stage3(rho, rho_u, rho0, rho0_u0, grid);

%Retrieve u:
u = rho_u./rho;
end


%SSP-RK3 (stage 3)
function [rho, rho_u] = stage3(rho, rho_u, rho0, rho0_u0, grid)

%First stage calculation
[rho_star, rho_u_star] = euler(rho, rho_u, grid);
rho = (1/3)*rho0 + (2/3)*rho_star;
rho_u = (1/3)*rho0_u0 + (2/3)*rho_u_star;

end

%SSP-RK3 (stage 2)
function [rho, rho_u] = stage2(rho, rho_u, rho0, rho0_u0, grid)

%First stage calculation
[rho_star, rho_u_star] = euler(rho, rho_u, grid);
rho = (3/4)*rho0 + (1/4)*rho_star;
rho_u = (3/4)*rho0_u0 + (1/4)*rho_u_star;

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
gamma = sqrt(1+u.*u);
v = u./gamma; 
rho_n = rho;
rho_u_n = rho_u;
v_j = v;
v_j_plus1 = v(R);
v_j_minus1 = v(L);

%Central averages to get ...
pv_j_plus_half = 0.5*(rho(R).*v(R) + rho.*v);
pv_j_minus_half = 0.5*(rho(L).*v(L) + rho.*v);

%Coef + or -
coef_plus = ( 1./gamma - 1./gamma(R) )./(v_j - v_j_plus1);
coef_minus = ( 1./gamma(L) - 1./gamma )./(v_j_minus1 - v_j);

%Update rho, u
% Euler Update:
rho_n_plus1 = rho_n - c*(pv_j_plus_half - pv_j_minus_half);
rho_u_n_plus1 = rho_u_n + (c*coef_plus.*pv_j_plus_half ...
    - c*coef_minus.*pv_j_minus_half);

%Extract rho and u
rho = rho_n_plus1;
rho_u = rho_u_n_plus1;

end