function [grid] = diagnostics(rho,u,grid)

%Compute E
E = (1/2) * rho.*u.*u;
fprintf("KE TOTAL: %1.12f\n",sum(E)*grid.dx/grid.E0);
grid.E_vs_t(grid.iter) = sum(E);

% Error:
if (0)
    %3-stage
    % %for dt = (1/1)*0.98*grid.dx/50;
    %KE_old = 0.999999996715
    % %1/2 dt: O(3)
    %KE_new = 0.999999999586

    % 2-stage
    % KE_old = 1.000001270097;
    % %1/2 dt: O(2)
    % KE_new = 1.000000311896;
    %Error = 1/((1 - KE_old)/( 1 - KE_new));
end

% Run only at select iterations:
if (mod ( grid.iter, grid.Output_interval ) == 0 || grid.iter == grid.NT)

    % Clear the figure
    %clf()

    color = [ (grid.NT - grid.iter)/grid.NT,0,(grid.iter)/grid.NT];


    %Plot the diagnostic output comparison to fig1 GA Sod
    subplot(2,2,1)
    plot(grid.x,rho,"color",color)
    hold on
    title("Density")
    ylabel("Density")
    xlabel("x")

    subplot(2,2,2)
    plot(grid.x,u,"color",color)
    hold on
    title("Velocity")
    ylabel("Velocity")
    xlabel("x")

    subplot(2,2,3)
    plot(grid.x,E,"color",color)
    title("Kinetic Energy (x)")
    xlabel("x")
    ylabel("KE")
    hold on

    subplot(2,2,4)
    plot(grid.time_vec(1:grid.iter),grid.E_vs_t(1:grid.iter),"*")
    title("Kinetic Energy (t)")
    ylabel("KE")
    xlabel("t")

    pause(0.01)

end
