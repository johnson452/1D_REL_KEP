function [grid] = diagnostics(rho,u,grid)

%Compute E
gamma = sqrt(1+u.^2);
KE = (gamma - 1).*rho;
fprintf("KE TOTAL: %1.12f\n",sum(KE)*grid.dx/grid.E0);
grid.E_vs_t(grid.iter) = sum(KE);

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
    subplot(2,3,1)
    plot(grid.x,rho,"color",color)
    hold on
    title("Density")
    ylabel("Density")
    xlabel("x")

    subplot(2,3,2)
    plot(grid.x,u,"color",color)
    hold on
    title("Momentum")
    ylabel("Momentum")
    xlabel("x")

    subplot(2,3,3)
    plot(grid.x,KE,"color",color)
    title("Kinetic Energy (x)")
    xlabel("x")
    ylabel("KE")
    hold on

    subplot(2,3,4)
    plot(grid.time_vec(1:grid.iter),grid.E_vs_t(1:grid.iter),"*")
    title("Kinetic Energy (t)")
    ylabel("KE")
    xlabel("t")

    subplot(2,3,5)
    v = u./gamma;
    plot(grid.x,v,"color",color)
    title("Velocity")
    ylabel("v")
    xlabel("x")

    subplot(2,3,6)
    plot(grid.x,gamma,"color",color)
    title("gamma")
    ylabel("gamma")
    xlabel("x")

    pause(0.01)

end
