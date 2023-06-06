%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grant Johnson
% 6/2/2023
% 1D Euler Kep Based on Ammar's KEP Notes

%See readme for math
%built-in-periodic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

% Main function
[rho,u,grid] = make_grid();

%Make the diagnostic Figure
figure('units','normalized','outerposition',[0 0 0.5 0.5])

%%% Time loop %%%
while(grid.time < grid.t_max)
    
    %Call i/o and diagnostics
    grid = diagnostics(rho,u,grid);
    
    %Update the gridtime
    grid.time = grid.time + grid.dt;
    
    %Update the iterator
    grid.iter = grid.iter + 1;
 
    %Updater - updates all quantities simultaneosly
    % n -> n + 1 all quantities
    [rho,u] = push(rho,u,grid);

end
%%% End Time Loop %%%
%%% End main %%%