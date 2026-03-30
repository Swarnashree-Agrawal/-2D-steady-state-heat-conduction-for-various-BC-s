% This code is for 2D Steady State Heat Conduction equation when we have 
% Neumann and Dirichlet Boundary Conditions using Cell Centred Approach
% UNIFORM GRID & CELL CENTRED APPROACH %

clc
clear
close all

% Input Parameters
Nx = 81;    % Number of cells in X direction
Ny = 81;    % Number of cells in Y direction
Lx = 1;    % Length in x
Ly = 1;    % Length in y
dx = Lx / Nx;
dy = Ly / Ny;

xc = (dx/2) : dx : (Lx - dx/2);
yc = (dy/2) : dy : (Ly - dy/2);

% Boundary conditions
T_Top   = 100;    % Top boundary (Dirichlet)
T_Right = 50;     % Right boundary (Dirichlet)
Ti = (T_Right + T_Top) / 2; % Initial Temp
T = Ti * ones(Nx+2, Ny+2);      % Cell-centered temperature array

% Initial enforcement of BCs (set ghost cells consistent with BCs)
T(:,1)   = T(:,2);   % approximation dT/dx = 0   Left Ghost
T(1,:)   = T(2,:);   % approximation dT/dy = 0   Bottom Ghost
T(:,Nx+2) = 2*T_Right - T(:,Nx+1);  % Right Ghost
T(Ny+2,:) = 2*T_Top   - T(Ny+1,:);  % Top Ghost

tol = 1e-6;
error = 1;
iter = 0;
max_iter = 100000;

while (error > tol && iter <= max_iter)
    iter = iter + 1;
    T_old = T;
    errormax = 0;
    errorcomp = 0;
    for i = 2:Nx+1
        for j = 2:Ny+1
            T(i, j) = 0.25*(T(i-1, j) + T(i+1, j) + T(i, j-1) + T(i, j+1));
             diff=abs(T(i, j) - T_old(i,j));
            if diff>errorcomp
                errormax=diff;
            end
             errorcomp = max(errorcomp, errormax); % Updating maximum error
        end
    end


% Forcing BCs - Setting ghost cells for consistent BCs
T(:,1)   = T(:,2);   % approximation dT/dx = 0    
T(1,:)   = T(2,:);   % approximation dT/dy = 0
T(:,Nx+2) = 2*T_Right - T(:,Nx+1);  % right ghost
T(Ny+2,:) = 2*T_Top   - T(Ny+1,:);  % top ghost

    error = errormax;
    if error < tol
        fprintf("Tolerance met at %d iterations\n", iter);
        break
    end
end

    % Extracting interior temperature matrix for plotting and analysis
    T_extracted = T(2:Ny+1, 2:Nx+1); 

    % Adding Boundary values for plotting purpose
    xcc = [0, xc, Lx];
    ycc = [0, yc, Ly];
    Tc = Ti * ones(Nx+2, Ny+2);
    Tc(2:end-1,2:end-1) = T_extracted;

    Tc(:,1) = Tc(:,2);      % Left boundary (Neumann → zero gradient)
    Tc(1,:) = Tc(2,:);      % Bottom boundary (Neumann → zero gradient)
    Tc(:, Nx+2) = T_Right;  % Right boundary (Dirichlet)
    Tc(Ny+2, :) = T_Top;    % Top boundary (Dirichlet)
    Tc(Ny+2, Nx+2) = (T_Right+T_Top)/2; %Corner Point Crrection

    figure;
    levels = linspace(min(Tc(:)),max(Tc(:)),30);
    contourf(xcc,ycc,Tc,levels);
    colormap(jet(length(levels)-1));
    colorbar;
    axis equal       
    axis square 
    xlabel('x', 'FontSize', 14);
    ylabel('y', 'FontSize', 14);
    title('2d Steady-State Heat Conduction - Cell Centred Uniform')
    set(gca, 'FontSize', 16);  

    %Centreline Temperatures 
    mid_x = round(size(Tc,1)/2);   
    mid_y = round(size(Tc,2)/2); 
    T_x = Tc(mid_x, :); 
    T_y = Tc(:, mid_y);  
    T_mid = Tc(mid_x, mid_y);

    figure;
    plot(ycc, T_x, 'r-', 'LineWidth', 2); hold on;
    plot(xcc, T_y, 'b-', 'LineWidth', 2);
    xlabel('Position');
    ylabel('Temperature (°C)');
    axis equal       
    axis square 
    title('Numerical Centreline Temperatures - Cell Centred Uniform');
    legend('Along x (y = Ly/2)','Along y (x = Lx/2)','Location','north');
    grid on;

    [X, Y] = meshgrid(xc,yc);
    figure;
    axis equal       
    axis square       
    plot(X(:), Y(:), 'ko','MarkerSize',3)
    xlabel('Position');
    ylabel('Spacing');
    axis equal       
    axis square  
    title('Distribution');
    grid on;
