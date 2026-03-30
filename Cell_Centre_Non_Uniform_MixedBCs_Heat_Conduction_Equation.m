% This code is for 2D Steady State Heat Conduction equation when we have 
% Neumann and Dirichlet Boundary Conditions using Cell Centred approach
% NON UNIFORM GRID & CELL CENTRE BASED APPROACH %

clc;
clear; 
close all;

% Input Parameters
Nx_nodes = 81;         % Grid Points in X direction 
Ny_nodes = 81;         % Grid Points in Y direction
Lx = 1;         % Length in x
Ly = 1;         % Length in y

% Grid points calculations in x and y directions
n=0.75;

% Boundary temperatures
T_Right  = 50;         % Dirichlet on right face
T_Top    = 100;        % Dirichlet on top face
Ti = (T_Right+T_Top)*0.5;

% Number of cells in each direction - Calculation
    for i = 1:Nx_nodes
        s = (i-1)/(Nx_nodes-1);
        x_nodes(i) = (Lx/2) * (1 + sign(s-0.5) * abs(2*s-1)^n);
    end

    % For Y-direction
    for i = 1:Ny_nodes
        t = (i-1)/(Ny_nodes-1);
        y_nodes(i) = (Ly/2) * (1 + sign(t-0.5) * abs(2*t-1)^n);
    end

Nx_cc = length(x_nodes) - 1;   % this will be 80 if Nx_nodes = 81
Ny_cc = length(y_nodes) - 1;

% cell-center coordinates (midpoints)
x_cc = 0.5 * ( x_nodes(1:end-1) + x_nodes(2:end) );   % length Nx
y_cc = 0.5 * ( y_nodes(1:end-1) + y_nodes(2:end) );   % length Ny

% Initialization of Cell cetred Grid points
Tcc = Ti * ones(Ny_cc+2, Nx_cc+2); 

% Initial enforcement of BCs (set ghost cells consistent with BCs)
Tcc(:,1) = Tcc(:,2);   % left ghost column
Tcc(1,:) = Tcc(2,:);   % bottom ghost row
Tcc(:, Nx_cc+2) = 2*T_Right - Tcc(:, Nx_cc+1);   % right ghost
Tcc(Ny_cc+2, :) = 2*T_Top   - Tcc(Ny_cc+1, :);   % top ghost

tol = 1e-6;
error = 1;
iter = 0;
max_iter = 100000;

while (error > tol && iter <= max_iter)
    iter = iter + 1;
    T_old = Tcc;    % store previous
    errormax = 0;
    errorcomp = 0;
    
    % Loop over interior cell centers: j -> y index, i -> x index
    for j = 2:Ny_cc+1
        for i = 2:Nx_cc+1
            
            kx = i - 1;  % cell index in 1..Nx corresponding to ghost-array i
            ky = j - 1;  % cell index in 1..Ny corresponding to ghost-array j
            
            % compute distances between adjacent cell centers in x and y
            if kx == 1
                dx_w = x_cc(1) - (x_cc(1) - (x_cc(2)-x_cc(1))); 
            else
                dx_w = x_cc(kx) - x_cc(kx-1);
            end
            if kx == Nx_cc
                dx_e = x_cc(Nx_cc) - (x_cc(Nx_cc) - (x_cc(Nx_cc)-x_cc(Nx_cc-1))); 
            else
                dx_e = x_cc(kx+1) - x_cc(kx);
            end
            
            if ky == 1
                dy_s = y_cc(1) - (y_cc(1) - (y_cc(2)-y_cc(1))); 
            else
                dy_s = y_cc(ky) - y_cc(ky-1);
            end
            if ky == Ny_cc
                dy_n = y_cc(Ny_cc) - (y_cc(Ny_cc) - (y_cc(Ny_cc)-y_cc(Ny_cc-1))); 
            else
                dy_n = y_cc(ky+1) - y_cc(ky);
            end
            
            % aE acts on T(i+1), aW acts on T(i-1)
            aE = 2 / ( dx_e * (dx_w + dx_e) );
            aW = 2 / ( dx_w * (dx_w + dx_e) );
            % aN acts on T(j+1), aS acts on T(j-1)
            aN = 2 / ( dy_n * (dy_s + dy_n) );
            aS = 2 / ( dy_s * (dy_s + dy_n) );
            % center coefficient (sum of neighbors)
            aP = aE + aW + aN + aS;
            
            Tcc(j,i) = ( aE * Tcc(j, i+1) + aW * Tcc(j, i-1) + aN * Tcc(j+1, i) + aS * Tcc(j-1, i) ) / aP;
            diff=abs(Tcc(i, j) - T_old(i,j));
            if diff>errorcomp
                errormax=diff;
            end
             errorcomp = max(errorcomp, errormax); % Updating maximum error
        end
    end
    
    % Forcing BCs - Setting ghost cells for consistent BCs
    Tcc(:,1) = Tcc(:,2);    % left ghost column
    Tcc(1,:) = Tcc(2,:);    % bottom ghost row
    Tcc(:, Nx_cc+2) = 2*T_Right - Tcc(:, Nx_cc+1);
    Tcc(Ny_cc+2, :) = 2*T_Top   - Tcc(Ny_cc+1, :);
    
  error = errormax;
    if error < tol
        fprintf("Tolerance met at %d iterations\n", iter);
        break
    end
end

% Extracting interior temperature matrix for plotting and analysis
    T_extracted = Tcc(2:Ny_cc+1, 2:Nx_cc+1); 

    % Adding Boundary values for plotting purpose
    xc = [0, x_cc, Lx];
    yc = [0, y_cc, Ly];
    Tc = Ti * ones(Nx_cc+2, Ny_cc+2);
    Tc(2:end-1,2:end-1) = T_extracted;

    Tc(:,1) = Tc(:,2);      % Left boundary (Neumann → zero gradient)
    Tc(1,:) = Tc(2,:);      % Bottom boundary (Neumann → zero gradient)
    Tc(:, Nx_cc+2) = T_Right;  % Right boundary (Dirichlet)
    Tc(Ny_cc+2, :) = T_Top;    % Top boundary (Dirichlet)
    Tc(Ny_cc+2, Nx_cc+2) = (T_Right+T_Top)/2; %Corner Point Crrection

    figure;
    levels = linspace(min(Tc(:)),max(Tc(:)),30);
    contourf(xc,yc,Tc,levels);
    colormap(jet(length(levels)-1));
    colorbar;
    axis equal       
    axis square 
    xlabel('x', 'FontSize', 14);
    ylabel('y', 'FontSize', 14);
    title('2d Steady-State Heat Conduction - Cell Centred Non Uniform')
    set(gca, 'FontSize', 16);  

    %Centreline Temperatures 
    mid_x = round(size(Tc,1)/2);   
    mid_y = round(size(Tc,2)/2); 
    T_x = Tc(mid_x, :); 
    T_y = Tc(:, mid_y);  
    T_mid = Tc(mid_x, mid_y);

    figure;
    plot(yc, T_x, 'r-', 'LineWidth', 2); hold on;
    plot(xc, T_y, 'b-', 'LineWidth', 2);
    xlabel('Position');
    ylabel('Temperature (°C)');
    axis equal       
    axis square 
    title('Numerical Centreline Temperatures - Cell Centred Non Uniform');
    legend('Along x (y = Ly/2)','Along y (x = Lx/2)','Location','north');
    grid on;


    [X, Y] = meshgrid(x_cc,y_cc);

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
