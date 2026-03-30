% This code is for 2D Steady State Heat Conduction equation when we have 
% Neumann and Dirichlet Boundary Conditions using Node based approach
% UNIFORM GRID & NODE BASED APPROACH %

clc;
clear;
close all;

% Input Parameters
Nx = 81;         % Grid Points in X direction 
Ny = 81;         % Grid Points in Y direction
Lx = 1;         % Length in x
Ly = 1;         % Length in y
dx = Lx/(Nx-1); % Grid Spacing in X direction
dy = Ly/(Ny-1); % Grid Spacing in Y direction

tic
% Grid points calculations in x and y directions
x = linspace(0, Lx, Nx); 
y = linspace(0, Ly, Ny);

%Boundary conditions
T_Top   = 100;    % Top boundary in deg 
T_Right = 50;    % Right boundary in deg

Ti = (T_Right+T_Top)/2; %Initial Temp
T = Ti*ones(Nx,Ny); %Initialising grid points

% Corner points BCs
T(1,1) = T_Top;
% T(Nx,1) = ;
T(1, Nx) = (T_Right+T_Top)/2;
T(Ny,Nx) = T_Right;

T(:,end) = T_Right;       % right boundary (Dirichlet)
T(end,:) = T_Top;         % top boundary (Dirichlet)
T(:,1) = T(:,2);          % left insulated: mirror column 2 into column 1
T(1,:) = T(2,:);          % bottom insulated: mirror row 2 into row 1
T(1,1) = 0.5 * (T(1,2) + T(2,1));


tol=1e-6;
error=1;

iter=0;
max_iter=100000;
while (error > tol && iter<=max_iter)
    iter=iter+1;
    % disp(iter);
    T_old = T; % Store the old temperature values
 
    errormax=0;
    errorcomp=0;
    for i = 2:Nx-1
        for j = 2:Ny-1
            T(i, j) = 0.25*(T(i-1,j)+T(i+1,j)+T(i,j-1) +T(i,j+1)); % Heat eqn in algebraic form
            diff=abs(T(i, j) - T_old(i,j));
            if diff>errorcomp
                errormax=diff;
            end
             errorcomp = max(errorcomp, errormax); % Updating maximum error
       end
    end 

% Boundary Conditions 
    % Forcing BCs - Setting ghost cells for consistent BCs
    T(:,1)   = T(:,2);   % approximation dT/dx = 0    
    T(1,:)   = T(2,:);   % approximation dT/dy = 0
    T(:,end) = T_Right;  
    T(end,:) = T_Top;
    T(Ny,Nx) = 0.5*(T_Right+T_Top);
    T(1,1) = 0.5 * (T(1,2) + T(2,1));

    error=errormax;
   % disp(error)
    % error = 1e-5
    if error < tol
        sprintf("Tolerence Met at %d", iter)
        break
    end
end


%Centreline Temperatures 
mid_x = round(Nx/2);   % Rounding to get an integer when Nx is odd
mid_y = round(Ny/2);   % Rounding to get an integer when Ny is odd

T_x = T(mid_x, :); 
T_y = T(:, mid_y);  
T_mid = T(mid_x, mid_y)

toc

figure;
plot(y, T_x, 'r-', 'LineWidth', 2); hold on;
plot(x, T_y, 'b-', 'LineWidth', 2);
xlabel('Position');
ylabel('Temperature (°C)');
axis equal       
axis square 
title('Numerical Centreline Temperatures - Node Based Uniform');
legend('Along x (y = Ly/2)','Along y (x = Lx/2)','Location','north');
grid on;

T;
figure;
levels = linspace(min(T(:)),max(T(:)),30);
contourf(x,y,T,levels);
colormap(jet(length(levels)-1));
colorbar;
axis equal       
axis square 
xlabel('x', 'FontSize', 14);
ylabel('y', 'FontSize', 14);
title('2d Steady-State Heat Conduction - Node Based Uniform')
 
[X, Y] = meshgrid(x,y);

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

