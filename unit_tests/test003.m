% load x_coordinates, y_coordinates, u_velocity and v_velocity csv files from parent directory and contourf plot u_velocity field and v_velocity field

clear; clc; close all;

x = readmatrix('x_coordinates_0.csv');
y = readmatrix('y_coordinates_0.csv');

for i = 0:9
    u = readmatrix(strcat('u_velocity_', num2str(i), '.csv'));
    v = readmatrix(strcat('v_velocity_', num2str(i), '.csv'));

    % set the values of u and v to NAN within a radius of 1 around the origin
    u(abs(x + 1i*y) < 1) = NaN;
    v(abs(x + 1i*y) < 1) = NaN;

    % figure;
    % contourf(x, y, v, 200, 'LineColor', 'none');
    % shading flat;
    % colorbar;
    % colormap jet;
    % title('Transverse Field Contour Plot');
    % xlabel('x');
    % ylabel('y');
    % clim([-1,1]);
    % axis equal;

    figure(i+100);
    contourf(x, y, u, 200, 'LineColor', 'none');
    shading flat;
    colorbar;
    colormap jet;
    title(sprintf('Streamwise velocity contour plot Timestep %d', i));
    xlabel('x');
    ylabel('y');
    clim([-10,10]);
    axis equal;
    % 
    % % quiver plot the velocity field
    % figure;
    % quiver(x, y, u, v);
    % title('Quiver Plot');
    % xlabel('x');
    % ylabel('y');   
    % axis equal;
end

for ii = 1:10000
    jj = mod(ii, 10);
    figure(jj+100);
    pause(0.5);
end