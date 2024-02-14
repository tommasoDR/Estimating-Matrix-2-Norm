% Script that plots the function in three dimension

A = rand(3,2) * 4;
Q = transpose(A) * A;

% Define the function
f = @(x1, x2) -(x1.^2 * Q(1,1) + x1.*x2 * Q(1,2) + x2.^2 * Q(2,2) ...
    + x1.*x2 * Q(2,1)) / (x1.^2 + x2.^2);

% Create a grid of x1 and x2 values
[x1, x2] = meshgrid(-10:0.1:10, -10:0.1:10);

% Evaluate the function for each pair of x1 and x2
z = f(x1, x2);

% Create a surface plot
figure;
fsurf(f, [-100 100 -100 100]);

% Add labels and title
xlabel('x1');
ylabel('x2');
zlabel('f(x)');
title('Surface Plot of f(x)');
