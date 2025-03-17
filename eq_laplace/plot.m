clc; clear; close all;

% Carregar os dados do arquivo .dat
dados = load('matriz.dat'); 

% Separar colunas
X = dados(:, 1);  % Índice X
Y = dados(:, 2);  % Índice Y
U = dados(:, 3);  % Valor da solução U(x,y)

% Número de pontos na malha (NX e NY)
NX = max(X);  % Último valor de X indica NX
NY = max(Y);  % Último valor de Y indica NY

% Criar matriz de solução U(X, Y)
U_matrix = reshape(U, [NX, NY]);  % Reformata os dados para matriz

% Criar grids de coordenadas X e Y
[X_grid, Y_grid] = meshgrid(1:NX, 1:NY);

% Plotar a solução como um mapa de calor (heatmap)
figure;
contourf(X_grid, Y_grid, U_matrix', 50, 'LineColor', 'none');
colorbar;
xlabel('X');
ylabel('Y');
title('Solução da Equação de Laplace');
colormap jet; % Mapa de cores mais visível

saveas(gcf, 'heatmap_2D.png');