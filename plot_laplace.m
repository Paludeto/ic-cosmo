% Lê os dados do arquivo .dat
data = load('matriz.dat'); 

% Separa as colunas x, y e temperatura T
x = data(:,1);  
y = data(:,2);  
T = data(:,3);  

% Determina o tamanho da malha
Nx = numel(unique(x));  
Ny = numel(unique(y));  

% Reorganiza os dados corretamente para um heatmap
T_grid = reshape(T, Ny, Nx);  

% Plota um heatmap 2D
figure;
imagesc(unique(x), unique(y), T_grid);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('X'); ylabel('Y');
title('Mapa de Temperatura');

% Salva o gráfico no mesmo diretório
saveas(gcf, 'heatmap_2D.png');