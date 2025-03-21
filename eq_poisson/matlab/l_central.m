clc; clear; close all;

dados = load('dat/l_central.dat');
x = dados(:, 1);
y = dados(:, 2);

% Função analítica
f = @(xCoord, yCoord) yCoord .* (xCoord.^2 - xCoord .* yCoord);

figure;

% Solução numérica (x, y)
plot(x, y, 'b--', 'LineWidth', 2);
hold on;

% Solução analítica em (nX / 2) + 1
plot(x, f(0.5, x), 'r-', 'LineWidth', 1);

xlabel("Coluna");
ylabel("Linha");
title("Solução numérica vs. analítica na linha central");
grid on;

% Legenda para identificar cada curva
legend("Numérica", "Analítica", 'Location', 'best');

% Salva a figura em arquivo
saveas(gcf, "img/l_central.png");