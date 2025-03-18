clc; clear; close all;

dados = load('residuo.dat');

iteracoes = dados(:, 1);
erro = dados(:, 2);

figure;
semilogy(iteracoes, erro, 'b-', 'LineWidth', 2);
xlabel('Iterações');
ylabel('Resíduo');
title('Resíduo x Iterações');
grid on;

saveas(gcf, 'residuo.png');