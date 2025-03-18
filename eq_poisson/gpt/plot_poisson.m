function plot_poisson()
    % PLOTA OS RESULTADOS DO CÓDIGO FORTRAN DA EQUAÇÃO DE POISSON.
    % Certifique-se de que os arquivos:
    %   - solucao2d.dat
    %   - residuo.dat
    %   - corte_vertical.dat
    % estão no mesmo diretório em que este script é executado.

    % --------------------------
    % 1) Mapa de calor da solução 2D
    % --------------------------
    % Cada linha de solucao2d.dat tem (x, y, u)
    data2d = load('solucao2d.dat');
    xVals = unique(data2d(:,1));
    yVals = unique(data2d(:,2));
    Nx = length(xVals);
    Ny = length(yVals);

    % Transformar em grade/matriz (X, Y, U) para visualização
    % Obs: a ordem de reshape depende de como salvamos no Fortran
    % Mas se gravamos "i" em loop interno e "j" no loop externo,
    % então ao fazermos reshape(..., Nx, Ny), cada linha "j" estará
    % "contínua". Se inverter, troque Nx <-> Ny.
    [X, Y] = meshgrid(xVals, yVals);
    U = reshape(data2d(:,3), Nx, Ny);  % pode precisar Nx,Ny ou Ny,Nx

    figure;
    imagesc(xVals, yVals, U');  % Usamos U' pois (i,j) ~ (x,y)
    set(gca, 'YDir', 'normal'); % Ajustar para eixo y "para cima"
    colorbar;
    title('Mapa de Calor da Solução Numérica');
    xlabel('x');
    ylabel('y');

    % --------------------------
    % 2) Gráfico do resíduo x iteração
    % --------------------------
    resData = load('residuo.dat');  % (iter, res)
    figure;
    plot(resData(:,1), resData(:,2), '-o');
    grid on;
    xlabel('Iterações');
    ylabel('Resíduo');
    title('Convergência do Resíduo (Gauss-Seidel)');

    % --------------------------
    % 3) Corte vertical em x ~ 0.5
    % --------------------------
    cutData = load('corte_vertical.dat');  % (y, u_num, u_ex)
    figure;
    plot(cutData(:,1), cutData(:,2), 'o-'); hold on;
    plot(cutData(:,1), cutData(:,3), 'x-');
    grid on;
    xlabel('y');
    ylabel('u');
    legend('Solução Numérica','Solução Exata','Location','Best');
    title('Corte Vertical em x = 0.5');
end
