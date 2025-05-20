% plot3d_full_save.m
% ==================
% Carrega T_volume.dat, monta volume, faz 3 plots (fatias, isosuperfície e scatter do cubo)
% e salva cada figura em PNG.

%% 1) Load dos dados
data = load('dat/T_volume.dat');  % ajuste o caminho se necessário
coords = data(:,1:3);
Tvals  = data(:,4);

%% 2) Detecta dimensões da malha
nX = max(coords(:,1));
nY = max(coords(:,2));
nZ = max(coords(:,3));

%% 3) Parâmetros físicos e grid para fatias/isosuperfície
L  = 0.1;  
x = linspace(0,L,nX);
y = linspace(0,L,nY);
z = linspace(0,L,nZ);
[X,Y,Z] = meshgrid(x,y,z);

%% 4) Prepara T no formato [y x z] p/ meshgrid
T = zeros(nX,nY,nZ);
for idx = 1:size(coords,1)
    i = coords(idx,1);
    j = coords(idx,2);
    k = coords(idx,3);
    T(i,j,k) = Tvals(idx);
end
Tm = permute(T, [2,1,3]);

%% 5) Plot 1: slices em x=L/2, y=L/2, z=L/2
fig1 = figure('Name','Fatias de Temperatura','NumberTitle','off');
hs = slice(X,Y,Z,Tm, L/2, L/2, L/2);
set(hs,'EdgeColor','none','FaceColor','interp')
colormap(jet), colorbar
xlabel('x (m)'), ylabel('y (m)'), zlabel('z (m)')
title('Fatias em x=L/2, y=L/2, z=L/2')
daspect([1 1 1]), view(3)
saveas(fig1, 'img/fatias_temperatura.png')

%% 6) Plot 2: isosuperfície na média de T
isoVal = mean(Tvals);
fig2 = figure('Name','Isosuperfície Média','NumberTitle','off');
p = patch(isosurface(X,Y,Z,Tm, isoVal));
isonormals(X,Y,Z,Tm,p)
set(p,'FaceColor','red','EdgeColor','none','FaceAlpha',0.6)
colormap(jet), colorbar
xlabel('x (m)'), ylabel('y (m)'), zlabel('z (m)')
title(sprintf('Isosuperfície T = %.2f', isoVal))
camlight('headlight'), lighting gouraud
daspect([1 1 1]), view(3)
saveas(fig2, 'img/isosuperficie_temperatura.png')

%% 7) Plot 3: scatter3 de todo o cubo
% converte índices (i,j,k) em coordenadas físicas
dx = L/(nX-1);
dy = L/(nY-1);
dz = L/(nZ-1);
Xp = (coords(:,1)-1)*dx;
Yp = (coords(:,2)-1)*dy;
Zp = (coords(:,3)-1)*dz;

fig3 = figure('Name','Scatter 3D do Cubo','NumberTitle','off');
scatter3(Xp, Yp, Zp, 20, Tvals, 'filled')
colormap(jet), colorbar
xlabel('x (m)'), ylabel('y (m)'), zlabel('z (m)')
title('Scatter3D – Temperatura no Cubo Completo')
daspect([1 1 1]), view(3)
saveas(fig3, 'img/scatter_temperatura.png')