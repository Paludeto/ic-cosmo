data = load('dat/T_volume.dat');
scatter3(data(:,1), data(:,2), data(:,3), 20, data(:,4), 'filled');
colorbar;
title("Distribuição de temperatura ao longo da superfície")
saveas(gcf, "img/scatter3.png")