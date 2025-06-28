


% Cargar datos
load('2.mat'); % Sacados de los bubles en "AnalisisProfundo"

% Parámetros de los ejes
x_vals = 1:6;              % Cantidad de estimulos
y_vals = 5:5:60;           % Tiempo de ventana (segundos)
labels_y = string(y_vals);

% Crear colormap blanco → carmesi oscuro
rojo_carmesi_oscuro = [0.55, 0, 0.11];
cmap = [linspace(1, rojo_carmesi_oscuro(1), 256)', ...
        linspace(1, rojo_carmesi_oscuro(2), 256)', ...
        linspace(1, rojo_carmesi_oscuro(3), 256)'];

% Rango común (minimo y máximo global)
vmin = min([recuento_ind(:); recuento_par(:)], [], 'omitnan');
vmax = max([recuento_ind(:); recuento_par(:)], [], 'omitnan');

% Crear figura
f=figure;
f.Position = [100 100 1400 500];

% --- Subplot 1: recuento_ind ---
subplot(1,2,1)
h1 = heatmap(x_vals, labels_y, recuento_ind, ...
    'MissingDataColor', [0 0 0], ...           % negro para NaN
    'MissingDataLabel', '', ...
    'Colormap', cmap, ...
    'ColorLimits', [vmin, vmax]);
h1.Title = 'Muestras independientes';
h1.XLabel = 'Cantidad de estimulos';
h1.YLabel = 'Tiempo de ventana (s)';
h1.FontSize = 14;  % Tamaño de fuente aumentado

% --- Subplot 2: recuento_par ---
subplot(1,2,2)
h2 = heatmap(x_vals, labels_y, recuento_par, ...
    'MissingDataColor', [0 0 0], ...
    'MissingDataLabel', '', ...
    'Colormap', cmap, ...
    'ColorLimits', [vmin, vmax]);
h2.Title = 'Muestras pareadas';
h2.XLabel = 'Cantidad de estimulos';
h2.YLabel = 'Tiempo de ventana (s)';
h2.FontSize = 14;  % Tamaño de fuente aumentado
