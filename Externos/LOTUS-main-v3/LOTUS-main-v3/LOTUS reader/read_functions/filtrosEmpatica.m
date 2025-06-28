function [dat_filt] = filtrosEmpatica(dat)
%FILTROSEMPATICA Filtra la señal GSR y obtiene las señales tónica y fásica
%   Este filtro aplica un filtro exponencial y un filtro de media para
%   separar la señal en sus componentes tónica y fásica.

% Verificar si dat es una estructura
if ~isstruct(dat)
    error('Error: La entrada debe ser una estructura.');
end

% Verificar si 'EDA' existe
if isfield(dat, 'EDA') 

    dat_filt = dat;

    x = dat.EDA{1}.data;
    y1 = [0; NaN(length(x)-1, 1)];
    y2 = [0; NaN(length(x)-1, 1)];

    % Aplicamos filtro exponencial a GSR (fs = 4 Hz , fc = 1.60 Hz)
    % que representa el filtro analogico RC del circuito de la fásica.
    alpha = 0.715365;
    for i=2:length(x)
        if ~isnan(x(i)) && ~isnan(y1(i-1))
            y1(i) = alpha*x(i) + (1-alpha)*y1(i-1);
        elseif ~isnan(x(i)) && isnan(y1(i-1))
            y1(i) = x(i);
        end
    end
    % Aplicamos un segundo filtro exponencial a y1 que representa el
    % filtro digital implementado en PSoC a la señal.
    for i=2:length(y1)
        if ~isnan(y1(i)) && ~isnan(y2(i-1))
            y2(i) = alpha*y1(i) + (1-alpha)*y2(i-1);
        elseif ~isnan(x(i)) && isnan(y1(i-1))
            y2(i) = y1(i);
        end
    end

    % Obtenemos la señal tónica de la GSR con un filtro de media de frecuencia
    % de corte igual a fc = 0.05 Hz con una fs = 4 Hz.
    % |H(w)| = (1/N) * | sin(N·w/2) / sin(w/2) |
    % w = 2·pi·fc/fs
    % N = 35 para una frecuencia de corte (-3dB) en 0.05 Hz
    N = 35;
    b = (1/N)*ones(1,N);
    a = 1;
    TONIC = filter(b, a, y2);
    PHASIC = y2 - TONIC;

    % Modificar la tabla EDA con la columna data por la señal GSR filtrada
    tabla_EDA = dat_filt.EDA{1};
    tabla_EDA.data = y2;
    dat_filt.EDA{1} = tabla_EDA;

    % Añadir dat_filt.TONIC que contiene la señal tónica en la columna data
    tabla_TONIC = tabla_EDA;
    tabla_TONIC.data = TONIC;
    dat_filt.TONIC{1} = tabla_TONIC;

    % Añadir dat_filt.PHASIC que contiene la señal fásica en la columna data
    tabla_PHASIC = tabla_EDA;
    tabla_PHASIC.data = PHASIC;
    dat_filt.PHASIC{1} = tabla_PHASIC;

end

end

