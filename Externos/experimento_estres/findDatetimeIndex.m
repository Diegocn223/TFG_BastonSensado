function [index_start, index_end] = findDatetimeIndex(datetime_init, datetime_end, datetime_signal)
% findDatetimeIndex Encuentra los índices de la ventana de tiempo en el vector datetime_signal
%   Esta función toma dos valores datetime (`datetime_init` y `datetime_end`) y un vector de tiempo (`datetime_signal`) y devuelve los índices correspondientes a la ventana de tiempo definida por `datetime_init` y `datetime_end`.
%   El índice de inicio (`index_start`) es el índice del primer valor en `datetime_signal` que es mayor que `datetime_init`, y el índice de fin (`index_end`) es el índice del último valor en `datetime_signal` que es menor que `datetime_end`.
%
% Entradas:
%   - datetime_init: Tiempo de inicio de la ventana en formato datetime.
%   - datetime_end: Tiempo final de la ventana en formato datetime.
%   - datetime_signal: Vector de tiempo en formato datetime de la señal a recortar.
%
% Salidas:
%   - index_start: Índice de inicio de la ventana calculado para la señal.
%   - index_end: Índice de fin de la ventana calculado para la señal.

% Verificar que datetime_init, datetime_end y datetime_signal sean de tipo datetime
    if ~isdatetime(datetime_init) || ~isdatetime(datetime_end) || ~isdatetime(datetime_signal)
        error('Todos los inputs deben ser del tipo datetime.');
    end
    
    % Verificar que datetime_signal sea un vector de datetime
    if numel(datetime_signal) < 2
        error('datetime_signal debe ser un vector con al menos 2 elementos.');
    end
    
    datetime_init.TimeZone = datetime_signal.TimeZone;  % Copia la zona horaria
    datetime_end.TimeZone = datetime_signal.TimeZone;

    datetime_init.Format = datetime_signal.Format;  % Copia el formato de visualización
    datetime_end.Format = datetime_signal.Format;

    % Encontrar el índice del datetime más cercano que sea mayor que datetime_init
    index_start = find(datetime_signal > datetime_init, 1, 'first');
    
    % Si no se encuentra un datetime mayor que datetime_init, lanzar un error con mensaje personalizado
    if isempty(index_start)
        error('No se han encontrado datos en la ventana de tiempo especificada.');
    end
    
    % Encontrar el índice del datetime más cercano que sea menor que datetime_end
    index_end = find(datetime_signal < datetime_end, 1, 'last');
    
    % Si no se encuentra un datetime menor que datetime_end, lanzar un error con mensaje personalizado
    if isempty(index_end)
        error('No se han encontrado datos en la ventana de tiempo especificada.');
    end
end
