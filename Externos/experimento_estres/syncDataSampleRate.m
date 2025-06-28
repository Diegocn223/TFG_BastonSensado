function [empaticaDataExperiment_sync,caneDataExperiment_sync] = syncDataSampleRate(empaticaDataExperiment,caneDataExperiment)
% syncDataSampleRate Sincroniza las señales de CaneSense y Empatica ajustando la tasa de muestreo
%   Esta función toma como entrada dos estructuras de datos: `empaticaDataExperiment` y `caneDataExperiment`, que contienen las señales correspondientes a diferentes experimentos (Relax1, Music, Relax2, Arithmetic).
%   La función realiza un ajuste en la tasa de muestreo de las señales: submuestrea la señal PPG de Empatica de 64 Hz a 32 Hz, y submuestrea las señales GSR, TONIC y PHASIC de CaneSense de 32 Hz a 4 Hz. Los datos ajustados se devuelven en las estructuras de salida `empaticaDataExperiment_out` y `caneDataExperiment_out`.
%
%   Ejemplo de uso:
%   [empaticaDataExperiment_out, caneDataExperiment_out] = syncDataSampleRate(empaticaDataExperiment, caneDataExperiment);
%   En este ejemplo, `empaticaDataExperiment` y `caneDataExperiment` son las estructuras de datos originales, y la función devuelve dos nuevas estructuras con las señales submuestreadas.
%
% Entradas:
%   - empaticaDataExperiment: Estructura de datos que contiene las señales de Empatica (PPG) para cada experimento (Relax1, Music, Relax2, Arithmetic).
%   - caneDataExperiment: Estructura de datos que contiene las señales de CaneSense (GSR, TONIC, PHASIC) para cada experimento (Relax1, Music, Relax2, Arithmetic).
%
% Salidas:
%   - empaticaDataExperiment_out: Estructura de datos de Empatica con la señal PPG submuestreada a 32 Hz.
%   - caneDataExperiment_out: Estructura de datos de CaneSense con las señales GSR, TONIC y PHASIC submuestreadas a 4 Hz.

if isempty(caneDataExperiment) || isempty(empaticaDataExperiment)
    error('Error: Ningúna variable de entrada puede estar vacía.');
end

if ~isstruct(empaticaDataExperiment) || ~isstruct(caneDataExperiment)
    error('Error: Formato de las variables de entrada incorrecto.');
end

if numel(fieldnames(empaticaDataExperiment)) ~= 4 || numel(fieldnames(caneDataExperiment)) ~= 4
    error('Error: Los datos struct deben contener exactamente 4 campos (Relax1, Music, Relax2, Arithmetic).');
end

experimentNames = {'Relax1', 'Music', 'Relax2', 'Arithmetic'};
caneDataExperiment_sync = caneDataExperiment;
empaticaDataExperiment_sync = empaticaDataExperiment;

% Submuestreo:
% - factor 2 de la señal PPG de la Empatica de 64 Hz a 32 Hz.
% - factor 8 de la señal GSR, TONIC, PHASIC del bastón de 32 Hz a 4 Hz.
for i = 1:numel(experimentNames)
    empaticaDataExperiment_sync.(experimentNames{i}).PPG = empaticaDataExperiment.(experimentNames{i}).PPG(1:2:end, :);
    caneDataExperiment_sync.(experimentNames{i}).GSR = caneDataExperiment.(experimentNames{i}).GSR(1:8:end, :);
    caneDataExperiment_sync.(experimentNames{i}).TONIC = caneDataExperiment.(experimentNames{i}).TONIC(1:8:end, :);
    caneDataExperiment_sync.(experimentNames{i}).PHASIC = caneDataExperiment.(experimentNames{i}).PHASIC(1:8:end, :);   
end

end