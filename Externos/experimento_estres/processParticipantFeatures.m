function [empaticaFeatures, caneFeatures] = processParticipantFeatures(empaticaData, caneData, protocolTime_datetime, winTime, winOverlapTime, participant_ID)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % Obtenemos las señales correspondientes a lso diferentes experimentos.
    [caneDataExperiment, empaticaDataExperiment] = getExperimentalSignal(empaticaData, caneData, protocolTime_datetime);

    % Re-muestreo de las señales para que tengan la misma Fs
    [empaticaDataExperiment_sync,caneDataExperiment_sync] = syncDataSampleRate(empaticaDataExperiment,caneDataExperiment);

    % Extracción de características
    [empaticaFeatures, caneFeatures] = calculateFeatures(empaticaDataExperiment_sync, caneDataExperiment_sync, participant_ID, winTime, winOverlapTime);
    
end