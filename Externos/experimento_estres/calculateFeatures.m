function [empaticaFeatureTable, caneFeatureTable] = calculateFeatures(empaticaDataExperiment, caneDataExperiment, participant_ID, winTime, winOverlapTime)
% calculateFeatures Extrae y organiza características fisiológicas de señales Empatica y CaneSense
%   Esta función toma como entrada dos estructuras de datos: `empaticaDataExperiment` y `caneDataExperiment`, que contienen las señales correspondientes a diferentes experimentos (Relax1, Music, Relax2, Arithmetic). 
%   A partir de estas señales, la función segmenta los datos en ventanas temporales, y calcula un conjunto de características fisiológicas (features) de interés, tanto para la señal PPG como para las señales EDA (GSR, TONIC, PHASIC). 
%
%   Ejemplo de uso:
%   [empaticaFeatureTable, caneFeatureTable] = calculateFeatures(empaticaDataExperiment, caneDataExperiment, 1, 15, 3.75);
%   En este ejemplo, los datos se dividen en ventanas de 15 segundos con solapamiento de 3.75 segundos, y se devuelven las tablas de características para Empatica y CaneSense.
%
% Entradas:
%   - empaticaDataExperiment: Estructura de datos que contiene las señales de Empatica (PPG, GSR, TONIC, PHASIC) para cada experimento (Relax1, Music, Relax2, Arithmetic).
%   - caneDataExperiment: Estructura de datos que contiene las señales de CaneSense (PPG, GSR, TONIC, PHASIC) para cada experimento.
%   - participant_ID: Número que identifica al sujeto del experimento.
%   - winTime: Duración de cada ventana en segundos.
%   - winOverlapTime: Tiempo de solapamiento entre ventanas en segundos.
%
% Salidas:
%   - empaticaFeatureTable: Tabla con las características extraídas de las señales de Empatica, estructurada por ventanas, con identificadores y etiquetas.
%   - caneFeatureTable: Tabla con las características extraídas de las señales de CaneSense, estructurada por ventanas, con identificadores y etiquetas.


%% Validaciones
    if ~isnumeric(participant_ID) || ~isscalar(participant_ID) || participant_ID < 0 || participant_ID > 99 || mod(participant_ID,1) ~= 0
        error('El valor de participant_ID debe ser un número entero entre 0 y 99 (inclusive).');
    end
    if ~isnumeric(winTime) || ~isscalar(winTime) || winTime <= 0
        error('El valor de winTime debe ser un número real positivo mayor que 0.');
    end

    if ~isnumeric(winOverlapTime) || ~isscalar(winOverlapTime) || winOverlapTime <= 0
        error('El valor de winOverlapTime debe ser un número real positivo mayor que 0.');
    end

    if winTime <= winOverlapTime
        error('El valor de winTime debe ser mayor que winOverlapTime.');
    end

    if isempty(empaticaDataExperiment) || ~isstruct(empaticaDataExperiment)
        error('empaticaDataExperiment debe ser un struct no vacío.');
    end
    if numel(fieldnames(empaticaDataExperiment)) ~= 4
        error('empaticaDataExperiment debe contener exactamente 4 campos (Relax1, Music, Relax2, Arithmetic).');
    end

    if isempty(caneDataExperiment) || ~isstruct(caneDataExperiment)
        error('caneDataExperiment debe ser un struct no vacío.');
    end
    if numel(fieldnames(caneDataExperiment)) ~= 4
        error('caneDataExperiment debe contener exactamente 4 campos (Relax1, Music, Relax2, Arithmetic).');
    end

    experimentNames = {'Relax1', 'Music', 'Relax2', 'Arithmetic'};
    experimentNumber = {0, 1, 2, 3};
    labelWeight = {0, 1, 0, 1};
    
    Fs_PPG = 32; % Hz
    Fs_GSR = 4; % Hz
    
    winLengthPPG = winTime * Fs_PPG;
    winLengthPPGOverlap = round(winOverlapTime * Fs_PPG);
    
    winLengthGSR = winTime * Fs_GSR;
    winLengthGSROverlap = round(winOverlapTime * Fs_GSR);



    for i=1:numel(experimentNames)
        % Para asegurarnos que las dos señales (empatica y cane) tienen el mismo número de elementos:
        %% PPG signal
        ppg_length = min([length(empaticaDataExperiment.(experimentNames{i}).PPG.data), length(caneDataExperiment.(experimentNames{i}).PPG.data)]);
        
        if ppg_length < winLengthPPG
            error('La longitud de la señal PPG de Empática/Bastón debe ser mayor que la longitud de la ventana winTime'); 
        end
        
        %% EDA signal
        gsr_length = min([length(empaticaDataExperiment.(experimentNames{i}).GSR.data), length(caneDataExperiment.(experimentNames{i}).GSR.data)]);
        
        if gsr_length < winLengthGSR
            error('La longitud de la señal GSR de Empática/Bastón debe ser mayor que la longitud de la ventana winTime'); 
        end

        % Verificamos el hipotético para que la señal de gsr y ppg no sea
        % lo suficientemente diferente como para que el número de ventanas
        % varíe
        
        if (ppg_length / Fs_PPG) < (gsr_length / Fs_GSR)
            gsr_length = gsr_length - 1;
        elseif (gsr_length / Fs_GSR) < (ppg_length / Fs_PPG)
            ppg_length = ppg_length - 1;
        end

        canePPG = caneDataExperiment.(experimentNames{i}).PPG.data(1:ppg_length);
        empaticaPPG = empaticaDataExperiment.(experimentNames{i}).PPG.data(1:ppg_length);

        caneGSR = caneDataExperiment.(experimentNames{i}).GSR.data(1:gsr_length);
        caneTONIC = caneDataExperiment.(experimentNames{i}).TONIC.data(1:gsr_length);
        canePHASIC = caneDataExperiment.(experimentNames{i}).PHASIC.data(1:gsr_length);
        empaticaGSR = empaticaDataExperiment.(experimentNames{i}).GSR.data(1:gsr_length);
        empaticaTONIC = empaticaDataExperiment.(experimentNames{i}).TONIC.data(1:gsr_length);
        empaticaPHASIC = empaticaDataExperiment.(experimentNames{i}).PHASIC.data(1:gsr_length);
    
        %% GENERAL FEATURES
        % Número de ventanas (PPG y EDA tienen el mismo número de ventanas)
        numberOfWins = floor( (ppg_length - winLengthPPG) / winLengthPPGOverlap + 1);
        % Experiment Number: Relax1 (0) - Music (1) - Relax2 (2) - Calculo (3)
        experimentID = repmat(experimentNumber{i}, numberOfWins, 1);
        % ID del participante
        subjectID = repmat(participant_ID, numberOfWins, 1);
        % Label: Estrés (1) - No-estrés (0)
        label = repmat(labelWeight{i}, numberOfWins, 1);
        
        %% PPG FEATURES
        % Feature 1: Max PPG value
        ppgMaxValue = NaN(numberOfWins,2);
        % Feature 2: Min PPG value
        ppgMinValue = NaN(numberOfWins,2);
        
        %% Añadir características... %%
        % Feature 3: R to R interval
        RR = NaN(numberOfWins,2);
        % Feature 4: Heart Rate
        HRV = NaN(numberOfWins,2);
        % Feature 5: Standar Deviation of RR
        SDNN = NaN(numberOfWins,2);
        % Feature 6: Root Mean Square of Sucesives RR interval diference
        RMSSD = NaN(numberOfWins,2);
        % Feature 7: Percent of RR larger than 20ms
        pNN20 = NaN(numberOfWins,2);
        % Feature 8: Percent of RR larger than 50ms
        pNN50 = NaN(numberOfWins,2);
        % Feature 9: Sucesive Differences Standardar Deviation
        SDSD = NaN(numberOfWins,2);
        % Feature 10: Baevsky's Stress Index
        SI = NaN(numberOfWins,2);
        % Feature 11: Low Frecuency Power
        LF = NaN(numberOfWins,2);
        % Feature 12: High Frecuency Power
        HF = NaN(numberOfWins,2);
        % Feature 13: Low Frecuency to High Frecuency Rate
        SI = NaN(numberOfWins,2);
        
        %% EDA FEATURES
        % Feature 3: Min GSR value
        gsrMinValue = NaN(numberOfWins,2);
        % Feature 4: Max TONIC value
        tonicMaxValue = NaN(numberOfWins,2);
        % Feature 5: STD of PHASIC
        phasicSTD = NaN(numberOfWins,2);
        %% Añadir características... %%
    
        %% PPG FEATURE EXTRACTION
        countWin = 0;
        for j=1:winLengthPPGOverlap:ppg_length - winLengthPPG + 1
    
            % Calculamos el número de la ventana
            countWin = countWin + 1;
            
            % Obtenemos la señal en la ventana correspondiente
            ppgCaneWin = canePPG(j:j+winLengthPPG-1);
            ppgEmpaticaWin = empaticaPPG(j:j+winLengthPPG-1);
    
            % Feature 1: Valor máximo de la ventana
            ppgMaxValue(countWin,1) = max(ppgCaneWin);
            ppgMaxValue(countWin,2) = max(ppgEmpaticaWin);
    
            % Feature 2: Valor mínimo de la ventana
            ppgMinValue(countWin,1) = min(ppgCaneWin);
            ppgMinValue(countWin,2) = min(ppgEmpaticaWin);
            
            %% Completar ... %%      
        
        end
        
        %% EDA FEATURE EXTRACTION
        countWin = 0;
        for j=1:winLengthGSROverlap:gsr_length - winLengthGSR + 1
    
            % Calculamos el número de la ventana
            countWin = countWin + 1;
            
            % Obtenemos la señal en la ventana correspondiente
            % GSR
            gsrCaneWin = caneGSR(j:j+winLengthGSR-1);
            gsrEmpaticaWin = empaticaGSR(j:j+winLengthGSR-1);
            % TONIC
            tonicCaneWin = caneTONIC(j:j+winLengthGSR-1);
            tonicEmpaticaWin = empaticaTONIC(j:j+winLengthGSR-1);
            % PHASIC
            phasicCaneWin = canePHASIC(j:j+winLengthGSR-1);
            phasicEmpaticaWin = empaticaPHASIC(j:j+winLengthGSR-1);
    
            % Feature 3: Valor mínimo de la ventana GSR
            gsrMinValue(countWin,1) = min(gsrCaneWin);
            gsrMinValue(countWin,2) = min(gsrEmpaticaWin);
    
            % Feature 4: Valor máximo de la ventana TONIC
            tonicMaxValue(countWin,1) = max(tonicCaneWin);
            tonicMaxValue(countWin,2) = max(tonicEmpaticaWin);
            
            % Feature 5: Valor STD de la ventana PHASIC
            phasicSTD(countWin,1) = std(phasicCaneWin);
            phasicSTD(countWin,2) = std(phasicEmpaticaWin);
            
            %% Completar ... %%      
     
        end
    
        %% FEATURE ORGANIZATION TABLE
        % 1. Creamos una tabla TEMPORAL de las características de PPG de la Empatica y del bastón
        caneFeatureTable_temp = table( ...
            ppgMaxValue(:,1), ...
            ppgMinValue(:,1), ...
            gsrMinValue(:,1), ...
            tonicMaxValue(:,1), ...
            phasicSTD(:,1), ...
            ... % <------- Completar con todas las caracteristicas
            subjectID, ...
            experimentID, ...
            label, ...
            'VariableNames',{'ppgMaxValue','ppgMinValue','gsrMinValue','tonicMaxValue','phasicSTD', ...
            ... % <------- Completar con todas las caracteristicas
            'subjectID', 'experimentID', 'label'});
    
        empaticaFeatureTable_temp = table( ...
            ppgMaxValue(:,2), ...
            ppgMinValue(:,2), ...
            gsrMinValue(:,2), ...
            tonicMaxValue(:,2), ...
            phasicSTD(:,2), ...
            ... % <------- Completar con todas las caracteristicas
            subjectID, ...
            experimentID, ...
            label, ...
            'VariableNames',{'ppgMaxValue','ppgMinValue','gsrMinValue','tonicMaxValue','phasicSTD', ...
            ... % <------- Completar con todas las caracteristicas
            'subjectID', 'experimentID', 'label'});
    
        % 2. Si es la primera iteración, inicializar las tablas finales
        if i == 1
            caneFeatureTable = caneFeatureTable_temp;
            empaticaFeatureTable = empaticaFeatureTable_temp;
        else
        % 3. Si las tablas finales han sido inicializadas, se va almacenando los datos
            caneFeatureTable = [caneFeatureTable; caneFeatureTable_temp];
            empaticaFeatureTable = [empaticaFeatureTable; empaticaFeatureTable_temp]; 
        end
    end
    % Por último, añado una columna que simplemente es el número de ventana
    caneFeatureTable.WinSample = (1:height(caneFeatureTable))';
    caneFeatureTable = movevars(caneFeatureTable, 'WinSample', 'Before', 1);

    empaticaFeatureTable.WinSample = (1:height(empaticaFeatureTable))';
    empaticaFeatureTable = movevars(empaticaFeatureTable, 'WinSample', 'Before', 1);

end