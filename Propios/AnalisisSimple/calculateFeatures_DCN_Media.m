function [empaticaFeatureTable, caneFeatureTable] = calculateFeatures_DCN_viejoConMedia(empaticaDataExperiment, caneDataExperiment, participant_ID, winTime, winOverlapTime)
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

        % Extraemos PPG del experimento i
        canePPG = caneDataExperiment.(experimentNames{i}).PPG.data(1:ppg_length);
        empaticaPPG = empaticaDataExperiment.(experimentNames{i}).PPG.data(1:ppg_length);
        % Normalizamos
        canePPG = (canePPG-min(canePPG))/(max(canePPG)-min(canePPG));
        empaticaPPG = (empaticaPPG-min(empaticaPPG))/(max(empaticaPPG)-min(empaticaPPG));
    
        %% GENERAL FEATURES
        % Number of windows in total
        numberOfWins = floor( (ppg_length - winLengthPPG) / winLengthPPGOverlap + 1);
        % Experiment Number: Relax1 (0) - Music (1) - Relax2 (2) - Calculo (3)
        experimentID = repmat(experimentNumber{i}, numberOfWins, 1);
        % Participant ID
        subjectID = repmat(participant_ID, numberOfWins, 1);
        % Label: Estrés (1) - No-estrés (0)
        label = repmat(labelWeight{i}, numberOfWins, 1);
        
        %% PPG FEATURES

        % Magnitude Caracteristics
        % Feature 1: Max PPG value
        ppgMaxValue = NaN(numberOfWins,2);
        % Feature 2: Min PPG value
        ppgMinValue = NaN(numberOfWins,2); 
        % Feature 3: Mean PPG value
        ppgMean = NaN(numberOfWins,2);
        % Feature 4: Median PPG value
        ppgMedian = NaN(numberOfWins,2);
        % Feature 5: Standar Deviation PPG value
        ppgSD = NaN(numberOfWins,2);
        % Feature 6: InterQuartilic Range PPG value
        ppgIQR = NaN(numberOfWins,2);
        % Feature 7: Mean Absolute Deviation PPG value
        ppgMAD_Mean = NaN(numberOfWins,2);
        % Feature 8: Median Absolute Deviation PPG value
        ppgMAD_Median = NaN(numberOfWins,2);
        % Feature 9: absolute Mean Difference PPG value
        ppgMD = NaN(numberOfWins,2);
        
        % Time Caracteristics
        % Feature 10:  Mean IBI
        MeanIBI = NaN(numberOfWins,2);
        % Feature 11: Heart Rate
        HR = NaN(numberOfWins,2);
        % Feature 12: Standar Deviation of IBI
        SDNN = NaN(numberOfWins,2);
        % Feature 13: Root Mean Square of Sucesives IBI interval Diference
        RMSSD = NaN(numberOfWins,2);
        % Feature 14: Sucesive Differences Standardar Deviation
        SDSD = NaN(numberOfWins,2);
        % Feature 15: Percent of IBI larger than 20ms
        pNN20 = NaN(numberOfWins,2);
        % Feature 16: Percent of IBI larger than 50ms
        pNN50 = NaN(numberOfWins,2);
        % Feature 17: Baevsky's Stress Index
        SI = NaN(numberOfWins,2);

        % Frecuency Caracteristics
        % Feature 18: Low Frecuency Power
        LF = NaN(numberOfWins,2);
        % Feature 19: High Frecuency Power
        HF = NaN(numberOfWins,2);
        % Feature 20: Low Frecuency to High Frecuency Rate
        LFHF = NaN(numberOfWins,2);
    
        %% PPG FEATURE EXTRACTION
        countWin = 0;
        for j=1:winLengthPPGOverlap:ppg_length - winLengthPPG + 1
    
            % Window number
            countWin = countWin + 1;
            
            % Signal inside the window
            ppgCaneWin = canePPG(j:j+winLengthPPG-1);
            ppgEmpaticaWin = empaticaPPG(j:j+winLengthPPG-1);

            % Magnitude Caracteristics
            % Feature 1: Max PPG value
            ppgMaxValue(countWin,1) = max(ppgCaneWin);
            ppgMaxValue(countWin,2) = max(ppgEmpaticaWin);
            % Feature 2: Min  PPG value
            ppgMinValue(countWin,1) = min(ppgCaneWin);
            ppgMinValue(countWin,2) = min(ppgEmpaticaWin);
            % Feature 3: Mean PPG value
            ppgMean(countWin,1) = mean(ppgCaneWin);
            ppgMean(countWin,2) = mean(ppgEmpaticaWin);
            % Feature 4: Median PPG value
            ppgMedian(countWin,1) = median(ppgCaneWin);
            ppgMedian(countWin,2) = median(ppgEmpaticaWin);
            % Feature 5: Standar Deviation PPG value
            ppgSD(countWin,1) = std(ppgCaneWin);
            ppgSD(countWin,2) = std(ppgEmpaticaWin);
            % Feature 6: InterQuartilic Range PPG value
            ppgIQR(countWin,1) = iqr(ppgCaneWin);
            ppgIQR(countWin,2) = iqr(ppgEmpaticaWin);
            % Feature 7: Mean Absolute Deviation PPG value
            ppgMAD_Mean(countWin,1) = mad(ppgCaneWin);
            ppgMAD_Mean(countWin,2) = mad(ppgEmpaticaWin);
            % Feature 8: Median Absolute Deviation PPG value
            ppgMAD_Median(countWin,1) = mad(ppgCaneWin,1);
            ppgMAD_Median(countWin,2) = mad(ppgEmpaticaWin,1);
            % Feature 9: absolute Mean Difference PPG value
            ppgMD(countWin,1) = mean(mean(abs(ppgCaneWin(:) - ppgCaneWin(:)')));
            ppgMD(countWin,2) = mean(mean(abs(ppgEmpaticaWin(:) - ppgEmpaticaWin(:)')));

            % Time Caracteristics
            % Calculamos los cruces por cero (semiperiodo)
            zeroCane=[];
            zeroEmpatica=[];
            cane_umbral=mean(ppgCaneWin);
            empatica_umbral=mean(ppgCaneWin);
            for idx = 2:length(ppgCaneWin)-1
                % Paso por cero
                if (ppgCaneWin(idx) > cane_umbral && ppgCaneWin(idx+1) < cane_umbral) || (ppgCaneWin(idx) < cane_umbral && ppgCaneWin(idx+1) > cane_umbral)
                    zeroCane = cat(1,zeroCane,idx);
                end
            end
            for idx = 2:length(ppgEmpaticaWin)-1
                % Paso por cero
                if (ppgEmpaticaWin(idx) > empatica_umbral && ppgEmpaticaWin(idx+1) < empatica_umbral) || (ppgEmpaticaWin(idx) < empatica_umbral && ppgEmpaticaWin(idx+1) > empatica_umbral)
                    zeroEmpatica = cat(1,zeroEmpatica,idx);
                end
            end

            % IBI individual
            timestampCane       = caneDataExperiment.(experimentNames{i}).PPG.TimeStampPosix;
            timestampEmpatica   = empaticaDataExperiment.(experimentNames{i}).PPG.TimeStampPosix;
            timestampCaneWin        = timestampCane(j:j+winLengthPPG-1);
            timestampEmpaticaWin    = timestampEmpatica(j:j+winLengthPPG-1);
            IBI_Cane = 2*(timestampCaneWin(zeroCane(2:end))         -timestampCaneWin(zeroCane(1:end-1)));
            IBI_Empatica = 2*(timestampEmpaticaWin(zeroEmpatica(2:end)) -timestampEmpaticaWin(zeroEmpatica(1:end-1)));
            % Normalizamos%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %habrá que extraer el IBI de toda la gráfica a la vez
            IBI_Cane = (IBI_Cane-min(IBI_Cane))/(max(IBI_Cane)-min(IBI_Cane));
            IBI_Empatica = (IBI_Empatica-min(IBI_Empatica))/(max(IBI_Empatica)-min(IBI_Empatica));

            % Feature 10: Mean IBI
            MeanIBI(countWin,1) = mean(IBI_Cane);
            MeanIBI(countWin,2) = mean(IBI_Empatica);
            % Feature 11: Heart Rate
            HR(countWin,1) = 60/MeanIBI(countWin,1);
            HR(countWin,2) = 60/MeanIBI(countWin,2);
            % Feature 12: Standar Deviation of IBI
            SDNN(countWin,1) = sqrt(sum((IBI_Cane-MeanIBI(countWin,1)))^2/(winLengthPPG-1));
            SDNN(countWin,2) = sqrt(sum((IBI_Empatica-MeanIBI(countWin,2)))^2/(winLengthPPG-1));
            % Feature 13: Root Mean Square of Sucesives IBI interval diference
            RMSSD(countWin,1) = sqrt(sum((IBI_Cane(2:end)-IBI_Cane(1:end-1)).^2)/(winLengthPPG-1));
            RMSSD(countWin,2) = sqrt(sum((IBI_Empatica(2:end)-IBI_Empatica(1:end-1)).^2)/(winLengthPPG-1));
            % Feature 14: Sucesive Differences Standardar Deviation
            SucDif_Cane=abs(IBI_Cane(2:end)-IBI_Cane(1:end-1));
            SDSD(countWin,1) = sqrt(sum((SucDif_Cane-mean(SucDif_Cane)).^2)/(winLengthPPG-2));
            SucDif_Empatica=abs(IBI_Empatica(2:end)-IBI_Empatica(1:end-1));
            SDSD(countWin,2) = sqrt(sum((SucDif_Empatica-mean(SucDif_Empatica)).^2)/(winLengthPPG-2));
            % Feature 15: Percent of IBI larger than 20ms
            pNN20(countWin,1) = sum(IBI_Cane>20)/winLengthPPG*100;
            pNN20(countWin,2) = sum(IBI_Empatica>20)/winLengthPPG*100;
            % Feature 16: Percent of IBI larger than 50ms
            pNN50(countWin,1) = sum(IBI_Cane>50)/winLengthPPG*100;
            pNN50(countWin,2) = sum(IBI_Empatica>50)/winLengthPPG*100;
            % Feature 17: Baevsky's Stress Index
            % Calculamos el hsitograma primero: incremeto de 50ms
            t_intervalos=0.05;
            if(length(IBI_Cane)>2)
                caneintervalos=(min(IBI_Cane):t_intervalos:max(IBI_Cane))+t_intervalos/2;% Puntos centrales de los rangos
                canehist = histcounts(IBI_Cane, caneintervalos);
                [caneAMo_pre, canehist_idx] = max(canehist);
                caneMo = caneintervalos(canehist_idx);
                caneAMo = caneAMo_pre/winLengthPPG;
                caneMxDMn = max(IBI_Cane) - min(IBI_Cane);
                SI(countWin,1) = (caneAMo/(2*caneMo*caneMxDMn))*100; % Para el bastón 
            else
                SI(countWin,1) = 0;
            end
            if(length(IBI_Empatica)>2)
                empaticaintervalos=(min(IBI_Empatica):t_intervalos:max(IBI_Empatica))+t_intervalos/2;% Puntos centrales de los rangos
                empaticahist = histcounts(IBI_Empatica, empaticaintervalos);
                [empaticaAMo_pre, empaticahistidx] = max(empaticahist);
                empaticaMo = (empaticaintervalos(empaticahistidx) + empaticaintervalos(empaticahistidx+1)) / 2;
                empaticaAMo = empaticaAMo_pre/winLengthPPG;
                empaticaMxDMn = max(IBI_Empatica) - min(IBI_Empatica);
                SI(countWin,2) = (empaticaAMo/(2*empaticaMo*empaticaMxDMn))*100; % Para la empatica
            else
                SI(countWin,2) = 0;
            end

            % Frecuency Caracteristics
            VL_cutoff = 0.04 / (Fs_PPG / 2); % VL (0-0.04Hz)
            L_cutoff = 0.15 / (Fs_PPG / 2);  % L (0.04-0.15Hz)
            H_cutoff = 0.4 / (Fs_PPG / 2);  % H (0.15-0.4Hz)
            % Diseño del filtro Butterworth de orden 4
            [b_low, a_low] = butter(4, [VL_cutoff, L_cutoff], 'bandpass');
            [b_high, a_high] = butter(4, [L_cutoff, H_cutoff], 'bandpass');
            % Feature 18: Low Frecuency Power
            LF(countWin,1) = mean(abs(filter(b_low, a_low, ppgCaneWin)).^2) / Fs_PPG;
            LF(countWin,2) = mean(abs(filter(b_low, a_low, ppgEmpaticaWin)).^2) / Fs_PPG;
            % Feature 19: High Frecuency Power
            HF(countWin,1) = mean(abs(filter(b_high, a_high, ppgCaneWin)).^2) / Fs_PPG;
            HF(countWin,2) = mean(abs(filter(b_high, a_high, ppgEmpaticaWin)).^2) / Fs_PPG;
            % Feature 20: Low Frecuency to High Frecuency Power Rate
            LFHF(countWin,1) = LF(countWin,1)/HF(countWin,1);
            LFHF(countWin,2) = LF(countWin,2)/HF(countWin,2);
        end     
    
        %% FEATURE ORGANIZATION TABLE
        % 1. Creamos una tabla TEMPORAL de las características de PPG de la Empatica y del bastón
        caneFeatureTable_temp = table( ...
            ppgMaxValue(:,1), ...
            ppgMinValue(:,1), ...
            ppgMean(:,1), ...
            ppgMedian(:,1), ...
            ppgSD(:,1), ...
            ppgIQR(:,1), ...
            ppgMAD_Mean(:,1), ...
            ppgMAD_Median(:,1), ...
            ppgMD(:,1), ...
            ...
            MeanIBI(:,1), ...
            HR(:,1), ...
            SDNN(:,1), ...
            RMSSD(:,1), ...
            SDSD(:,1), ...
            pNN20(:,1), ...
            pNN50(:,1), ...
            SI(:,1), ...
            ...
            LF(:,1), ...
            HF(:,1), ...
            LFHF(:,1), ...
            ...
            subjectID, ...
            experimentID, ...
            label, ...
            'VariableNames',{'ppgMaxValue','ppgMinValue', ...
            ... % <------- Completar con todas las caracteristicas
            'ppgMean','ppgMedian', 'ppgSD', 'ppgIQR', 'ppgMAD_Mean', 'ppgMAD_Median', 'ppgMD', ...
            ...
            'MeanIBI','HR','SDNN','RMSSD','SDSD','pNN20','pNN50','SI','LF','HF','LFHF',...
            ...
            'subjectID', 'experimentID', 'label'});
    
        empaticaFeatureTable_temp = table( ...
            ppgMaxValue(:,2), ...
            ppgMinValue(:,2), ...
            ...
            ... % <------- Completar con todas las caracteristicas
            ppgMean(:,2), ...
            ppgMedian(:,2), ...
            ppgSD(:,2), ...
            ppgIQR(:,2), ...
            ppgMAD_Mean(:,2), ...
            ppgMAD_Median(:,2), ...
            ppgMD(:,2), ...
            ...
            MeanIBI(:,2), ...
            HR(:,2), ...
            SDNN(:,2), ...
            RMSSD(:,2), ...
            SDSD(:,2), ...
            pNN20(:,2), ...
            pNN50(:,2), ...
            SI(:,2), ...
            ...
            LF(:,2), ...
            HF(:,2), ...
            LFHF(:,2), ...
            ...
            subjectID, ...
            experimentID, ...
            label, ...
            'VariableNames',{'ppgMaxValue','ppgMinValue', ...
            ... % <------- Completar con todas las caracteristicas
            'ppgMean','ppgMedian', 'ppgSD', 'ppgIQR', 'ppgMAD_Mean', 'ppgMAD_Median', 'ppgMD', ...
            ...
            'MeanIBI','HR','SDNN','RMSSD','SDSD','pNN20','pNN50','SI','LF','HF','LFHF',...
            ...
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