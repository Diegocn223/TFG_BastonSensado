function [empaticaFeatureTable, caneFeatureTable] = calculateFeatures_DCN_2(empaticaDataExperiment, caneDataExperiment, participant_ID, winTime, winOverlapTime)
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

        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Normalizamos
        canePPG = canePPG*(max(canePPG)-min(canePPG))-min(canePPG);
        empaticaPPG = empaticaPPG*(max(empaticaPPG)-min(empaticaPPG))-min(empaticaPPG);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}

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
        MeanRR = NaN(numberOfWins,2);
        % Feature 4: Heart Rate
        HR = NaN(numberOfWins,2);
        % Feature 5: Standar Deviation of RR
        SDNN = NaN(numberOfWins,2);
        % Feature 6: Root Mean Square of Sucesives RR interval diference
        RMSSD = NaN(numberOfWins,2);
        % Feature 7: Sucesive Differences Standardar Deviation
        SDSD = NaN(numberOfWins,2);
        % Feature 8: Percent of RR larger than 20ms
        pNN20 = NaN(numberOfWins,2);
        % Feature 9: Percent of RR larger than 50ms
        pNN50 = NaN(numberOfWins,2);
        % Feature 10: Baevsky's Stress Index
        SI = NaN(numberOfWins,2);
        % Feature 11: Low Frecuency Power
        LF = NaN(numberOfWins,2);
        % Feature 12: High Frecuency Power
        HF = NaN(numberOfWins,2);
        % Feature 13: Low Frecuency to High Frecuency Rate
        LFHF = NaN(numberOfWins,2);
        
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

cane_umbral_array=[];
empatica_umbral_array=[];
zeros_cane_array=[];
zeros_cane_time=[];
zeros_empatica_array=[];
zeros_empatica_time=[];
cane_umbral_time=[];
empatica_umbral_time=[];




cane_bueno_array=[];
cane_bueno_time=[];

indices_max_pos=[];
indices_max_neg=[];
indices_min_pos=[];
indices_min_neg=[];
indices_true=[];
muestra_mala=0;

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

            % Time Variables
            % Calculamos los cruces por cero (semiperiodo)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Quizas haya que sustituir por el metodo avanzado
            zeroCane=[];
            zeroEmpatica=[];
            cane_umbral=mean(ppgCaneWin);
            empatica_umbral=mean(ppgCaneWin);
            for idx = 2:length(ppgCaneWin)-1
                % Paso por cero
                if (ppgCaneWin(idx) > cane_umbral && ppgCaneWin(idx+1) < cane_umbral) || (ppgCaneWin(idx) < cane_umbral && ppgCaneWin(idx+1) > cane_umbral)
                    zeroCane = cat(1,zeroCane,idx);
                end
% Maximo
if ppgCaneWin(idx-1) < ppgCaneWin(idx) && ppgCaneWin(idx) > ppgCaneWin(idx+1)
    if ppgCaneWin(idx) > cane_umbral
        indices_max_pos = cat(1,indices_max_pos,idx);
        if ~muestra_mala
            indices_true = cat(1,indices_true,idx);
        end
        muestra_mala=0;
    else
        indices_max_neg = cat(1,indices_max_neg,idx);
        muestra_mala = 1;
    end
end
% Minimo
if ppgCaneWin(idx-1) > ppgCaneWin(idx) && ppgCaneWin(idx) < ppgCaneWin(idx+1)

        indices_min_pos = cat(1,indices_min_pos,idx);
    if ppgCaneWin(idx) > cane_umbral
        muestra_mala = 1;
    else
        indices_min_neg = cat(1,indices_min_neg,idx);
        if ~muestra_mala
            indices_true = cat(1,indices_true,idx);
        end
        muestra_mala=0;
    end
end
plot(ppgCaneWin)
            end
            for idx = 2:length(ppgEmpaticaWin)-1
                % Paso por cero
                if (ppgEmpaticaWin(idx) > empatica_umbral && ppgEmpaticaWin(idx+1) < empatica_umbral) || (ppgEmpaticaWin(idx) < empatica_umbral && ppgEmpaticaWin(idx+1) > empatica_umbral)
                    zeroEmpatica = cat(1,zeroEmpatica,idx);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("_________________________")

            % R to R interval individual
            timestampCane       = caneDataExperiment.(experimentNames{i}).PPG.TimeStampPosix;
            timestampEmpatica   = empaticaDataExperiment.(experimentNames{i}).PPG.TimeStampPosix;
            timestampCaneWin        = timestampCane(j:j+winLengthPPG-1);
            timestampEmpaticaWin    = timestampEmpatica(j:j+winLengthPPG-1);
            RR_Cane = 2*(timestampCaneWin(zeroCane(2:end))         -timestampCaneWin(zeroCane(1:end-1)));
            RR_Empatica = 2*(timestampEmpaticaWin(zeroEmpatica(2:end)) -timestampEmpaticaWin(zeroEmpatica(1:end-1)));
         
cane_bueno_array=[cane_bueno_array; ppgCaneWin(indices_min_pos)];
cane_bueno_time=[cane_bueno_time; timestampCaneWin(indices_min_pos)];

cane_umbral_array=[cane_umbral_array cane_umbral];
empatica_umbral_array=[empatica_umbral_array empatica_umbral];
zeros_cane_array=[zeros_cane_array; ppgCaneWin(zeroCane)];
zeros_cane_time=[zeros_cane_time; timestampCaneWin(zeroCane)];
zeros_empatica_array=[zeros_empatica_array; ppgEmpaticaWin(zeroEmpatica)];
zeros_empatica_time=[zeros_empatica_time; timestampEmpaticaWin(zeroEmpatica)];
cane_umbral_time=[cane_umbral_time timestampCaneWin(end/2)];
empatica_umbral_time=[empatica_umbral_time timestampEmpaticaWin(end/2)];

            % Feature 3: R to R interval
            MeanRR(countWin,1) = mean(RR_Cane);
            MeanRR(countWin,2) = mean(RR_Empatica);
            % Feature 4: Heart Rate
            HR(countWin,1) = 60/MeanRR(countWin,1);
            HR(countWin,2) = 60/MeanRR(countWin,2);
            % Feature 5: Standar Deviation of RR
            SDNN(countWin,1) = sqrt(sum((RR_Cane-MeanRR(countWin,1)))^2/(winLengthPPG-1));
            SDNN(countWin,2) = sqrt(sum((RR_Empatica-MeanRR(countWin,2)))^2/(winLengthPPG-1));
            % Feature 6: Root Mean Square of Sucesives RR interval diference
            RMSSD(countWin,1) = sqrt(sum((RR_Cane(2:end)-RR_Cane(1:end-1)).^2)/(winLengthPPG-1));
            RMSSD(countWin,2) = sqrt(sum((RR_Empatica(2:end)-RR_Empatica(1:end-1)).^2)/(winLengthPPG-1));
            % Feature 7: Sucesive Differences Standardar Deviation
            SucDif_Cane=abs(RR_Cane(2:end)-RR_Cane(1:end-1));
            SDSD(countWin,1) = sqrt(sum((SucDif_Cane-mean(SucDif_Cane)).^2)/(winLengthPPG-2));
            SucDif_Empatica=abs(RR_Empatica(2:end)-RR_Empatica(1:end-1));
            SDSD(countWin,2) = sqrt(sum((SucDif_Empatica-mean(SucDif_Empatica)).^2)/(winLengthPPG-2));
            % Feature 8: Percent of RR larger than 20ms
            pNN20(countWin,1) = sum(RR_Cane>20)/winLengthPPG*100;
            pNN20(countWin,2) = sum(RR_Empatica>20)/winLengthPPG*100;
            % Feature 9: Percent of RR larger than 50ms
            pNN50(countWin,1) = sum(RR_Cane>50)/winLengthPPG*100;
            pNN50(countWin,2) = sum(RR_Empatica>50)/winLengthPPG*100;
            % Feature 10: Baevsky's Stress Index
            % Calculamos el hsitograma primero: incremeto de 50ms
            t_intervalos=0.05;
            if(length(RR_Cane)>2)
                caneintervalos=(min(RR_Cane):t_intervalos:max(RR_Cane))+t_intervalos/2;% Puntos centrales de los rangos
                canehist = histcounts(RR_Cane, caneintervalos);
                [caneAMo_pre, canehist_idx] = max(canehist);
                caneMo = caneintervalos(canehist_idx);
                caneAMo = caneAMo_pre/winLengthPPG;
                caneMxDMn = max(RR_Cane) - min(RR_Cane);
                SI(countWin,1) = (caneAMo/(2*caneMo*caneMxDMn))*100; % Para el bastón 
            else
                SI(countWin,1) = 0;
            end
            if(length(RR_Empatica)>2)
                empaticaintervalos=(min(RR_Empatica):t_intervalos:max(RR_Empatica))+t_intervalos/2;% Puntos centrales de los rangos
                empaticahist = histcounts(RR_Empatica, empaticaintervalos);
                [empaticaAMo_pre, empaticahistidx] = max(empaticahist);
                empaticaMo = (empaticaintervalos(empaticahistidx) + empaticaintervalos(empaticahistidx+1)) / 2;
                empaticaAMo = empaticaAMo_pre/winLengthPPG;
                empaticaMxDMn = max(RR_Empatica) - min(RR_Empatica);
                SI(countWin,2) = (empaticaAMo/(2*empaticaMo*empaticaMxDMn))*100; % Para la empatica
            else
                SI(countWin,2) = 0;
            end

            % Frecuency Variables
            % PPG
            VL_cutoff = 0.04 / (Fs_PPG / 2); % VL (0-0.04Hz)
            L_cutoff = 0.15 / (Fs_PPG / 2);  % L (0.04-0.15Hz)
            H_cutoff = 0.4 / (Fs_PPG / 2);  % H (0.15-0.4Hz)
            % Diseño del filtro Butterworth de orden 4
            [b_low, a_low] = butter(4, [VL_cutoff, L_cutoff], 'bandpass');
            [b_high, a_high] = butter(4, [L_cutoff, H_cutoff], 'bandpass');
            % Feature 11: Low Frecuency Power
            LF(countWin,1) = mean(abs(filter(b_low, a_low, ppgCaneWin)).^2) / Fs_PPG;
            LF(countWin,2) = mean(abs(filter(b_low, a_low, ppgEmpaticaWin)).^2) / Fs_PPG;
            % Feature 12: High Frecuency Power
            HF(countWin,1) = mean(abs(filter(b_high, a_high, ppgCaneWin)).^2) / Fs_PPG;
            HF(countWin,2) = mean(abs(filter(b_high, a_high, ppgEmpaticaWin)).^2) / Fs_PPG;
            % Feature 13: Low Frecuency to High Frecuency Power Rate
            LFHF(countWin,1) = LF(countWin,1)/HF(countWin,1);
            LFHF(countWin,2) = LF(countWin,2)/HF(countWin,2);
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
            MeanRR(:,1), ...
            HR(:,1), ...
            SDNN(:,1), ...
            RMSSD(:,1), ...
            SDSD(:,1), ...
            pNN20(:,1), ...
            pNN50(:,1), ...
            SI(:,1), ...
            LF(:,1), ...
            HF(:,1), ...
            LFHF(:,1), ...
            ...
            subjectID, ...
            experimentID, ...
            label, ...
            'VariableNames',{'ppgMaxValue','ppgMinValue','gsrMinValue','tonicMaxValue','phasicSTD', ...
            ... % <------- Completar con todas las caracteristicas
            'MeanRR','HR','SDNN','RMSSD','SDSD','pNN20','pNN50','SI','LF','HF','LFHF',...
            ...
            'subjectID', 'experimentID', 'label'});
    
        empaticaFeatureTable_temp = table( ...
            ppgMaxValue(:,2), ...
            ppgMinValue(:,2), ...
            gsrMinValue(:,2), ...
            tonicMaxValue(:,2), ...
            phasicSTD(:,2), ...
            ... % <------- Completar con todas las caracteristicas
            MeanRR(:,2), ...
            HR(:,2), ...
            SDNN(:,2), ...
            RMSSD(:,2), ...
            SDSD(:,2), ...
            pNN20(:,2), ...
            pNN50(:,2), ...
            SI(:,2), ...
            LF(:,2), ...
            HF(:,2), ...
            LFHF(:,2), ...
            ...
            subjectID, ...
            experimentID, ...
            label, ...
            'VariableNames',{'ppgMaxValue','ppgMinValue','gsrMinValue','tonicMaxValue','phasicSTD', ...
            ... % <------- Completar con todas las caracteristicas
            'MeanRR','HR','SDNN','RMSSD','SDSD','pNN20','pNN50','SI','LF','HF','LFHF',...
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