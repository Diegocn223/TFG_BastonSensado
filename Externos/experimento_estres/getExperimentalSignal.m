function [caneDataExperiment, empaticaDataExperiment] = getExperimentalSignal(empaticaData, caneData, protocolTime_datetime)
% getExperimentalSignal Extrae y organiza los datos de CaneSense y Empatica para cada experimento
%   Esta función toma los datos de CaneSense (`caneData`) y Empatica (`empaticaData`), así como el vector `protocolTime_datetime` que marca el inicio y fin de cada uno de los cuatro experimentos de la sesión.
%   La función organiza y extrae los datos correspondientes a cada experimento (Relax 1, Estrés Música, Relax 2, Estrés Aritmética) y los almacena en dos estructuras de salida: `caneDataExperiment` y `empaticaDataExperiment`.
%
%   Ejemplo de uso:
%   protocolTime_datetime = [
%       datetime('16:30:35'), datetime('16:30:35') + minutes(4);  % Relax 1
%       datetime('16:36:12'), datetime('16:36:12') + minutes(3);  % Estrés Música
%       datetime('16:42:01'), datetime('16:42:01') + minutes(4);  % Relax 2
%       datetime('16:49:37'), datetime('16:49:37') + minutes(4)   % Estrés Aritmética
%   ];
%
%   [caneDataExperiment, empaticaDataExperiment] = getExperimentalSignal(empaticaData, caneData, protocolTime_datetime, participant_ID);
%   En este ejemplo, `empaticaData` y `caneData` son estructuras de datos cargadas previamente, y la función devuelve dos estructuras: `caneDataExperiment` y `empaticaDataExperiment` con los datos organizados por cada experimento.
%
% Entradas:
%   - empaticaData: Estructura de datos del reloj Empatica que contiene diferentes tipos de señales (GSR, TONIC, PHASIC, PPG, SYSTP, TAGS).
%   - caneData: Tabla de datos de la señal del bastón inteligente (CaneSense App), con las fechas de los registros.
%   - protocolTime_datetime: Array de 4 filas y 2 columnas de tiempos en formato datetime, donde la primera columna es la fecha de inicio y la segunda columna es la fecha de fin de cada uno de los experimentos (Relax 1, Estrés Música, Relax 2, Estrés Aritmética).
%
% Salidas:
%   - caneDataExperiment: Estructura que contiene los datos del bastón para cada uno de los experimentos.
%   - empaticaDataExperiment: Estructura que contiene los datos de Empatica para cada uno de los experimentos.

% Verificación de que caneData es una tabla y tiene al menos una columna
    if ~istable(caneData)
        error('caneData debe ser una tabla.');
    elseif width(caneData) < 1
        error('caneData debe tener al menos una columna.');
    end

    % Verificación de que empaticaData es un struct y tiene al menos un campo
    if ~isstruct(empaticaData)
        error('empaticaData debe ser un struct.');
    elseif numel(fieldnames(empaticaData)) < 1
        error('empaticaData debe tener al menos un campo.');
    end
    
    % Verificación de que protocolTime_datetime es un array de 4 filas y 2 columnas de tipo datetime
    if ~ismatrix(protocolTime_datetime) || size(protocolTime_datetime, 1) ~= 4 || size(protocolTime_datetime, 2) ~= 2
        error('protocolTime_datetime debe ser un array de 4 filas y 2 columnas.');
    elseif ~all(isdatetime(protocolTime_datetime(:)))
        error('Todos los elementos de protocolTime_datetime deben ser objetos datetime.');
    end


    % Creamos un struct para poner los datos de la Empatica semejante a
    % caneData
    empaticaDataPreprocess = struct();
    empaticaDataFieldNames = {'EDA', 'TONIC', 'PHASIC', 'BVP', 'SystP', 'Tags'};
    empaticaDataPreprocessFieldNames = {'GSR', 'TONIC', 'PHASIC', 'PPG', 'SYSTP', 'TAGS'};
    
    % Crear las tablas para cada campo (excepto TAGS que es diferente)
    for i = 1:numel(empaticaDataFieldNames)
        if strcmp(empaticaDataPreprocessFieldNames{i}, 'TAGS')
            if ~isempty(empaticaData.Tags{1})
                empaticaDataPreprocess.(empaticaDataPreprocessFieldNames{i}) = table(empaticaData.Tags{1}.Tag, str2double(empaticaData.Tags{1}.UTC_FileTime{1}), ...
                    datetime(str2double(empaticaData.Tags{1}.UTC_FileTime{1}), 'ConvertFrom', 'posixtime', 'TimeZone', 'Europe/Madrid', 'Format', 'dd-MM-yyyy HH:mm:ss'), ...
                    'VariableNames', {'Tag', 'TimeStampPosix', 'TimeStampDate'});
            else
                empaticaDataPreprocess.(empaticaDataPreprocessFieldNames{i}) = table([],[],[], ...
                    'VariableNames', {'Tag', 'TimeStampPosix', 'TimeStampDate'});
            end
        else
            empaticaDataPreprocess.(empaticaDataPreprocessFieldNames{i}) = table(empaticaData.(empaticaDataFieldNames{i}){1}.data, ...
                empaticaData.(empaticaDataFieldNames{i}){1}.unix_time ./ 1e6, empaticaData.(empaticaDataFieldNames{i}){1}.local_time, ...
                'VariableNames', {'data', 'TimeStampPosix', 'TimeStampDate'});
        end
    end

    % Dado que empaticaData.SystP.unix_time está en nanosegundos (no en
    % micro), lo volvemos a dividir entre 1000.
    empaticaDataPreprocess.SYSP.TimeStampPosix = empaticaDataPreprocess.SYSTP.TimeStampPosix ./ 1e3;
    
    % A continuación vamos a construir las ventanas de tiempo de los
    % experimentos realizados conociendo protocolTime_datetime

    relax1Datetime_init = protocolTime_datetime(1,1);
    relax1Datetime_end = protocolTime_datetime(1,2);

    stressMusicDatetime_init = protocolTime_datetime(2,1);
    stressMusicDatetime_end = protocolTime_datetime(2,2);

    relax2Datetime_init = protocolTime_datetime(3,1);
    relax2Datetime_end = protocolTime_datetime(3,2);

    stressArithmeticDatetime_init = protocolTime_datetime(4,1);
    stressArithmeticDatetime_end = protocolTime_datetime(4,2);
    

    %% Buscamos los índices de datos del BASTÓN para cada experimento
    % Relax 1
    [index_start, index_end] = findDatetimeIndex(relax1Datetime_init, relax1Datetime_end,  caneData.TimeStampDate);
    caneDataRelax1 = caneData(index_start:(index_end+1), :);

    % Estrés Música   
    [index_start, index_end] = findDatetimeIndex(stressMusicDatetime_init, stressMusicDatetime_end,  caneData.TimeStampDate);
    caneDataMusic = caneData(index_start:(index_end+1), :);

    % Relax 2
    [index_start, index_end] = findDatetimeIndex(relax2Datetime_init, relax2Datetime_end,  caneData.TimeStampDate);
    caneDataRelax2 = caneData(index_start:(index_end+1), :);

    % Estrés Aritmética  
    [index_start, index_end] = findDatetimeIndex(stressArithmeticDatetime_init, stressArithmeticDatetime_end,  caneData.TimeStampDate);
    caneDataArithmetic = caneData(index_start:(index_end+1), :);

    caneDataExperiment_table = struct();
    caneDataExperiment_table.Relax1 = caneDataRelax1;
    caneDataExperiment_table.Music = caneDataMusic;
    caneDataExperiment_table.Relax2 = caneDataRelax2;
    caneDataExperiment_table.Arithmetic = caneDataArithmetic;

    %% Buscamos los índices de datos de EMPATICA para cada experimento
    % Relax 1
        % GSR, TONIC, PHASIC
        [index_start, index_end] = findDatetimeIndex(relax1Datetime_init, relax1Datetime_end,  empaticaDataPreprocess.GSR.TimeStampDate);
        empaticaGSRDataRelax1 = empaticaDataPreprocess.GSR(index_start:(index_end+1),:);
        empaticaTONICDataRelax1 = empaticaDataPreprocess.TONIC(index_start:(index_end+1),:);
        empaticaPHASICDataRelax1 = empaticaDataPreprocess.PHASIC(index_start:(index_end+1),:);
        % PPG
        [index_start, index_end] = findDatetimeIndex(relax1Datetime_init, relax1Datetime_end,  empaticaDataPreprocess.PPG.TimeStampDate);
        empaticaPPGDataRelax1 = empaticaDataPreprocess.PPG(index_start:(index_end+1),:);
        % SYSTP
        [index_start, index_end] = findDatetimeIndex(relax1Datetime_init, relax1Datetime_end,  empaticaDataPreprocess.SYSTP.TimeStampDate);
        empaticaSYSTPDataRelax1 = empaticaDataPreprocess.SYSTP(index_start:(index_end+1),:);
        % TAGS
        if ~isempty(empaticaDataPreprocess.TAGS.TimeStampDate)
            [index_start, index_end] = findDatetimeIndex(relax1Datetime_init, relax1Datetime_end,  empaticaDataPreprocess.TAGS.TimeStampDate);
            empaticaTAGSDataRelax1 = empaticaDataPreprocess.TAGS(index_start:(index_end+1),:);
        else
            empaticaTAGSDataRelax1 = [];
        end
        
        relax1 = struct();
        relax1.GSR = empaticaGSRDataRelax1;
        relax1.TONIC = empaticaTONICDataRelax1;
        relax1.PHASIC = empaticaPHASICDataRelax1;
        relax1.PPG = empaticaPPGDataRelax1;
        relax1.SYSTP = empaticaSYSTPDataRelax1;
        relax1.TAGS = empaticaTAGSDataRelax1;
        
    % Estrés Música
        % GSR
        [index_start, index_end] = findDatetimeIndex(stressMusicDatetime_init, stressMusicDatetime_end,  empaticaDataPreprocess.GSR.TimeStampDate);
        empaticaGSRDataStressMusic = empaticaDataPreprocess.GSR(index_start:(index_end+1),:);
        empaticaTONICDataStressMusic = empaticaDataPreprocess.TONIC(index_start:(index_end+1),:);
        empaticaPHASICDataStressMusic = empaticaDataPreprocess.PHASIC(index_start:(index_end+1),:);
        % PPG
        [index_start, index_end] = findDatetimeIndex(stressMusicDatetime_init, stressMusicDatetime_end,  empaticaDataPreprocess.PPG.TimeStampDate);
        empaticaPPGDataStressMusic = empaticaDataPreprocess.PPG(index_start:(index_end+1),:);
        % SYSTP
        [index_start, index_end] = findDatetimeIndex(stressMusicDatetime_init, stressMusicDatetime_end,  empaticaDataPreprocess.SYSTP.TimeStampDate);
        empaticaSYSTPDataStressMusic = empaticaDataPreprocess.SYSTP(index_start:(index_end+1),:);
        % TAGS
        if ~isempty(empaticaDataPreprocess.TAGS.TimeStampDate)
            [index_start, index_end] = findDatetimeIndex(stressMusicDatetime_init, stressMusicDatetime_end,  empaticaDataPreprocess.TAGS.TimeStampDate);
            empaticaTAGSDataStressMusic = empaticaDataPreprocess.TAGS(index_start:(index_end+1),:);
        else
            empaticaTAGSDataStressMusic = [];
        end

        music = struct();
        music.GSR = empaticaGSRDataStressMusic;
        music.TONIC = empaticaTONICDataStressMusic;
        music.PHASIC = empaticaPHASICDataStressMusic;
        music.PPG = empaticaPPGDataStressMusic;
        music.SYSTP = empaticaSYSTPDataStressMusic;
        music.TAGS = empaticaTAGSDataStressMusic;

    % Relax 2
        % GSR
        [index_start, index_end] = findDatetimeIndex(relax2Datetime_init, relax2Datetime_end,  empaticaDataPreprocess.GSR.TimeStampDate);
        empaticaGSRDataRelax2 = empaticaDataPreprocess.GSR(index_start:(index_end+1),:);
        empaticaTONICDataRelax2 = empaticaDataPreprocess.TONIC(index_start:(index_end+1),:);
        empaticaPHASICDataRelax2 = empaticaDataPreprocess.PHASIC(index_start:(index_end+1),:);
        % PPG
        [index_start, index_end] = findDatetimeIndex(relax2Datetime_init, relax2Datetime_end,  empaticaDataPreprocess.PPG.TimeStampDate);
        empaticaPPGDataStressRelax2 = empaticaDataPreprocess.PPG(index_start:(index_end+1),:);
        % SYSTP
        [index_start, index_end] = findDatetimeIndex(relax2Datetime_init, relax2Datetime_end,  empaticaDataPreprocess.SYSTP.TimeStampDate);
        empaticaSYSTPDataRelax2 = empaticaDataPreprocess.SYSTP(index_start:(index_end+1),:);
        % TAGS
        if ~isempty(empaticaDataPreprocess.TAGS.TimeStampDate)
            [index_start, index_end] = findDatetimeIndex(relax2Datetime_init, relax2Datetime_end,  empaticaDataPreprocess.TAGS.TimeStampDate);
            empaticaTAGSDataRelax2 = empaticaDataPreprocess.TAGS(index_start:(index_end+1),:);
        else
            empaticaTAGSDataRelax2 = [];
        end

        relax2 = struct();
        relax2.GSR = empaticaGSRDataRelax2;
        relax2.TONIC = empaticaTONICDataRelax2;
        relax2.PHASIC = empaticaPHASICDataRelax2;
        relax2.PPG = empaticaPPGDataStressRelax2;
        relax2.SYSTP = empaticaSYSTPDataRelax2;
        relax2.TAGS = empaticaTAGSDataRelax2;

    % Estrés Aritmética
        % GSR
        [index_start, index_end] = findDatetimeIndex(stressArithmeticDatetime_init, stressArithmeticDatetime_end,  empaticaDataPreprocess.GSR.TimeStampDate);
        empaticaGSRDataStressArithmetic = empaticaDataPreprocess.GSR(index_start:(index_end+1),:);
        empaticaTONICStressArithmetic = empaticaDataPreprocess.TONIC(index_start:(index_end+1),:);
        empaticaPHASICStressArithmetic = empaticaDataPreprocess.PHASIC(index_start:(index_end+1),:);
        % PPG
        [index_start, index_end] = findDatetimeIndex(stressArithmeticDatetime_init, stressArithmeticDatetime_end,  empaticaDataPreprocess.PPG.TimeStampDate);
        empaticaPPGDataStressArithmetic = empaticaDataPreprocess.PPG(index_start:(index_end+1),:);
        % SYSTP
        [index_start, index_end] = findDatetimeIndex(stressArithmeticDatetime_init, stressArithmeticDatetime_end,  empaticaDataPreprocess.SYSTP.TimeStampDate);
        empaticaSYSTPStressArithmetic = empaticaDataPreprocess.SYSTP(index_start:(index_end+1),:);
        % TAGS
        if ~isempty(empaticaDataPreprocess.TAGS.TimeStampDate)
            [index_start, index_end] = findDatetimeIndex(stressArithmeticDatetime_init, stressArithmeticDatetime_end,  empaticaDataPreprocess.TAGS.TimeStampDate);
            empaticaTAGSStressArithmetic = empaticaDataPreprocess.TAGS(index_start:(index_end+1),:);
        else
            empaticaTAGSStressArithmetic = [];
        end

        arithmetic = struct();
        arithmetic.GSR = empaticaGSRDataStressArithmetic;
        arithmetic.TONIC = empaticaTONICStressArithmetic;
        arithmetic.PHASIC = empaticaPHASICStressArithmetic;
        arithmetic.PPG = empaticaPPGDataStressArithmetic;
        arithmetic.SYSTP = empaticaSYSTPStressArithmetic;
        arithmetic.TAGS = empaticaTAGSStressArithmetic;

        empaticaDataExperiment = struct();
        empaticaDataExperiment.Relax1 = relax1;
        empaticaDataExperiment.Music = music;
        empaticaDataExperiment.Relax2 = relax2;
        empaticaDataExperiment.Arithmetic = arithmetic;

        %% Modificamos la estructura de datos del Bastón
        % Para que tenga una forma similar a la de Empatica
        GSR_relax1 = [caneDataExperiment_table.Relax1(:,'GSR'), caneDataExperiment_table.Relax1(:,'TimeStampPosix'), caneDataExperiment_table.Relax1(:,'TimeStampDate')];
        TONIC_relax1 = [caneDataExperiment_table.Relax1(:,'TONIC'), caneDataExperiment_table.Relax1(:,'TimeStampPosix'), caneDataExperiment_table.Relax1(:,'TimeStampDate')];
        PHASIC_relax1 = [caneDataExperiment_table.Relax1(:,'PHASIC'), caneDataExperiment_table.Relax1(:,'TimeStampPosix'), caneDataExperiment_table.Relax1(:,'TimeStampDate')];
        PPG_relax1 = [caneDataExperiment_table.Relax1(:,'PPG'), caneDataExperiment_table.Relax1(:,'TimeStampPosix'), caneDataExperiment_table.Relax1(:,'TimeStampDate')];
        FSR_relax1 = [caneDataExperiment_table.Relax1(:,'FSR'), caneDataExperiment_table.Relax1(:,'TimeStampPosix'), caneDataExperiment_table.Relax1(:,'TimeStampDate')];
        
        GSR_relax1.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};
        TONIC_relax1.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};
        PHASIC_relax1.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};
        PPG_relax1.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};
        FSR_relax1.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};

        Relax1 = struct();
        Relax1.GSR = GSR_relax1;
        Relax1.TONIC = TONIC_relax1;
        Relax1.PHASIC = PHASIC_relax1;
        Relax1.PPG = PPG_relax1;
        Relax1.FSR = FSR_relax1;

        GSR_music = [caneDataExperiment_table.Music(:,'GSR'), caneDataExperiment_table.Music(:,'TimeStampPosix'), caneDataExperiment_table.Music(:,'TimeStampDate')];
        TONIC_music = [caneDataExperiment_table.Music(:,'TONIC'), caneDataExperiment_table.Music(:,'TimeStampPosix'), caneDataExperiment_table.Music(:,'TimeStampDate')];
        PHASIC_music = [caneDataExperiment_table.Music(:,'PHASIC'), caneDataExperiment_table.Music(:,'TimeStampPosix'), caneDataExperiment_table.Music(:,'TimeStampDate')];
        PPG_music = [caneDataExperiment_table.Music(:,'PPG'), caneDataExperiment_table.Music(:,'TimeStampPosix'), caneDataExperiment_table.Music(:,'TimeStampDate')];
        FSR_music = [caneDataExperiment_table.Music(:,'FSR'), caneDataExperiment_table.Music(:,'TimeStampPosix'), caneDataExperiment_table.Music(:,'TimeStampDate')];

        GSR_music.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};
        TONIC_music.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};
        PHASIC_music.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};
        PPG_music.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};
        FSR_music.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};

        Music = struct();
        Music.GSR = GSR_music;
        Music.TONIC = TONIC_music;
        Music.PHASIC = PHASIC_music;
        Music.PPG = PPG_music;
        Music.FSR = FSR_music;

        GSR_relax2 = [caneDataExperiment_table.Relax2(:,'GSR'), caneDataExperiment_table.Relax2(:,'TimeStampPosix'), caneDataExperiment_table.Relax2(:,'TimeStampDate')];
        TONIC_relax2 = [caneDataExperiment_table.Relax2(:,'TONIC'), caneDataExperiment_table.Relax2(:,'TimeStampPosix'), caneDataExperiment_table.Relax2(:,'TimeStampDate')];
        PHASIC_relax2 = [caneDataExperiment_table.Relax2(:,'PHASIC'), caneDataExperiment_table.Relax2(:,'TimeStampPosix'), caneDataExperiment_table.Relax2(:,'TimeStampDate')];
        PPG_relax2 = [caneDataExperiment_table.Relax2(:,'PPG'), caneDataExperiment_table.Relax2(:,'TimeStampPosix'), caneDataExperiment_table.Relax2(:,'TimeStampDate')];
        FSR_relax2 = [caneDataExperiment_table.Relax2(:,'FSR'), caneDataExperiment_table.Relax2(:,'TimeStampPosix'), caneDataExperiment_table.Relax2(:,'TimeStampDate')];

        GSR_relax2.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};
        TONIC_relax2.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};
        PHASIC_relax2.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};
        PPG_relax2.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};
        FSR_relax2.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};

        Relax2 = struct();
        Relax2.GSR = GSR_relax2;
        Relax2.TONIC = TONIC_relax2;
        Relax2.PHASIC = PHASIC_relax2;
        Relax2.PPG = PPG_relax2;
        Relax2.FSR = FSR_relax2;

        GSR_arithmetic = [caneDataExperiment_table.Arithmetic(:,'GSR'), caneDataExperiment_table.Arithmetic(:,'TimeStampPosix'), caneDataExperiment_table.Arithmetic(:,'TimeStampDate')];
        TONIC_arithmetic = [caneDataExperiment_table.Arithmetic(:,'TONIC'), caneDataExperiment_table.Arithmetic(:,'TimeStampPosix'), caneDataExperiment_table.Arithmetic(:,'TimeStampDate')];
        PHASIC_arithmetic = [caneDataExperiment_table.Arithmetic(:,'PHASIC'), caneDataExperiment_table.Arithmetic(:,'TimeStampPosix'), caneDataExperiment_table.Arithmetic(:,'TimeStampDate')];
        PPG_arithmetic = [caneDataExperiment_table.Arithmetic(:,'PPG'), caneDataExperiment_table.Arithmetic(:,'TimeStampPosix'), caneDataExperiment_table.Arithmetic(:,'TimeStampDate')];
        FSR_arithmetic = [caneDataExperiment_table.Arithmetic(:,'FSR'), caneDataExperiment_table.Arithmetic(:,'TimeStampPosix'), caneDataExperiment_table.Arithmetic(:,'TimeStampDate')];

        GSR_arithmetic.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};
        TONIC_arithmetic.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};
        PHASIC_arithmetic.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};
        PPG_arithmetic.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};
        FSR_arithmetic.Properties.VariableNames = {'data', 'TimeStampPosix', 'TimeStampDate'};

        Arithmetic = struct();
        Arithmetic.GSR = GSR_arithmetic;
        Arithmetic.TONIC = TONIC_arithmetic;
        Arithmetic.PHASIC = PHASIC_arithmetic;
        Arithmetic.PPG = PPG_arithmetic;
        Arithmetic.FSR = FSR_arithmetic;

        caneDataExperiment = struct();
        caneDataExperiment.Relax1 = Relax1;
        caneDataExperiment.Music = Music;
        caneDataExperiment.Relax2 = Relax2;
        caneDataExperiment.Arithmetic = Arithmetic;
        
end
