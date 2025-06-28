function dataStruct = loadAllData(parado)


    if(parado==1)
        basePath = "C:\Users\Usuario\Escritorio\TFG\TODOS\Participantes\Toma de datos (hasta 24-05-25)";
        carpeta_data='signal';
    elseif (parado==2)
        basePath = "C:\Users\Usuario\Escritorio\TFG\TODOS\Participantes\Toma de datos - filtrado (11-06-25)";
        carpeta_data='filt_signal';
    else
        basePath = "C:\Users\Usuario\Escritorio\TFG\TODOS\Participantes\Toma de datos - Experimento cinta (13-06-25)";
        carpeta_data='filt_signal';
    end
    % Inicializa el struct de salida
    dataStruct = struct();

    % Obtener carpetas que empiezan con 'S' y contienen números
    dirInfo = dir(fullfile(basePath, 'S*'));
    dirInfo = dirInfo([dirInfo.isdir]);  % Solo carpetas

    for k = 1:length(dirInfo)
        folderName = dirInfo(k).name;

        % Extraer número de participante
        tokens = regexp(folderName, '^S(\d+)$', 'tokens');
        if isempty(tokens)
            continue;
        end
        participantID = tokens{1}{1};
        structField = ['S', participantID];

        % Rutas
        signalFolder = fullfile(basePath, folderName, carpeta_data);


        % Nombres esperados
        caneFile = fullfile(signalFolder, sprintf('S%s_caneDataExperiment.mat', participantID));
        empaticaFile = fullfile(signalFolder, sprintf('S%s_empaticaDataExperiment.mat', participantID));
        audioFile = fullfile(fullfile(basePath, folderName), sprintf('S%s_audioEventVector.mat', participantID));

        % Inicializa subestructura
        substruct = struct();
        % Cargar caneData
        if isfile(caneFile)
            temp = load(caneFile);
            fn = fieldnames(temp);
            substruct.caneData = temp.(fn{1});  % Asume que hay solo una variable en el .mat
        end
        
        % Cargar empaticaData
        if isfile(empaticaFile)
            temp = load(empaticaFile);
            fn = fieldnames(temp);
            substruct.empaticaData = temp.(fn{1});
        end
        
        % Cargar audioEventVector
        if isfile(audioFile)
            temp = load(audioFile);
            fn = fieldnames(temp);
            substruct.audioEventVector = temp.(fn{1});
        end
        
        % Asignar al struct principal
        dataStruct.(structField) = substruct;
    end
end
