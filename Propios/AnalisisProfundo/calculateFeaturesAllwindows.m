function [tableStruct] = calculateFeaturesAllwindows(data,tiempo_ventana, idx_estimulo, participantesRechazados, AllSignal,tipoExp)
%INPUTS
%{
    data = datos procesados por la función 'loadAllData.m'
%}
%OUTPUTS
%{
    tableempaticapreevent = tabla de características de la empatica de las ventanas de pre_evento
    tableempaticapostevent = tabla de características de la empatica de las ventanas de post_evento
    tablecanepreevent = tabla de características del bastón de las ventanas de pre_evento
    tablecanepostevent = tabla de características del bastón de las ventanas de post_evento
    tableSYSPpreevent
    tableSYSPpostevet
%}
%INICIALIZACIONES
rechazo = participantesRechazados;
%rechazo = [];
tableempaticapreevent = table();
tableempaticapostevent = table();
tablecanepreevent = table();
tablecanepostevent = table();
tableSYSPpreevent = table();
tableSYSPpostevent = table();
participante = fieldnames(data);


for i = 1:length(participante)
    nombre = participante{i};
    % Extraer número del participante
    num_part = str2double(extractAfter(nombre, 'S'));    
    % Verificar si está en la lista de exclusión
    if ismember(num_part, rechazo)
        continue  % Saltar este participante
    end
    % Obtener subestructura de datos
    data_struct = data.(nombre);
    [preevent_features_cane, preevent_features_empatica, preevent_SYSPfeatures_empatica, postevent_features_cane, postevent_features_empatica, postevent_SYSPfeatures_empatica] = calculateFeaturesM (data_struct, num_part, tiempo_ventana, idx_estimulo, AllSignal,tipoExp);
    tableempaticapreevent = [tableempaticapreevent;preevent_features_empatica];
    tableempaticapostevent = [tableempaticapostevent;postevent_features_empatica];
    tablecanepreevent = [tablecanepreevent;preevent_features_cane];
    tablecanepostevent = [tablecanepostevent;postevent_features_cane];
    tableSYSPpreevent = [tableSYSPpreevent;preevent_SYSPfeatures_empatica];
    tableSYSPpostevent = [tableSYSPpostevent;postevent_SYSPfeatures_empatica];

    disp([nombre,'Hecho']);
end

tableStruct = struct( ...
    'empatica_pre', tableempaticapreevent, ...
    'empatica_post', tableempaticapostevent, ...
    'cane_pre', tablecanepreevent, ...
    'cane_post', tablecanepostevent, ...
    'SYSP_pre', tableSYSPpreevent, ...
    'SYSP_post', tableSYSPpostevent ...
);

end