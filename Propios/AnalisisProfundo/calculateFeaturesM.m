function [preevent_features_cane, preevent_features_empatica, preevent_SYSPfeatures_empatica, postevent_features_cane, postevent_features_empatica, postevent_SYSPfeatures_empatica] = calculateFeaturesM (data_struct, indice, tiempo_ventana, idx_estimulo, AllSignal, tipoExp)
%INPUTS
%{
    -limit_preevent = matriz de 2 columnas y 6 filas que nos da la hora en datetime del inicio y fin de la ventana de pre-evento
    -preevent_features = matriz con las características para las 6 ventanas de pre-event
    -limit_postevent = matriz de 2 columnas y 6 filas que nos da la hora en datetime del inicio y fin de la ventana de post-evento
    -postevent_features = matriz con las características para las 6 ventanas de post-event
%}
%OUTPUTS
%{
    -data_struct = ej datos.S10
    -indice = nos dice que participante es para incluirlo en la tabla
%}
% PARAMETERS
fs = 32; % Para las PPG son 32Hz
vec_events = data_struct.audioEventVector.data;

% ALGORITHM
% 1. Encontrar las posiciones de los eventos
pos = find(vec_events == 1);
pos = [pos, pos - [0; pos(1:end-1)]];
pos = pos(pos(:,2) ~= 1);

% Dejamos solo el primero de los estimulos
if(AllSignal)
    pos=(tiempo_ventana(1):tiempo_ventana(1):180)*fs;
else
    pos=pos(idx_estimulo);
end


% 2. Crear las ventanas
limit_preevent = zeros(length(pos), 2);
limit_postevent = zeros(length(pos), 2);

for i = 1:length(pos)
    % Tiempos de evento
    t_event = data_struct.audioEventVector.TimeStampDate(pos(i));

    % Ventana pre-evento
    t_pre_start = t_event - seconds(tiempo_ventana(1));
    t_pre_end = t_event;

    % Ventana post-evento
    t_post_start = t_event;
    t_post_end = t_event + seconds(tiempo_ventana(2));

    % Índices más cercanos en la señal PPG para pre-evento
    [~, idxInicio_pre] = min(abs(data_struct.caneData.(tipoExp).PPG.TimeStampDate - t_pre_start));
    [~, idxFin_pre] = min(abs(data_struct.caneData.(tipoExp).PPG.TimeStampDate - t_pre_end));
    limit_preevent(i,:) = [idxInicio_pre, idxFin_pre];

    % Índices más cercanos en la señal PPG para post-evento
    [~, idxInicio_post] = min(abs(data_struct.caneData.(tipoExp).PPG.TimeStampDate - t_post_start));
    [~, idxFin_post] = min(abs(data_struct.caneData.(tipoExp).PPG.TimeStampDate - t_post_end));
    limit_postevent(i,:) = [idxInicio_post, idxFin_post];
end

%{
% Cogemos la primera mitad sin eventos como pre, la segunda mitad como post
t_muestras=length(data_struct.caneData.(tipoExp).PPG.TimeStampDate);
limit_preevent(1,:) = [1, floor(t_muestras/3)];
limit_postevent(1,:) = [floor(t_muestras/3) , t_muestras];
%}



%% Acondicionamos las señales
% Extraemos PPG del experimento i
canePPG = data_struct.caneData.(tipoExp).PPG.data;
empaticaPPG = data_struct.empaticaData.(tipoExp).PPG.data;
% Normalizamos
canePPG = (canePPG-min(canePPG))/(max(canePPG)-min(canePPG));
empaticaPPG = (empaticaPPG-min(empaticaPPG))/(max(empaticaPPG)-min(empaticaPPG));

% Suavizamos la señal para eliminiar falsos máximos
canePPG_smooth=smooth(canePPG,1/3*fs);
empaticaPPG_smooth=smooth(empaticaPPG,1/3*fs);

% Buscamos los maximos
caneidx_max=[];
empaticaidx_max=[];
for idx = 2:length(canePPG_smooth)-1
    % Max Cane
    if canePPG_smooth(idx-1) <= canePPG_smooth(idx) && canePPG_smooth(idx) > canePPG_smooth(idx+1)
        caneidx_max=[caneidx_max,idx];
    end
end
for idx = 2:length(empaticaPPG_smooth)-1
    % Max empatica
    if empaticaPPG_smooth(idx-1) <= empaticaPPG_smooth(idx) && empaticaPPG_smooth(idx) > empaticaPPG_smooth(idx+1)
        empaticaidx_max=[empaticaidx_max,idx];
    end
end

% IBI individual
% Sacamos tiempos del experimento i
timestampCane       = data_struct.caneData.(tipoExp).PPG.TimeStampPosix;
timestampEmpatica   = data_struct.empaticaData.(tipoExp).PPG.TimeStampPosix;
% Calculamos la diferencia entre máximos
CaneIBI = (timestampCane(caneidx_max(2:end)) - timestampCane(caneidx_max(1:end-1)));
EmpaticaIBI = (timestampEmpatica(empaticaidx_max(2:end)) - timestampEmpatica(empaticaidx_max(1:end-1)));
% Normalizamos
CaneIBInorm = (CaneIBI-min(CaneIBI))/(max(CaneIBI)-min(CaneIBI));
EmpaticaIBInorm = (EmpaticaIBI-min(EmpaticaIBI))/(max(EmpaticaIBI)-min(EmpaticaIBI));
% Devolvemos tamaño original
CaneIBI=[0; CaneIBI];
EmpaticaIBI=[0; EmpaticaIBI];
CaneIBInorm=[0; CaneIBInorm];
EmpaticaIBInorm=[0; EmpaticaIBInorm];

% IBI directo de la empatica
timestamp_SYSPIBI = data_struct.empaticaData.(tipoExp).SYSTP.TimeStampPosix/1000;
value_SYSPIBI = data_struct.empaticaData.(tipoExp).SYSTP.data;
% Normalizamos
value_SYSPIBInorm = (value_SYSPIBI-min(value_SYSPIBI))/(max(value_SYSPIBI)-min(value_SYSPIBI));

%{
% Maximos de Cane
figure
plot(timestampCane,canePPG)
hold on
plot(timestampCane,canePPG_smooth)
stem(timestampCane(caneidx_max),ones(length(caneidx_max),1))
title(['Cane: Paticipante', string(indice)])
%}
%{
% Maximos de Empatica + SYSP
figure
plot(timestampEmpatica,empaticaPPG)
hold on
plot(timestampEmpatica,empaticaPPG_smooth)
stem(timestampEmpatica(empaticaidx_max),ones(length(empaticaidx_max),1))
stem(timestamp_SYSPIBI,ones(length(timestamp_SYSPIBI),1),'k')
title(['Empatica: Paticipante', string(indice)])
%}
%{
% Selección de ventanas
figure
plot(timestampEmpatica,empaticaPPG)
hold on
stem(timestampEmpatica(limit_preevent(:,1)),ones(length(limit_preevent(:,1)),1),'g')
stem(timestampEmpatica(limit_postevent(:,1)),ones(length(limit_postevent(:,1)),1),'k')
stem(timestampEmpatica(limit_postevent(:,2)),ones(length(limit_postevent(:,1)),1),'r')
%}
%{
% Graficas de PPG
figure
subplot(2,1,1)
hold on
plot(timestampCane,canePPG, 'Color', [0 0 1])
plot(convertTo(data_struct.audioEventVector.TimeStampDate,"posixtime"),data_struct.audioEventVector.data,'Color', [0 0 0])
title("PPG_n:   Participante " + string(indice))
ax = gca;
ax.XTick = [];   % Oculta las etiquetas numéricas en X
ax.YTick = [];   % Oculta las etiquetas numéricas en Y
subplot(2,1,2)
hold on
plot(timestampEmpatica,empaticaPPG, 'Color', [0.8500, 0.3250, 0.0980])
plot(convertTo(data_struct.audioEventVector.TimeStampDate,"posixtime"),data_struct.audioEventVector.data,'Color', [0 0 0])
figure
% Graficas de IBI
subplot(2,1,1)
hold on
plot(timestampCane(caneidx_max),CaneIBInorm, 'Color', [0 0 1])
plot(convertTo(data_struct.audioEventVector.TimeStampDate,"posixtime"),data_struct.audioEventVector.data,'Color', [0 0 0])
ax = gca;
ax.XTick = [];   % Oculta las etiquetas numéricas en X
ax.YTick = [];   % Oculta las etiquetas numéricas en Y
title("IBI_n:   Participante " + string(indice))
subplot(2,1,2)
hold on
plot(timestampEmpatica(empaticaidx_max),EmpaticaIBInorm, 'Color', [0.8500, 0.3250, 0.0980])
plot(convertTo(data_struct.audioEventVector.TimeStampDate,"posixtime"),data_struct.audioEventVector.data,'Color', [0 0 0])
%}

%% Iteramos cada ventana
% Ventanas de pre-evento y cálculo de características
id = strcat("S", num2str(indice));
% Inicializar tablas acumuladoras
preevent_features_cane = table();
preevent_features_empatica = table();
preevent_SYSPfeatures_empatica = table();
for i = 1:size(limit_preevent,1)

    inicio = limit_preevent(i,1);
    final = limit_preevent(i,2);

    % Cane
    % Calculo de la ventana
    windows_preevent_canePPG = canePPG(inicio:final);
    windows_preevent_caneIBI = CaneIBI(ismember(timestampCane(caneidx_max),timestampCane(inicio:final)));
    windows_preevent_caneIBInorm = CaneIBInorm(ismember(timestampCane(caneidx_max),timestampCane(inicio:final)));
    % Calculo de características
    [tablewindow_cane, ~] = calculateWindowFeature (windows_preevent_canePPG,windows_preevent_caneIBI,windows_preevent_caneIBInorm, NaN, NaN);
    tablewindow_cane.("ID participante") = repmat(string(id), height(tablewindow_cane), 1);
    preevent_features_cane = [preevent_features_cane; tablewindow_cane];
    

    %Empatica
    % Calculo de la ventana
    windows_preevent_empaticaPPG = empaticaPPG(inicio:final);
    windows_preevent_empaticaIBI = EmpaticaIBI(ismember(timestampEmpatica(empaticaidx_max),timestampEmpatica(inicio:final)));
    windows_preevent_empaticaIBInorm = EmpaticaIBInorm(ismember(timestampEmpatica(empaticaidx_max),timestampEmpatica(inicio:final)));
    windows_preevent_SYSPIBI = value_SYSPIBI((timestamp_SYSPIBI>timestampEmpatica(inicio)) & (timestamp_SYSPIBI<timestampEmpatica(final)));
    windows_preevent_SYSPIBInorm = value_SYSPIBInorm((timestamp_SYSPIBI>timestampEmpatica(inicio)) & (timestamp_SYSPIBI<timestampEmpatica(final)));
    % Calculo de características
    [tablewindow_empatica, SYSPIBItable] = calculateWindowFeature (windows_preevent_empaticaPPG,windows_preevent_empaticaIBI,windows_preevent_empaticaIBInorm, windows_preevent_SYSPIBI, windows_preevent_SYSPIBInorm);
    tablewindow_empatica.("ID participante") = repmat(string(id), height(tablewindow_empatica), 1);
    SYSPIBItable.("ID participante") = repmat(string(id), height(SYSPIBItable), 1);
    preevent_features_empatica = [preevent_features_empatica; tablewindow_empatica];
    preevent_SYSPfeatures_empatica = [preevent_SYSPfeatures_empatica; SYSPIBItable];

    %plot(timestampEmpatica(inicio:final),empaticaPPG(inicio:final),'g')
end
%Ventanas de Post-evento
% Inicializar tablas acumuladoras
postevent_features_cane = table();
postevent_features_empatica = table();
postevent_SYSPfeatures_empatica = table();

for i = 1:size(limit_postevent,1)

    inicio = limit_postevent(i,1);
    final = limit_postevent(i,2);

    %Cane
    % Calculo de la ventana
    windows_postevent_canePPG = canePPG(inicio:final);
    windows_postevent_caneIBI = CaneIBI(ismember(timestampCane(caneidx_max),timestampCane(inicio:final)));
    windows_postevent_caneIBInorm = CaneIBInorm(ismember(timestampCane(caneidx_max),timestampCane(inicio:final)));
    % Calculo de características
    tablewindow_cane = calculateWindowFeature (windows_postevent_canePPG,windows_postevent_caneIBI,windows_postevent_caneIBInorm, NaN, NaN);
    tablewindow_cane.("ID participante") = repmat(string(id), height(tablewindow_cane), 1);
    postevent_features_cane = [postevent_features_cane; tablewindow_cane];

    %Empatica
    % Calculo de la ventana
    windows_postevent_empaticaPPG = empaticaPPG(inicio:final);
    windows_postevent_empaticaIBI = EmpaticaIBI(ismember(timestampEmpatica(empaticaidx_max),timestampEmpatica(inicio:final)));
    windows_postevent_empaticaIBInorm = EmpaticaIBInorm(ismember(timestampEmpatica(empaticaidx_max),timestampEmpatica(inicio:final)));
    windows_postevent_SYSPIBI = value_SYSPIBI((timestamp_SYSPIBI>timestampEmpatica(inicio)) & (timestamp_SYSPIBI<timestampEmpatica(final)));
    windows_postevent_SYSPIBInorm = value_SYSPIBInorm((timestamp_SYSPIBI>timestampEmpatica(inicio)) & (timestamp_SYSPIBI<timestampEmpatica(final)));
   % Calculo de características
    [tablewindow_empatica, SYSPIBItable] = calculateWindowFeature (windows_postevent_empaticaPPG,windows_postevent_empaticaIBI,windows_postevent_empaticaIBInorm, windows_postevent_SYSPIBI, windows_postevent_SYSPIBInorm);
    tablewindow_empatica.("ID participante") = repmat(string(id), height(tablewindow_empatica), 1);
    SYSPIBItable.("ID participante") = repmat(string(id), height(SYSPIBItable), 1);
    postevent_features_empatica = [postevent_features_empatica; tablewindow_empatica];
    postevent_SYSPfeatures_empatica = [postevent_SYSPfeatures_empatica; SYSPIBItable];


    %plot(timestampEmpatica(inicio:final),empaticaPPG(inicio:final),'r')
end



end