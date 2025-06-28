function [audio,secuencia,fs]=genAudio(total_time, initial_silence,silence_duration,event_duration)

% Inicializacion

% Cargamos las carpetas con los audios
dir_bueno='Audios\ListaBueno\';
dir_malo='Audios\ListaMalo\';
lista_bueno=dir([dir_bueno,'*mp3']);
lista_malo=dir([dir_malo,'*mp3']);
% Y reordenamos aleatoriamente la lista mala
lista_malo = lista_malo(randperm(length(lista_malo)));
iter_malo = 1;

% Sacamos uno calmado aleatoriamente, será la base sobre la que suene la
% fuente de estrés
[bueno, fs_bueno] = audioread([dir_bueno,lista_bueno(randi(length(lista_bueno))).name]);
fs=fs_bueno;
% Si no es estereo, lo duplicamos
if(width(bueno)==1)
    bueno = cat(2,bueno,bueno);
end

% Parametros: (los tachados no estan implementados, pero se podria)
t_total    = total_time;              % Tiempo total
t_inic      = initial_silence;               % Segundos iniciales de relajación
%n          = 7;                % Numero de muestras estresantes
t_min_sil   = silence_duration(1);               % Segundos minimos después de la muestra
t_max_sil   = silence_duration(min(2, end));     % Segundos maximos después de la muestra
t_min_act   = event_duration(1);                 % Segundos minimos que dura la muestra
t_max_act   = event_duration(min(2, end));       % Segundos maximos que dura la muestra 


% Recortamos el audio al total (3min)
bueno= bueno(1:t_total*fs,:);

% Iniciamos el audio malo con el silencio inicial
audio_malo = zeros(t_inic*fs_bueno,2);
% Marcamos con ceros este silencio en el vector guia
secuencia = zeros(t_inic*fs_bueno,1);


% Bucle de creacion del audio malo
% Hasta que tengan la misma longitud
while(length(audio_malo) < length(bueno))

    % Cargamos un audio malo por orden de la lista aleatoria
    [malo, fs_malo] = audioread([dir_malo,lista_malo(iter_malo).name]);
    iter_malo = iter_malo + 1;
    % Si el audio malo tiene una frecuencia distinta, lo remuestreamos
    if fs_malo ~= fs_bueno
        malo = resample(malo, fs_bueno, fs_malo);
    end
    % Si no es estereo, lo duplicamos
    if(width(malo)==1)
        malo = cat(2,malo,malo);
    end

    % Incluimos el audio
    t_activo = randi([t_min_act,t_max_act]);
    long_malo = min(t_activo*fs_bueno,length(malo));
    audio_malo = cat(1,audio_malo,malo(1:long_malo,:));
    secuencia = cat(1,secuencia,ones(long_malo,1));

    % Incluimos un silencio previo (entre muestras)
    t_silencio = randi([t_min_sil,t_max_sil]);
    audio_malo = cat(1,audio_malo,zeros(t_silencio*fs_bueno,2));
    secuencia = cat(1,secuencia,zeros(t_silencio*fs_bueno,1));
    
end

% Recortamos el exceso para que queden de la misma longitud
% (normalmente esto genera que al final quede relajado)
audio_malo = audio_malo(1:length(bueno),:);
secuencia = secuencia(1:length(bueno));

% Juntamos los dos audios y representamos
audio=bueno+audio_malo;
%audio=mean(bueno,audio_malo);

% Representamos los sustos
plot((1:length(bueno))/fs_bueno,secuencia,'b')
ylabel('Estimulación');
title('Secuencia de Estímulos');
% Calcular los ticks del eje X en minutos y segundos
duracion_total = length(bueno) / fs_bueno;
xt = 0:30:duracion_total;  % Crear ticks cada 30 segundos
xticklabels_str = arrayfun(@(x) sprintf('%02d:%02d', floor(x/60), mod(x,60)), xt, 'UniformOutput', false);
% Aplicar las etiquetas personalizadas al eje X
xticks(xt);
xticklabels(xticklabels_str);
xtickangle(-45);

yticks([0,1]);

end