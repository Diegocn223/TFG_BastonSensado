%% Proceso a seguir para cada toma de datos

% Primero de todo, indicar el número del audio a reproducir
n_audio = input('Introduce el número del audio: ');



% Carga o generación del audio

% Damos formato al nombre
path_audio = sprintf('%s%02d.mat', 'Audios\Audio', n_audio);

% Comparamos si ya existe o no
if (exist(path_audio,'file'))   
    load(path_audio);
        disp(['Audio encontrado y cargado: ',path_audio]);
else
    %[audio,secuencia,fs]=genAudio(total_time,initial_silence,silence_duration,event_duration)
     [audio,secuencia,fs]=genAudio(    180   ,       60      ,     [5,25]     ,    [5]    ); 

    save(path_audio,'audio','secuencia','fs');
        disp(['Audio creado y guardado: ',path_audio]);
end

plot((1:length(audio))/fs,secuencia,'b')

clear path_audio




% Reproducción del audio

% Seleccionamos la nueva carpeta de datos medido
carpeta_audio = sprintf('%s%02d', 'Datos\Audio', n_audio);
mkdir(carpeta_audio);
    disp(['Carpeta "', carpeta_audio,'" lista.']);
    disp(' ');
    disp(['Si desea grabar con BitBrain, introduzca esta dirección en la grabadora: (',fullfile(pwd,carpeta_audio),').']);

respuesta = "";
while ~strcmpi(respuesta, "PLAY")
    disp('Comience a grabar. Cuando esté listo, introduzca "PLAY" en este terminal de matlab:');
    respuesta = input('', 's');  % 's' indica que se espera una cadena de texto
end

% Todo este bloque para coger el mismo nombre
d = dir(carpeta_audio);
archivo_datos = fullfile(carpeta_audio, sprintf('%s%02d', '\Prueba', length(d(~[d.isdir]))+1));% Toma la última carpeta por orden de nombre

% Se reproduce
    disp(' ');
    disp('Reproduciendo...')
[vec_time,vec_secuencia]=playAudio(audio,secuencia,fs);


% Guardamos los datos de reproduccion
save([archivo_datos,'.mat'],'vec_time','vec_secuencia');
    disp('Datos de reproduccion guardados');
    disp('Fin');

clear archivo_datos carpeta_audio d nombres ultima_carpeta respuesta
