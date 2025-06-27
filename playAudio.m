function [vec_time,vec_secuencia]=playAudio(audio,secuencia,fs)

% Se reproduce
sound(audio,fs)

% Creamos el vector de tiempos a 32Hz
time_inic=datetime('now', 'TimeZone', 'Europe/Madrid', 'Format', 'dd-MM-yyyy HH:mm:ss.SSSSSS');
vec_time=time_inic:seconds(1/32):(time_inic+seconds(180)-seconds(1/32));

% Iniciamos el vector a rellenar con el tamaño anterior
vec_secuencia = NaN(length(vec_time));

% Pintamos la secuencia original remuestreada
figure
plot(vec_time,secuencia(floor(1:fs/32:length(secuencia))),'b')
hold on

for cursor = 1:length(vec_time) % Generamos vectores a la frecuencia de la música

    % Remuesteamos la secuancia a 32Hz 
    % (en caso de no ser diezmado estricto, nos aseguramos de coger el punto anterior --> floor)
    vec_secuencia(cursor) = secuencia(floor(cursor*(fs/32)));

    % Representamos el nuevo punto, avance del cursor a 32Hz
    plot(vec_time(cursor),vec_secuencia(cursor),'r*');

    % Relentizamos a 32Hz mediante una pausa

    t_restante = seconds(vec_time(cursor) - datetime('now', 'TimeZone', 'Europe/Madrid', 'Format', 'dd-MM-yyyy HH:mm:ss.SSSSSS'));
    if t_restante > 0
        pause(t_restante);
    end
end

end