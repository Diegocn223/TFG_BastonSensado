function [vec_time,vec_secuencia]=playAudio(audio,secuencia,fs)

% Se reproduce
sound(audio,fs)

% Pintamos la secuencia
figure
plot((1:length(audio))/fs,secuencia,'b')


t_inicial = tic;  % Iniciar temporizador
vec_time = [];
vec_secuencia = [];
cont_reproduccion = 10000;
for cursor = 1:(fs/128):length(audio) % Generamos vectores a la frecuencia de la música
    vec_time = [vec_time,datetime('now', 'TimeZone', 'Europe/Madrid', 'Format', 'dd-MM-yyyy HH:mm:ss.SSSSSS')];
    vec_secuencia = [vec_secuencia, secuencia(floor(cursor))];
    
    cont_reproduccion = cont_reproduccion + 1;
    if(cont_reproduccion>100)
        plot((1:length(audio))/fs,secuencia,'b')
        hold on
        plot((1:floor(cursor))/fs,secuencia(1:floor(cursor)),'r-*');
        hold off
        cont_reproduccion = 0;
    end

    % Calcular el tiempo transcurrido y ajustar la pausa
    t_esperado = cursor/fs;
    t_transcurrido = toc(t_inicial);
    t_restante = t_esperado - t_transcurrido;
    
    if t_restante > 0
        pause(t_restante);
    end
end

end