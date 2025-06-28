%% 0: solo tablas   1: tablas y guarda las gráficas en las carpetas
Pinta=0;

%% 0: Movimiento (filtrado)   1:Parado    2:Parado (filtrado)
Parado=0;
dataStruct = loadAllData(Parado); % Parado: 0 (movimiento), 1 (parado)

recuento_par = NaN(12,6);
recuento_ind = NaN(12,6);

sumatorio_caracs = 0;
sumatorio_caracs_dif = 0;

%% Duración de la ventana [s] (Si se pone un vector de ventanas es para el mapa de calor)
inc_t = 5;
i_vec=inc_t:inc_t:60;
%i_vec=5;


%% Array de eventos a considerar (Si se pone un vector de ventanas es para el mapa de calor)
j_vec=[[1 0 0 0 0 0];[1 2 0 0 0 0];[1 2 3 0 0 0];[1 2 3 4 0 0];[1 2 3 4 5 0];[1 2 3 4 5 6]]';
%j_vec=[1 2 3 4 5]';

%% Seleccionar los participantes a eleminar (sustituir array o crear uno nuevo en la parte inferior de esta zona)
% Selección de paticipantes
todos=[2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19];
ninguno=[];
todos_19=[2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18];
mean5p10s=[2 3 6 8 10 15 16 18]; % Eliminar tantos quiere decir que no es tan decisivo
cane6e2s=[2 3 6 10 11 18];
empatica1e60s = [5 6 11 14 16]; %Se probó a eliminar varios sujetos y eliminar esta combinación fue la mejor
descarte_GSR= [3 8 10 13 14 15 16 18];

% Segun parecido entre MiAlgoritmo y SYSP
buenos=[9 6 19 15];
medios=[3 2 18 17 16 14 13 12 10];
malos=[7 5 4 11];


% 1: Bueno  2: Regular  3: Malo (pico muy desarado, los que he descartado)
cane_clas = [1 1 2 2 1 3 2 2 2 3 2 2 3 2 2 3 2 2];
empatica_clas = [2 1 2 1 1 2 2 1 2 3 2 3 1 1 2 3 2 1];

if Parado~=0
    datos_usados = union(todos(cane_clas==3 | empatica_clas==3),ninguno); 
    datos_usados =[];   % <<<----- Poner array de vectores a eliminar
else
    datos_usados = ninguno;                                                 % <<<-----Poner array de vectores a eliminar
end
participantesRechazados=datos_usados;

%% Funcion
contador=0;
for i = i_vec
tiempo_ventana=[i i];
for j = j_vec
    if(i>(120-sum(j~=0)*(15+5)))
        continue %Recorte de la ventana
    end
idx_estimulo=[nonzeros(j)];
%AllSignal=0; %1:toda la musica, 0:solo los eventos
Exp=["Relax1", "Music", "Relax2", "Arithmetic"];
tipoExp=Exp(2);% No está implementado el resto de señales

% Extraemos las características solo de los eventos (AllSignal = 0)
[CaracEvent] = calculateFeaturesAllwindows(dataStruct, tiempo_ventana, idx_estimulo, participantesRechazados,0,tipoExp);

% Pintamos
if(Pinta)
    % Extraemos las características de todo el audio (AllSignal = 1)
    [CaracAll] = calculateFeaturesAllwindows(dataStruct, tiempo_ventana, idx_estimulo, participantesRechazados,1,tipoExp);
    if(Parado)
        pinta(CaracAll,CaracEvent,18-length(participantesRechazados),tiempo_ventana,dataStruct)
    else
        pinta(CaracAll,CaracEvent,10,tiempo_ventana,dataStruct)
    end
end

[tablasNormalidad] = estudionormalidad(CaracEvent);

[tablasSignificancia]=AnalisisSignificancia(CaracEvent, tablasNormalidad);

recuento_ind(i/inc_t,sum(j~=0)) = sum([tablasSignificancia.cane.Significancia;tablasSignificancia.SYSP.Significancia;tablasSignificancia.empatica.Significancia]);
recuento_par(i/inc_t,sum(j~=0)) = sum([tablasSignificancia.cane.Significancia_Dif;tablasSignificancia.SYSP.Significancia_Dif;tablasSignificancia.empatica.Significancia_Dif]);

sumatorio_caracs = sumatorio_caracs + sum([tablasSignificancia.cane.Significancia,tablasSignificancia.empatica.Significancia,logical([0;0;0;0;0;0;0;0;0;tablasSignificancia.SYSP.Significancia;0;0;0])]');
sumatorio_caracs_dif = sumatorio_caracs_dif + sum([tablasSignificancia.cane.Significancia_Dif,tablasSignificancia.empatica.Significancia_Dif,logical([0;0;0;0;0;0;0;0;0;tablasSignificancia.SYSP.Significancia_Dif;0;0;0])]');

contador=contador+1;
 

end
end

porcentaje_caracs=ceil(sumatorio_caracs./([2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 2 2 2]*contador)*100);
porcentaje_caracs_dif=ceil(sumatorio_caracs_dif./([2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 2 2 2]*contador)*100);

clearvars -except tablasNormalidad tablasSignificancia CaracAll CaracEvent recuento_ind recuento_par porcentaje_caracs porcentaje_caracs_dif