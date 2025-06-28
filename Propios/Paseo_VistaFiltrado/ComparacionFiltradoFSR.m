figure
hold on
load('C:\Users\Usuario\Escritorio\TFG\TODOS\Participantes\Toma de datos - Experimento cinta (13-06-25)\S20\raw data\11_06_2025 - 16_56_grabacion_data.mat')
%plot(caneData.TimeStampDate,caneData.PPG,'b')
load('C:\Users\Usuario\Escritorio\TFG\TODOS\Participantes\Toma de datos - Experimento cinta (13-06-25)\S20\raw data\filt_11_06_2025 - 16_56_grabacion_data.mat')
%plot(caneData.TimeStampDate,caneData.PPG,'r')
load('C:\Users\Usuario\Escritorio\TFG\TODOS\Participantes\Toma de datos - Experimento cinta (13-06-25)\S20\raw data\RAFAV2-3YK9L1J11W_06-11-25.mat')
plot(dat.BVP{1}.local_time,dat.BVP{1}.data)

plot(caneData.TimeStampDate,(caneData.FSR-mean(caneData.FSR))*10,'g')


figurePa obtener 
hold on
load('C:\Users\Usuario\Escritorio\TFG\TODOS\Participantes\Toma de datos - Experimento cinta (13-06-25)\S3\raw data\04_06_2025 - 11_33_grabacion_data.mat')
load('C:\Users\Usuario\Escritorio\TFG\TODOS\Participantes\Toma de datos - Experimento cinta (13-06-25)\S3\raw data\filt_04_06_2025 - 11_33_grabacion_data.mat')
load('C:\Users\Usuario\Escritorio\TFG\TODOS\Participantes\Toma de datos - Experimento cinta (13-06-25)\S3\raw data\RAFAV2-3YK9L1J11W_06-04-25.mat')
plot(dat.BVP{1}.local_time,dat.BVP{1}.data)


plot(caneData.TimeStampDate,(caneData.FSR-mean(caneData.FSR))*10,'g')
