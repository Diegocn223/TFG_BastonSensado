function [tablageneral, tablaIBI] = calculateWindowFeature (windows_PPG,windows_IBI,windows_IBInorm, windows_SYSPIBI, windows_SYSPIBInorm) 

tabla={table(), table()};
   
for i = 1:2

    if i==1
        IBI = windows_IBI;
        IBInorm = windows_IBInorm;
    else
        windows_PPG=NaN;
        IBI = windows_SYSPIBI;
        IBInorm = windows_SYSPIBInorm;
    end


    % Magnitude Caracteristics
    % Feature 1: Max PPG value
    ppgMaxValue = max(windows_PPG);
    % Feature 2: Min  PPG value
    ppgMinValue = min(windows_PPG);
    % Feature 3: Mean PPG value
    ppgMean = mean(windows_PPG);
    % Feature 4: Median PPG value
    ppgMedian = median(windows_PPG);
    % Feature 5: Standar Deviation PPG value
    ppgSD = std(windows_PPG);
    % Feature 6: InterQuartilic Range PPG value
    ppgIQR = iqr(windows_PPG);
    % Feature 7: Mean Absolute Deviation PPG value
    ppgMAD_Mean = mad(windows_PPG);
    % Feature 8: Median Absolute Deviation PPG value
    ppgMAD_Median = mad(windows_PPG,1);
    % Feature 9: absolute Mean Difference PPG value
    ppgMD = mean(mean(abs(windows_PPG(:) - windows_PPG(:)')));


    % Time Caracteristics
    % Feature 10: Mean IBI
    MeanIBI = mean(IBInorm);
    % Feature 11: Heart Rate
    % (sin normalizar para obtener valores reales)
    HR = 60/mean(IBI);
    % Feature 12: Standar Deviation of IBI
    SDNN = sqrt(sum((IBInorm-MeanIBI).^2)/(length(IBInorm)-1));
    % Feature 13: Root Mean Square of Sucesives IBI interval diference
    RMSSD = sqrt(sum((IBInorm(2:end)-IBInorm(1:end-1)).^2)/(length(IBInorm)-1));
    % Feature 14: Sucesive Differences Standardar Deviation
    SucDif_Cane=abs(IBInorm(2:end)-IBInorm(1:end-1));
    SDSD = sqrt(sum((SucDif_Cane-mean(SucDif_Cane)).^2)/(length(IBInorm)-2));
    % Feature 15_1: Percent of IBI larger than 20ms
    % (sin normalizar para obtener valores reales)
    pNN20 = mean((IBI(2:end)-IBI(1:end-1))>0.020)*100;
    % Feature 15_2: Percent of IBI larger than 20ms
    % (sin normalizar para obtener valores reales)
    pNN50 = mean((IBI(2:end)-IBI(1:end-1))>0.050)*100;
    % Feature 16: Percent of IBI larger than 50ms
    % (sin normalizar para obtener valores reales)
    pIBI750 = mean(IBI>0.750)*100;
    % Feature 17: Baevsky's Stress Index
    % Calculamos el hsitograma primero: incremeto de 50ms
    t_intervalos=0.05;
    if(length(IBInorm)>2)
        intervalos_hist=[(min(IBInorm):t_intervalos:max(IBInorm)), max(IBInorm)+t_intervalos];% Ponemos un punto extra para captar el final
        if(length(intervalos_hist)<2)
            % Si es muy pequeño, lo dividimos en dos regiones
            incremento_t =(max(IBInorm)-min(IBInorm))/4;
            intervalos_hist = [min(IBInorm)+incremento_t max(IBInorm)-incremento_t]; % Puntos centralos a los rangos
        end
        histograma = histcounts(IBInorm, intervalos_hist);
        [caneAMo_pre, canehist_idx] = max(histograma);
        caneMo = intervalos_hist(canehist_idx);
        if caneMo ==0
            caneMo = t_intervalos;
        end
        caneAMo = caneAMo_pre/length(IBInorm);
        caneMxDMn = max(IBInorm) - min(IBInorm);
        SI = (caneAMo/(2*caneMo*caneMxDMn))*100; % Para el bastón 
        if(SI==inf)
            SI = NaN;
        end
    else
        SI = NaN;
    end

    % Frecuency Caracteristics
    fs=32; % Frecuencia de la PPG
    VL_cutoff = 0.04 / (fs / 2); % VL (0-0.04Hz)
    L_cutoff = 0.15 / (fs / 2);  % L (0.04-0.15Hz)
    H_cutoff = 0.4 / (fs / 2);  % H (0.15-0.4Hz)
    % Diseño del filtro Butterworth de orden 4
    [b_low, a_low] = butter(4, [VL_cutoff, L_cutoff], 'bandpass');
    [b_high, a_high] = butter(4, [L_cutoff, H_cutoff], 'bandpass');
    % Feature 18: Low Frecuency Power
    LF = mean(abs(filter(b_low, a_low, windows_PPG)).^2) / fs;
    % Feature 19: High Frecuency Power
    HF = mean(abs(filter(b_high, a_high, windows_PPG)).^2) / fs;
    % Feature 20: Low Frecuency to High Frecuency Power Rate
    LFHF = LF/HF;
    
    %% Construcción de la tabla
    if(i==1)
        tabla{i} = table(...
            ppgMaxValue, ...
            ppgMinValue, ...
            ppgMean, ...
            ppgMedian, ...
            ppgSD, ...
            ppgIQR, ...
            ppgMAD_Mean, ...
            ppgMAD_Median, ...
            ppgMD, ...
            ...
            MeanIBI, ...
            HR, ...
            SDNN, ...
            RMSSD, ...
            SDSD, ...
            pNN20, ...
            pNN50, ...
            pIBI750, ...
            SI, ...
            ...
            LF, ...
            HF, ...
            LFHF...
            );
    else
        tabla{i} = table(...
            ...
            MeanIBI, ...
            HR, ...
            SDNN, ...
            RMSSD, ...
            SDSD, ...
            pNN20, ...
            pNN50, ...
            pIBI750, ...
            SI ...
            ...
            );
    end

end
tablageneral=tabla{1};
tablaIBI=tabla{2};

end