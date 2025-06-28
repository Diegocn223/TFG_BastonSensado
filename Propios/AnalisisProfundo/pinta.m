function pinta(CaracAll, CaracEvent, n_participantes, tiempo_ventana,dataStruct)



%% Creamos dos copias con las medianas e iqr (All)
mediaFeatures = struct();
iqrFeatures = struct();

n_ventanas=length(CaracAll.empatica_pre.ppgMaxValue)/n_participantes;

% Iniciamos struct con caracteristicas a 0
campos = fieldnames(CaracAll);
for i = 1:numel(campos)
    campo=campos{i};
    caracs=fieldnames(CaracAll.(campo));
    mediaFeatures.(campo) = table();
    iqrFeatures.(campo) = table();
    for j = 1:(numel(caracs)-4)
        carac=caracs{j};
        mediaFeatures.(campo).(carac)=zeros(n_ventanas,1);
        iqrFeatures.(campo).(carac)=zeros(n_ventanas,2);
    end
end

% Guardamos los valores
campos = fieldnames(CaracAll);
for i = 1:numel(campos)
    campo=campos{i};
    caracs=fieldnames(CaracAll.(campo));
    for k=1:n_ventanas
        for j = 1:(numel(caracs)-4)
            carac=caracs{j};
            mediaFeatures.(campo).(carac)(k)=mean(CaracAll.(campo).(carac)(k:n_ventanas:end));
            iqrFeatures.(campo).(carac)(k,:)=[quantile(CaracAll.(campo).(carac)(k:n_ventanas:end),0.25) quantile(CaracAll.(campo).(carac)(k:n_ventanas:end),0.75)];
        end
    end
end


    
    
    
    
    %% Iteramos por cada característica
    
    caracs=fieldnames(CaracAll.empatica_pre);
    for j = 1:(numel(caracs)-4)
        carac = caracs{j};
    
    
        % Cogemos cada característica de la ventanas
        cane.original=CaracAll.cane_pre.(carac);
        cane.data=mediaFeatures.cane_pre.(carac);
        cane.iqr=iqrFeatures.cane_pre.(carac);
        empatica.original=CaracAll.empatica_pre.(carac);
        empatica.data=mediaFeatures.empatica_pre.(carac);
        empatica.iqr=iqrFeatures.empatica_pre.(carac);
        if(j>9 && j<18)
            SYSP.original=CaracAll.SYSP_pre.(carac);
            SYSP.data=mediaFeatures.SYSP_pre.(carac);
            SYSP.iqr=iqrFeatures.SYSP_pre.(carac);
            %maximo=max(max([cane.data+cane.iqr/2,empatica.data+empatica.iqr/2,SYSP.data+SYSP.iqr]));
            %minimo=min(min([cane.data-cane.iqr/2,empatica.data-empatica.iqr/2,SYSP.data-SYSP.iqr]));
            maximo=max(max([cane.original,empatica.original,SYSP.original]));
            minimo=min(min([cane.original,empatica.original,SYSP.original]));
        else
            SYSP.original=[];
            SYSP.data=[];
            SYSP.iqr=[];
            maximo=max(max([cane.original,empatica.original]));
            minimo=min(min([cane.original,empatica.original]));
        end


   
        %% Ploteamos
        f = figure;
        f.Position = [0 0 960 1080];

        %% Gráficas de arriba (All)
        % Histograma
        subplot(22,15,[1,92])
        hold on
        histogram(cane.original,20,'FaceColor', 'blue',"FaceAlpha",0.3)
        histogram(empatica.original,20,'FaceColor', [0.8500, 0.3250, 0.0980],"FaceAlpha",0.3)
        if(~isempty(SYSP.original))
            histogram(SYSP.original,20,'FaceColor', 'magenta',"FaceAlpha",0.3)
        end
        camroll(90)
        set(gca, 'XTick', [], 'XTickLabel', []);
        if minimo == maximo
            xlim([minimo - 0.1, maximo + 0.1])
            if(isnan(minimo));minimo=0;end
            if(isnan(maximo));maximo=1;end
        else
            xlim([minimo, maximo])
            if(isnan(minimo));minimo=0;end
            if(isnan(maximo));maximo=1;end
        end
    
        % Valores por ventana
        subplot(22,15,[4,105])
        hold on

        % Cane
        y=cane.data;
        x=(0:length(cane.data)-1).*tiempo_ventana(1)+tiempo_ventana(1)/2;
        plot(x, y, 'Color', 'blue')
        % Rango intercuartilico con transparencia
        x_patch = [x, fliplr(x)];
        %x_patch = min(max(x_patch, min(x)), max(x));
        y_patch = [cane.iqr(:,1); fliplr(cane.iqr(:,2))]';

        %y_patch = min(max(y_patch, min(y)), max(y));
        % Dibujamos la "línea gruesa" con transparencia
        if(n_participantes>1)
            fill(x_patch, y_patch, ...
                'b', ...
                'FaceAlpha', 0.1, ...
                'EdgeColor', 'none');
        end
        % Empatica
        y=empatica.data;
        x=(0:length(empatica.data)-1).*tiempo_ventana(1)+tiempo_ventana(1)/2;
        plot(x, y, 'Color', [0.8500, 0.3250, 0.0980])
        % Rango intercuartilico con transparencia
        x_patch = [x, fliplr(x)];
        y_patch = [empatica.iqr(:,1); fliplr(empatica.iqr(:,2))]';
        % Dibujamos la "línea gruesa" con transparencia

        if(n_participantes>1)
            fill(x_patch, y_patch, ...
                [0.8500, 0.3250, 0.0980], ...
                'FaceAlpha', 0.1, ...
                'EdgeColor', 'none');
        end
        % SYSP
        if(~isempty(SYSP.data))
            y=SYSP.data;
            x=(0:length(SYSP.data)-1).*tiempo_ventana(1)+tiempo_ventana(1)/2;
            plot(x,y,'Color', 'magenta')
            % Rango intercuartilico con transparencia
            x_patch = [x, fliplr(x)];
            y_patch = [SYSP.iqr(:,1); fliplr(SYSP.iqr(:,2))]';
            % Dibujamos la "línea gruesa" con transparencia
            if(n_participantes>1)
                fill(x_patch, y_patch, ...
                    'm', ...
                    'FaceAlpha', 0.1, ...
                    'EdgeColor', 'none');
            end
        end

        % Evento
        if(n_participantes>1)
            plot([60 60],[minimo maximo],'Color', [0.5 0.5 0.5], 'LineWidth',2)
        else
            eventos=dataStruct.S19.audioEventVector;%Habría que poner el concreto
            timestamps = seconds(eventos.TimeStampDate - eventos.TimeStampDate(1));
            scaledTimestamps = (timestamps*180/timestamps(end))-5;
            plot(scaledTimestamps,eventos.data*(maximo-minimo)+minimo,'Color', [0.5 0.5 0.5])
        end

        %{
        plot(1,1,'Color',[0 1 0 0.3],LineWidth=5)
        plot(1,1,'Color',[1 0 0 0.3],LineWidth=5)
        legend(["Mediana (Bastón)","Rango intercuartilico (Bastón)","Mediana (Pulsera)","Rango intercuartilico (Pulsera)","Mediana (SYSP)","Rango intercuartilico (SYSP)", "Ventana Pre-Evento", "Ventana Post-Evento"],'Location',"northwest")
        %}
        title(carac, 'Interpreter', 'none')
        xlim([0,180])
        if minimo == maximo
            ylim([minimo - 0.1, maximo + 0.1])
        else
            ylim([minimo, maximo])
        end
    
        
        % Q-Q plot
        % Cane
        subplot(22,15,[122,154])
        p1=qqplot(cane.original);
        p1(1).Marker='.';
        p1(1).MarkerEdgeColor =  'blue';
        p1(2).Color = 'k';
        p1(3).Color = 'k';
        xlabel('');
        ylabel('');
        xticks([]);
        yticks([]);
        title('');
        % Empatica
        subplot(22,15,[132,164])
        p2=qqplot(empatica.original);
        p2(1).Marker='.';
        p2(1).MarkerEdgeColor =  [0.8500, 0.3250, 0.0980];
        p2(2).Color = 'k';
        p2(3).Color = 'k';
        xlabel('');
        ylabel('');
        xticks([]);
        yticks([]);
        title('');
        % SYSP
        if(~isempty(SYSP.data))
            subplot(22,15,[127,159])
            p2=qqplot(SYSP.original);
            p2(1).Marker='.';
            p2(1).MarkerEdgeColor =  'magenta';
            p2(2).Color = 'k';
            p2(3).Color = 'k';
            xlabel('');
            ylabel('');
            xticks([]);
            yticks([]);
            title('');
        end


        

        %% Graficas de abajo (Event)
        % Cane
        % Independiente
        subplot(22,15,[196,245])
        hold on
        pre_hist=CaracEvent.cane_pre.(carac);
        post_hist=CaracEvent.cane_post.(carac);
        inc_hist=(min(max(pre_hist),max(post_hist))-max(min(pre_hist),min(post_hist)))/20;
        vec_pre=min(pre_hist):inc_hist:max(pre_hist)+inc_hist;
        vec_post=min(post_hist):inc_hist:max(post_hist)+inc_hist;
        if(carac=="pNN20" || carac=="pNN750");vec_pre=-5:5:105;end
        if(carac=="pNN20" || carac=="pNN750");vec_post=-5:5:105;end
        if(length(vec_pre)>2)
            histogram(pre_hist,vec_pre,'FaceColor', 'green',"FaceAlpha",0.3)
        end
        if(length(vec_post)>2)
            histogram(post_hist,vec_post,'FaceColor', 'red',"FaceAlpha",0.3)
        end
        title("Bastón")
        % Pareado
        subplot(22,15,[271,320])
        hold on
        dif_hist=CaracEvent.cane_post.(carac)-CaracEvent.cane_pre.(carac);
        inc_hist_dif=inc_hist*(max(dif_hist)-min(dif_hist))/(max(max(pre_hist),max(post_hist))-min(min(pre_hist),min(post_hist)));
        vec_dif=min(dif_hist):inc_hist_dif:(max(dif_hist)+inc_hist);
        if(length(vec_dif)<2);vec_dif=0:5:100;end
        histogram(dif_hist,vec_dif,'FaceColor', 'blue')

        % Empatica
        % Independiente
        subplot(22,15,[206,255])
        hold on
        pre_hist=CaracEvent.empatica_pre.(carac);
        post_hist=CaracEvent.empatica_post.(carac);
        inc_hist=(min(max(pre_hist),max(post_hist))-max(min(pre_hist),min(post_hist)))/20;
        vec_pre=min(pre_hist):inc_hist:max(pre_hist)+inc_hist;
        vec_post=min(post_hist):inc_hist:max(post_hist)+inc_hist;
        if(carac=="pNN20" || carac=="pNN750");vec_pre=-5:5:105;end
        if(carac=="pNN20" || carac=="pNN750");vec_post=-5:5:105;end
        if(length(vec_pre)>2)
            histogram(pre_hist,vec_pre,'FaceColor', 'green',"FaceAlpha",0.3)
        end
        if(length(vec_post)>2)
            histogram(post_hist,vec_post,'FaceColor', 'red',"FaceAlpha",0.3)
        end
        title("Pulsera")
        % Pareado
        subplot(22,15,[281,330])
        hold on
        dif_hist=CaracEvent.empatica_post.(carac)-CaracEvent.empatica_pre.(carac);
        inc_hist_dif=inc_hist*(max(dif_hist)-min(dif_hist))/(max(max(pre_hist),max(post_hist))-min(min(pre_hist),min(post_hist)));
        vec_dif=min(dif_hist):inc_hist_dif:(max(dif_hist)+inc_hist);
        if(length(vec_dif)<2);vec_dif=0:5:100;end
        histogram(dif_hist,vec_dif,'FaceColor', [0.8500, 0.3250, 0.0980])

        % SYSP
        if(~isempty(SYSP.data))
            % Independiente
            subplot(22,15,[201,250])
            hold on
            pre_hist=CaracEvent.SYSP_pre.(carac);
            post_hist=CaracEvent.SYSP_post.(carac);
            inc_hist=(min(max(pre_hist),max(post_hist))-max(min(pre_hist),min(post_hist)))/20;
            vec_pre=min(pre_hist):inc_hist:max(pre_hist)+inc_hist;
            vec_post=min(post_hist):inc_hist:max(post_hist)+inc_hist;
            if(carac=="pNN20" || carac=="pNN750");vec_pre=-5:5:105;end
            if(carac=="pNN20" || carac=="pNN750");vec_post=-5:5:105;end
            if(length(vec_pre)>2)
                histogram(pre_hist,vec_pre,'FaceColor', 'green',"FaceAlpha",0.3)
            end
            if(length(vec_post)>2)
                histogram(post_hist,vec_post,'FaceColor', 'red',"FaceAlpha",0.3)
            end
            title("SYSP")
            % Pareado
            subplot(22,15,[276,325])
            hold on
            dif_hist=CaracEvent.SYSP_post.(carac)-CaracEvent.SYSP_pre.(carac);
            inc_hist_dif=inc_hist*(max(dif_hist)-min(dif_hist))/(max(max(pre_hist),max(post_hist))-min(min(pre_hist),min(post_hist)));
            vec_dif=min(dif_hist):inc_hist_dif:(max(dif_hist)+inc_hist);
            if(length(vec_dif)<2);vec_dif=0:5:100;end
            histogram(dif_hist,vec_dif,'FaceColor', 'm')
        end

        annotation('line', [0.16 0.88], [0.48 0.48], 'LineWidth', 2, 'Color', [0.7 0.7 0.7]);

        %% Cerramos
        saveas(gcf, fullfile('.\Graficas\Graf_png', ['Graf_', carac, '.png']));
        %print(fullfile('C:\Users\Usuario\Escritorio\TFG\MIO\Memoria\Graf', ['Graf_', carac, '.pdf']),'-dpdf','-bestfit')  % Recorta márgenes automáticamente
        print(fullfile('.\Graficas\Graf_eps', ['Graf_', carac, '.eps']),'-depsc')
        close;
        
    end



end

