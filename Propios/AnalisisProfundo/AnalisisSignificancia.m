function [tablasSignificancia]=AnalisisSignificancia(tablasCarac, tablasNormalidad)

tablasSignificancia = struct("cane", [],"empatica", [],"SYSP", []);

% Cane y Empatica
for ins = fieldnames(tablasNormalidad)'
    inst = ins{1};

    caracteristicas = tablasNormalidad.(inst).Caracteristica;
    
    h=NaN(length(caracteristicas),1);
    pvalor=NaN(length(caracteristicas),1);
    
    % Caracteristicas
    for i = 1:length(tablasNormalidad.(inst).Caracteristica)
        %% Dos muestras
        carac = caracteristicas{i};
        if(tablasNormalidad.(inst).NormalidadAmbos(i)=="Normales")
            % Ambos Normales
            [h(i),pvalor(i)] = ttest2(tablasCarac.([inst,'_pre']).(carac),tablasCarac.([inst,'_post']).(carac));
        else
            % No Normales
            [pvalor(i),h(i)] = ranksum(tablasCarac.([inst,'_pre']).(carac),tablasCarac.([inst,'_post']).(carac));
        end

        %% Pareado
        dif_carac = tablasCarac.([inst,'_post']).(carac) - tablasCarac.([inst,'_pre']).(carac);
        % Validar que hay suficientes datos
        if (numel(dif_carac) >= 3) && (max(dif_carac)~=min(dif_carac))
            [hnorm_dif_bin(i), pnorm_dif(i)] = swtest(dif_carac);

            if(hnorm_dif_bin(i)==0)
                hnorm_dif(i)="Normal";
            else
                hnorm_dif(i)="No normal";
            end
        else
            pnorm_dif(i) = NaN; hnorm_dif(i) = "No normal"; % No normal si hay pocos datos o constante
        end

        if(hnorm_dif(i)=="Normal")
            % Diferencia Normal
            [h_dif(i),pvalor_dif(i)] = ttest(tablasCarac.([inst,'_pre']).(carac),tablasCarac.([inst,'_post']).(carac));
            %[h_dif2(i),pvalor_dif2(i)] = ttest(dif_carac);
        else
            % Diferencia No Normal
            [pvalor_dif(i),h_dif(i)] = signrank(tablasCarac.([inst,'_pre']).(carac),tablasCarac.([inst,'_post']).(carac));
            %[pvalor_dif2(i),h_dif2(i)] = signrank(dif_carac);
        end



    end


    tablasSignificancia.(inst) = table(caracteristicas, tablasNormalidad.(inst).NormalidadAmbos, ...
        pvalor, logical(h), hnorm_dif', pvalor_dif', logical(h_dif'),...
    'VariableNames',{'Caracteristica','Normalidad','pValor','Significancia','Normalidad_Dif','pValor_Dif','Significancia_Dif'});
    
    disp(['Tabla de significancia de ', (inst)])
    disp(tablasSignificancia.(inst))

    clear pvalor h hnorm_dif pvalor_dif h_dif
end

end