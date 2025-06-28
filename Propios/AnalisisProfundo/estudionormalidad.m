function [tablasNormalidad] = estudionormalidad(tableStruct)
% Script para test de normalidad Shapiro-Wilk y Q-Q plots
% Variables necesarias: tablecanepreevent, tablecanepostevent,
% tableempaticapreevent, tableempaticapostevent

% Nivel de significación
alpha = 0.05;

% Definición de los grupos
grupos = {
    struct('pre', tableStruct.cane_pre, 'post', tableStruct.cane_post, 'nombre', 'cane'), ...
    struct('pre', tableStruct.empatica_pre, 'post', tableStruct.empatica_post, 'nombre', 'empatica'),...
    struct('pre', tableStruct.SYSP_pre, 'post', tableStruct.SYSP_post, 'nombre', 'SYSP')
};

tablasNormalidad = struct("cane",{NaN}, "SYSP", {NaN},"empatica",{NaN});

for g = 1:length(grupos)
    grupo = grupos{g};
    pre = grupo.pre;
    post = grupo.post;
    nombre = grupo.nombre;

    nombresVar = pre.Properties.VariableNames;
    numVars = numel(nombresVar);

    % Excluir la última variable (se asume que no es característica)
    numVarsAnalizar = numVars - 1;

    % Crear carpeta para guardar los Q-Q plots
    %{
    carpetaQQ = ['qqplots_', nombre];
    if ~exist(carpetaQQ, 'dir')
        mkdir(carpetaQQ);
    end
    %}

    % Inicializar resultados
    caracteristicas = nombresVar(1:numVarsAnalizar)';
    pValorPre = zeros(numVarsAnalizar, 1);
    pValorPost = zeros(numVarsAnalizar, 1);
    normalidadPre = strings(numVarsAnalizar, 1);
    normalidadPost = strings(numVarsAnalizar, 1);
    normalidadAmbos = strings(numVarsAnalizar, 1);

    for i = 1:numVarsAnalizar
        x = pre.(nombresVar{i});
        y = post.(nombresVar{i});
%{
        % Q-Q plot Pre
        figure;
        qqplot(x);
        title(['Q-Q Plot Pre: ', nombresVar{i}]);
        saveas(gcf, fullfile(carpetaQQ, ['qqplot_pre_', nombresVar{i}, '.png']));
        close;

        % Q-Q plot Post
        figure;
        qqplot(y);
        title(['Q-Q Plot Post: ', nombresVar{i}]);
        saveas(gcf, fullfile(carpetaQQ, ['qqplot_post_', nombresVar{i}, '.png']));
        close;
%}
        % Eliminar NaNs
        x = x(~isnan(x));
        y = y(~isnan(y));

        % Validar que hay suficientes datos
        if (numel(x) >= 3) && (max(x)~=min(x))
            [h_pre, p_pre] = swtest(x, alpha);
        else
            p_pre = NaN; h_pre = 1; % No normal si hay pocos datos o constante
        end

        if (numel(y) >= 3)  && (max(y)~=min(y))
            [h_post, p_post] = swtest(y, alpha);
        else
            p_post = NaN; h_post = 1;
        end

        pValorPre(i) = p_pre;
        pValorPost(i) = p_post;
        % "Normal" : "No normal";
        if(h_pre == 0) normalidadPre(i) = "normal"; else normalidadPre(i) ="no normal"; end
        if(h_post == 0) normalidadPost(i) = "normal"; else normalidadPost(i) ="no normal"; end

        if h_pre == 0 && h_post == 0
            normalidadAmbos(i) = 'Normales';
        else
            normalidadAmbos(i) = 'No normales';
        end
    end

    % Crear tabla de resultados
    tablaNormalidad = table(caracteristicas, pValorPre, normalidadPre, pValorPost, normalidadPost, normalidadAmbos, ...
        'VariableNames', {'Caracteristica', 'pValorPre', 'NormalidadPre', 'pValorPost', 'NormalidadPost', 'NormalidadAmbos'});

    % Mostrar en consola
    disp(['Resultados del test de normalidad Shapiro-Wilk para: ', nombre]);
    disp(tablaNormalidad);
%{
    % Guardar como CSV
    archivo = ['normalidad_', nombre, '.csv'];
    writetable(tablaNormalidad, archivo);
end
%}
    tablasNormalidad.(nombre) = tablaNormalidad;

disp('✅ Análisis de normalidad completado.');
end
