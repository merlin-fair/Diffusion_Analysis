clear; clc;

%%%%%%%%%%%%%%%%%%%%%   Paths
projectDir = '/MATLAB Drive/Linux/Tractography/matriz_labels/projectoDir';

controlDir = fullfile(projectDir, 'control_connectomes/');
patientDir = fullfile(projectDir, 'schz_connectomes/');

controlSubs = dir(fullfile(controlDir,'connectome_sift2_*'));
patientSubs = dir(fullfile(patientDir, 'connectome_sift2_*'));

%%%%%%%%%%%%%%%%%%%%%   Load connectomes 
mats = [];
group = [];
%group = group(:);
subjects = {};
%subjects = subjects(:);

for i = 1:length(controlSubs)
    expectedSize = []; % will store expected NxN size of connectomes
    sub = controlSubs(i).name;
    disp(['El nombre del archivo es: ', sub]); %prueba de gemini
    f = fullfile(controlDir, sub);
    mat = readmatrix(f);
    mats(:,:,end+1) = mat;

    [r, c] = size(mat);
    if r < 84 || c < 84
        warning('Sujeto %s tiene tamaño %dx%d. Ajustando a 84x84.', sub, r, c);
        mat_fixed = zeros(84, 84);
        % Copiar la información existente a la nueva matriz (asumiendo que los nodos faltantes son los últimos)
        mat_fixed(1:r, 1:c) = mat;
        mat = mat_fixed;
    end

    [filas, columnas] = size(mat);

    if filas == 84 && columnas == 84
        fprintf('Correcto: El sujeto %s tiene una matriz de 84x84.\n', sub);
    else
        % Esto es por si el parche no logró el tamaño esperado
        error('¡ALERTA! El sujeto %s tiene un tamaño de %dx%d después de la corrección.', sub, filas, columnas);
    end

    group(end+1,1) = 0; % controls
    subjects{end+1,1} = sub;
end

for i = 1:length(patientSubs)
    sub = patientSubs(i).name;
    disp(['El nombre del archivo es: ', sub]); %preba de gemini
    f = fullfile(patientDir, sub); 

    [r, c] = size(mat);
    if r < 84 || c < 84
        warning('Sujeto %s tiene tamaño %dx%d. Ajustando a 84x84.', sub, r, c);
        mat_fixed = zeros(84, 84);
        % Copiar la información existente a la nueva matriz (asumiendo que los nodos faltantes son los últimos)
        mat_fixed(1:r, 1:c) = mat;
        mat = mat_fixed;
    end

     mats(:,:,end+1) = mat;
    group(end+1,1) = 1; % 0 for control loop, 1 for patient loop
    subjects{end+1,1} = sub;

    [filas, columnas] = size(mat);

    if filas == 84 && columnas == 84
        fprintf('Correcto: El sujeto %s tiene una matriz de 84x84.\n', sub);
    else
        % Esto es por si el parche no logró el tamaño esperado
        error('¡ALERTA! El sujeto %s tiene un tamaño de %dx%d después de la corrección.', sub, filas, columnas);
    end

end 

nSub = size(mats,3);
nNodes = size(mats,1);

fprintf('Loaded %d subjects with %d nodes.\n', nSub, nNodes);

%%%%%%%%%%%%%%%%%%%%%   Basic QC
for s = 1:nSub
    mat = mats(:,:,s);

    if size(mat,1) ~= size(mat,2)
        error('Matrix for subject %s is not square.', subjects{s});
    end

    ta = mat';
    asymmetry = max(abs(mat(:) - ta(:)));

    if asymmetry > 1e-6
        warning('Matrix for subject %s is not perfectly symmetric.', subjects{s});
    end
end

%%%%%%%%%%%%%%%%%%%%%   Normalize each subject connectome
% Recommended first-pass normalization:
% divide each matrix by total connectivity strength.

matsNorm = zeros(size(mats));

for s = 1:nSub
    mat = mats(:,:,s);

    % Remove diagonal
    mat(1:nNodes+1:end) = 0;

    totalStrength = sum(mat(:));

    if totalStrength > 0
        matsNorm(:,:,s) = mat ./ totalStrength;
    else
        matsNorm(:,:,s) = mat;
    end
end 

%%%%%%%%%%%%%%%%%%%%% Optional log transform for statistical testing
% For SIFT2 connectomes, values can be heavy-tailed.
% log1p is useful, but after total normalization values are often small.
% Use either matsNorm directly or log1p(mats).
%
% Here we use normalized matrices.

data = matsNorm;

%%%%%%%%%%%%%%%%%%%%% Subject-level global metrics
globalStrength = zeros(nSub,1);
meanConnectivity = zeros(nSub,1);
density = zeros(nSub,1);
meanNodeStrength = zeros(nSub,1);

for s = 1:nSub
    mat = data(:,:,s);

    upperIdx = triu(true(nNodes), 1);
    edges = mat(upperIdx);

    globalStrength(s) = sum(edges);
    meanConnectivity(s) = mean(edges);
    density(s) = mean(edges > 0);

    nodeStrength = sum(mat, 2);
    meanNodeStrength(s) = mean(nodeStrength);
end

n = min([numel(subjects), numel(group), numel(globalStrength), ...
    numel(meanConnectivity), numel(density), numel(meanNodeStrength)]);
subjects = subjects(1:n);
group = group(1:n);
globalStrength = globalStrength(1:n);
meanConnectivity = meanConnectivity(1:n);
density = density(1:n);
meanNodeStrength = meanNodeStrength(1:n);

metricsTable = table(subjects, group, globalStrength, meanConnectivity, density, meanNodeStrength);

writetable(metricsTable, fullfile(projectDir, 'global_connectome_metrics.csv'));


%%%%%%%%%%%%%%%%%%%%% Group comparisons: global metrics
isCtrl = group == 0;
isPatient = group == 1;

fprintf('\nGlobal metric comparisons:\n');

metricNames = {'globalStrength', 'meanConnectivity', 'density', 'meanNodeStrength'};

for m = 1:length(metricNames)
    metric = metricsTable.(metricNames{m});

    [~,p,~,stats] = ttest2(metric(isCtrl), metric(isPatient), 'Vartype', 'unequal');

    fprintf('%s: t = %.3f, p = %.5f\n', metricNames{m}, stats.tstat, p);
end

%%%%%%%%%%%%%%%%%%%%% Group mean matrices
meanCtrl = mean(data(:,:,isCtrl), 3);
meanPatient = mean(data(:,:,isPatient), 3);
diffMat = meanPatient - meanCtrl;

writematrix(meanCtrl, fullfile(projectDir, 'mean_connectome_controls.csv'));
writematrix(meanPatient, fullfile(projectDir, 'mean_connectome_schizophrenia.csv'));
writematrix(diffMat, fullfile(projectDir, 'difference_connectome_scz_minus_ctrl.csv'));

%%%%%%%%%%%%%%%%%%%%% Visualize mean matrices
figure;
imagesc(meanCtrl);
axis image;
colorbar;
title('Mean connectome: controls');

figure;
imagesc(meanPatient);
axis image;
colorbar;
title('Mean connectome: schizophrenia');

figure;
maxAbs = max(abs(diffMat(:)));
imagesc(diffMat, [-maxAbs maxAbs]);
axis image;
colorbar;
title('Difference: schizophrenia - controls');





