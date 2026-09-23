clear; 
clearvars; 
clc;        
close all;  


%% ========================================================================
% ORIGINAL / MODIFIED: Paths and analysis settings
% =========================================================================
projectDir = %Path name

controlDir = fullfile(projectDir, %Floder's name);
patientDir = fullfile(projectDir,  %Floder's name);

participantsFile = fullfile(projectDir, 'participants.tsv');

mrtrixFsDefaultFile = ''; 
minEdgePresence = 0.50; 
fdrAlpha = 0.05; 

outputDir = fullfile(projectDir, 'connectome_group_analysis');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

controlSubs = dir(fullfile(controlDir, 'sub-*'));
patientSubs = dir(fullfile(patientDir, 'sub-*'));

if isempty(controlSubs) || isempty(patientSubs)
    error('No control and/or schizophrenia subject folders were found.');
end

clear mats group subjects
mats = [];
group = [];
subjects = {};
expectedNNodes = [];

Totalsubjects = [controlSubs; patientSubs]
%whos("Totalsubjects")
nControlSubjects = length(controlSubs);
mats = zeros(84, 84, numel(Totalsubjects));
k = 0;

% ---------------- Controls ----------------
for i = 1:length(controlSubs)
    k = k + 1;
    sub = controlSubs(i).name;
    subIDcontrol = erase(sub, '.csv');
    subjects{end+1,1} = subIDcontrol;
    group(end+1,1) = 0; 
    %disp(['El nombre del archivo es: ', sub]); 
    f = fullfile(controlDir, sub) 
    
    if ~isfile(f)
        error('Missing connectome for subject %s: %s', sub, f);
    end
    
    mat = readmatrix(f);
    [r, c] = size(mat);
    if r < 84 || c < 84 
        warning('Sujeto %s tiene tamaño %dx%d. Ajustando a 84x84.', sub, r, c);
        mat_fixed = zeros(84, 84);
        mat_fixed(1:r, 1:c) = mat;
        mat = mat_fixed; 
    end

    mats(:,:,k) = mat;  

    fprintf('i = %d | sujeto = %s | matrices en mats = %d\n', ...
        k, subIDcontrol, k);
    [filas, columnas] = size(mat); 

    if filas == 84 && columnas == 84
        fprintf('Correcto: El sujeto %s tiene una matriz de 84x84.\n', sub);
    else
        %Esto es por si el parche no logró el tamaño esperado
        error('¡ALERTA! El sujeto %s tiene un tamaño de %dx%d después de la corrección.', sub, filas, columnas);
    end 

    if size(mat,1) ~= size(mat,2)
        error('Connectome for subject %s is not square.', sub);
    end

    if any(~isfinite(mat(:)))
        error('Connectome for subject %s contains NaN or Inf values.', sub);
    end

    if isempty(expectedNNodes)
        expectedNNodes = size(mat,1);
    elseif size(mat,1) ~= expectedNNodes
        error('Connectome for subject %s has %d nodes; expected %d.', ...
            sub, size(mat,1), expectedNNodes);
    end
end

% ---------------- Schizophrenia ----------------
for i = 1:length(patientSubs)
    k = k + 1;
    sub = patientSubs(i).name;
    subIDpatient = erase(sub, '.csv');
     %disp(['El nombre del archivo es: ', sub]);
     subjects{end+1,1} = subIDpatient;
     group(end+1,1) = 1; % patients

    f = fullfile(patientDir, sub);

    if ~isfile(f)
        error('Missing connectome for subject %s: %s', sub, f);
    end

    mat = readmatrix(f);
    [r, c] = size(mat);
    if r < 84 || c < 84
        warning('Sujeto %s tiene tamaño %dx%d. Ajustando a 84x84.', sub, r, c);
        mat_fixed = zeros(84, 84);
        % Copiar la información existente a la nueva matriz (asumiendo que los nodos faltantes son los últimos)
        mat_fixed(1:r, 1:c) = mat;
        mat = mat_fixed;
    end
 
    mats(:,:,k) = mat; 
    fprintf('i = %d | sujeto = %s | matrices en mats = %d\n', ...
        k, subIDpatient, k);
   
    [filas, columnas] = size(mat); 


    if filas == 84 && columnas == 84
        fprintf('Correcto: El sujeto %s tiene una matriz de 84x84.\n', sub);
    else
        Esto es por si el parche no logró el tamaño esperado
        error('¡ALERTA! El sujeto %s tiene un tamaño de %dx%d después de la corrección.', sub, filas, columnas);
    end

    
    if size(mat,1) ~= size(mat,2)
        error('Connectome for subject %s is not square.', sub);
    end

    if any(~isfinite(mat(:)))
        error('Connectome for subject %s contains NaN or Inf values.', sub);
    end

    if isempty(expectedNNodes)
        expectedNNodes = size(mat,1);
    elseif size(mat,1) ~= expectedNNodes
        error('Connectome for subject %s has %d nodes; expected %d.', ...
            sub, size(mat,1), expectedNNodes);
    end
    
    if size(mat,1) ~= expectedNNodes
        error('Connectome for subject %s has %d nodes; expected %d.', ...
            sub, size(mat,1), expectedNNodes);
    end
    
end

nSub = size(mats,3);
nNodes = size(mats,1);

fprintf('Loaded %d subjects with %d nodes.\n', nSub, nNodes);
fprintf('  Controls:      %d\n', sum(group == 0));
fprintf('  Schizophrenia: %d\n', sum(group == 1));

%.............................................................................................

% NEW 2: Read participants.tsv and match ONLY the included subjects
% =========================================================================
% WHY:
% The full participants.tsv contains subjects who were not included
% in this study. Therefore demographics are matched AFTER the
% connectomes are loaded, using the exact subject folder names in 'subjects'.
%
% Age and gender are included as PRESPECIFIED covariates in the primary
% edge-wise and node-wise linear models. Their inclusion does not depend on
% whether the preliminary demographic comparison happens to be significant.

if ~isfile(participantsFile)
    error('participants.tsv not found: %s', participantsFile);
end

% Borrar el sufijo de los participantes: 

participants = readtable(participantsFile, ...
    'FileType', 'text', ...
    'Delimiter', '\t', ...
    'TextType', 'string', ...
    'TreatAsMissing', {'n/a','NA','NaN'});

requiredVars = {'participant_id','diagnosis','age','gender'};
for v = 1:numel(requiredVars)
    if ~ismember(requiredVars{v}, participants.Properties.VariableNames)
        error('participants.tsv must contain a column named "%s".', requiredVars{v});
    end
end
subjectStrings = string(subjects);
participantIDs = string(participants.participant_id);


% Match the participants table to the connectomes included here.
[found, demoIdx] = ismember(subjectStrings, participantIDs);
if any(~found)
    missingSubjects = strjoin(subjectStrings(~found), ', ');
    error('The following included subjects are missing from participants.tsv: %s', ...
        missingSubjects);
end

% Reorder participant information to exactly match connectome order.
if isnumeric(participants.age)
    age = double(participants.age(demoIdx));
else
    age = str2double(string(participants.age(demoIdx)));
end


genderString = upper(strtrim(string(participants.gender(demoIdx))));
diagnosisString = upper(strtrim(string(participants.diagnosis(demoIdx))));

if any(~isfinite(age))
    error('Missing or invalid age values were found for included subjects.');
end
if any(ismissing(genderString) | strlength(genderString) == 0)
    error('Missing gender values were found for included subjects.');
end
if any(ismissing(diagnosisString) | strlength(diagnosisString) == 0)
    error('Missing diagnosis values were found for included subjects.');
end

% NEW: Verify that the diagnosis recorded in participants.tsv agrees with
% the directory from which each connectome was loaded.
expectedDiagnosis = strings(nSub,1);
expectedDiagnosis(group == 0) = "CONTROL";
expectedDiagnosis(group == 1) = "SCHZ";


diagnosisString = string(diagnosisString);
expectedDiagnosis = string(expectedDiagnosis);

size(diagnosisString)
size(expectedDiagnosis)

mismatch = diagnosisString ~= expectedDiagnosis;
if any(mismatch)
    mismatchText = strjoin(subjectStrings(mismatch) + " (TSV=" + ...
        diagnosisString(mismatch) + ", folder=" + expectedDiagnosis(mismatch) + ")", ', ');
    error('Diagnosis mismatch between participants.tsv and group folders: %s', mismatchText);
end 

gender = categorical(genderString);
groupCat = categorical(group, [0 1], {'Control','Schizophrenia'});

isCtrl = (group == 0);
isPatient = (group == 1);

%.............................................................................................
% ----- Age comparison: Welch two-sample t-test -----
[~, pAge, ~, ageStats] = ttest2(age(isCtrl), age(isPatient), ...
    'Vartype', 'unequal');

fprintf('\nDEMOGRAPHIC CHECKS: ACTUAL INCLUDED CONNECTOME SAMPLE\n');
fprintf('----------------------------------------------------\n');
fprintf('Age: controls %.2f +/- %.2f years; schizophrenia %.2f +/- %.2f years\n', ...
    mean(age(isCtrl)), std(age(isCtrl)), ...
    mean(age(isPatient)), std(age(isPatient)));
fprintf('Age group comparison: Welch t = %.3f, p = %.5f\n', ...
    ageStats.tstat, pAge);
    
%.............................................................................................
% ----- Gender comparison -----
% Chi-square is used generally. For a 2x2 table with any expected count <5,
% Fisher's exact test is used instead.
[genderCounts, chi2Gender, pGenderChi] = crosstab(groupCat, gender);
expectedGender = sum(genderCounts,2) * sum(genderCounts,1) / sum(genderCounts(:));

if isequal(size(genderCounts), [2 2]) && any(expectedGender(:) < 5)
    [~, pGender] = fishertest(genderCounts);
    genderTestName = 'Fisher exact';
    genderStatistic = NaN;
else
    pGender = pGenderChi;
    genderTestName = 'Chi-square';
    genderStatistic = chi2Gender;
end

fprintf('Gender contingency table (rows = group; columns = gender categories):\n');
disp(genderCounts);
if strcmp(genderTestName, 'Chi-square')
    fprintf('Gender group comparison: chi-square = %.3f, p = %.5f\n', ...
        genderStatistic, pGender);
else
    fprintf('Gender group comparison: Fisher exact p = %.5f\n', pGender);
end
fprintf('Gender categories (column order): %s\n', ...
    strjoin(string(categories(gender)), ', '));

% Save exactly the demographic rows used in this analysis.
alignedDemo = table(subjectStrings, diagnosisString, groupCat, age, gender, ...
    'VariableNames', {'Subject','DiagnosisFromTSV','Group','Age','Gender'});
writetable(alignedDemo, fullfile(outputDir, 'participants_used_in_connectome_analysis.csv'));

%.............................................................................................
%% ========================================================================
% MODIFIED 3: QC, zero diagonal, and enforce exact symmetry
% =========================================================================
% WHY:
% Structural connectomes should be undirected/symmetric in this analysis.
% Small numerical asymmetries are reported, then the matrix is explicitly
% symmetrized. The diagonal is removed BEFORE any subsequent calculation.


matsRaw = zeros(size(mats));

for s = 1:nSub
    mat = mats(:,:,s);
    matT = mat';
  
    asymmetry = max(abs(mat(:) - matT(:)));
    if asymmetry > 1e-6
        warning('Matrix for subject %s is not perfectly symmetric (max = %.3g).', ...
            subjects{s}, asymmetry);
    end

    % Enforce exact symmetry.
    mat = (mat + mat') ./ 2;

    % Remove self-connections.
    mat(1:nNodes+1:end) = 0;

    % SIFT2 weights should not be negative.
    if any(mat(:) < 0)
        error('Negative connectivity weights found for subject %s.', subjects{s});
    end

    matsRaw(:,:,s) = mat;
end
%.............................................................................................

%% ========================================================================
% NEW / MODIFIED 4: Keep three distinct connectome representations
% =========================================================================
% WHY:
% The original script divided every subject by total connectivity strength
% and then used that normalized matrix for all analyses. This forces each
% subject's connectome to have the same total weight, so positive/negative
% differences become RELATIVE redistributions of connectivity.
%
% Here we retain:
%   matsRaw  = raw SIFT2 edge weights, diagonal removed
%   matsLog  = log1p(raw SIFT2), used for PRIMARY edge-wise statistics
%   matsNorm = raw SIFT2 / total matrix strength, used only as a SECONDARY
%              sensitivity analysis of relative connectivity.

matsLog = log1p(matsRaw);
matsNorm = zeros(size(matsRaw));

for s = 1:nSub
    mat = matsRaw(:,:,s);
    totalStrengthBothTriangles = sum(mat(:));

    if totalStrengthBothTriangles > 0
        matsNorm(:,:,s) = mat ./ totalStrengthBothTriangles;
    else
        matsNorm(:,:,s) = mat;
    end
end

%.............................................................................................
%% ========================================================================
% NEW 5: Read node labels directly from MRtrix fs_default.txt
% =========================================================================
% WHY:
% These connectomes use the FreeSurfer Desikan-Killiany parcellation after
% MRtrix labelconvert. Therefore the most reliable source of row/column names
% is the SAME MRtrix lookup table that defines the connectome node ordering:
%   share/mrtrix3/labelconvert/fs_default.txt
%
% This avoids maintaining a separate node_labels.csv and reduces the risk of
% accidentally applying labels from a different atlas or node ordering.

fsDefaultFile = %Path with the labels of Freesurfer. 
fprintf('\nCONNECTOME NODE LABELS\n');
fprintf('-----------------------\n');
fprintf('Using MRtrix lookup table:\n%s\n', fsDefaultFile);

fid = fopen(fsDefaultFile, 'r');
if fid < 0
    error('Could not open MRtrix lookup table: %s', fsDefaultFile);
end

lut = textscan(fid, '%f %s %s %f %f %f %f', ...
    'CommentStyle', '#', ...
    'MultipleDelimsAsOne', true);
fclose(fid);

nodeIndex = lut{1};
nodeLabels = string(lut{2});
nodeFreeSurferNames = string(lut{3});

% Se compara la cantidad de regiones leídas en el .txt contra el número de nodos de tus matrices
if numel(nodeLabels) ~= nNodes
    error(['MRtrix fs_default.txt contains %d nodes, but the connectomes contain %d. ' ...
        'This strongly suggests that the connectomes were generated with a different ' ...
        'parcellation / lookup table. Do not force these labels onto the matrices.'], ...
        numel(nodeLabels), nNodes);
end

if ~isequal(nodeIndex(:), (1:nNodes)')
    error(['Node indices in fs_default.txt are not exactly 1:nNodes. ' ...
        'Cannot safely assign labels to connectome rows/columns.']);
end

fprintf('Validated %d MRtrix Desikan-Killiany connectome labels.\n', nNodes);

% Save the exact node order used by all result tables for transparency.
nodeLookupUsed = table(nodeIndex, nodeLabels, nodeFreeSurferNames, ...
    'VariableNames', {'Node','MRtrixLabel','FreeSurferName'});
writetable(nodeLookupUsed, fullfile(outputDir, 'connectome_node_lookup_used.csv'));

%.............................................................................................

%% ========================================================================
% MODIFIED 6: Global metrics calculated from RAW matrices
% =========================================================================
% WHY:
% In the original script global metrics were calculated AFTER dividing each
% matrix by its total strength. That makes global strength and related
% quantities approximately fixed by construction and therefore unsuitable
% for group comparison. Global measures must be calculated BEFORE that
% normalization.

upperIdx = triu(true(nNodes), 1);

globalStrength = zeros(nSub,1);
meanConnectivity = zeros(nSub,1);
density = zeros(nSub,1);
meanNodeStrength = zeros(nSub,1);

for s = 1:nSub
    mat = matsRaw(:,:,s);
    edges = mat(upperIdx);

    globalStrength(s) = sum(edges);             % each undirected edge once
    meanConnectivity(s) = mean(edges);
    density(s) = mean(edges > 0);

    thisNodeStrength = sum(mat, 2);
    meanNodeStrength(s) = mean(thisNodeStrength);
end

metricsTable = table(subjectStrings, groupCat, age, gender, ...
    globalStrength, meanConnectivity, density, meanNodeStrength, ...
    'VariableNames', {'Subject','Group','Age','Gender', ...
    'GlobalStrength','MeanConnectivity','Density','MeanNodeStrength'});

writetable(metricsTable, fullfile(outputDir, 'global_connectome_metrics.csv'));

%.............................................................................................
%% ========================================================================
% MODIFIED 7: Global group comparisons
% =========================================================================
% WHY:
% We retain simple Welch tests for a transparent group comparison, but also
% report age/gender-adjusted linear-model group effects. Strongly skewed
% strength measures are log1p-transformed for the adjusted model.

fprintf('\nGLOBAL CONNECTOME METRICS\n');
fprintf('-------------------------\n');

% GlobalStrength, MeanConnectivity and MeanNodeStrength are algebraically
% related when every subject has the same number of nodes. Therefore only
% GlobalStrength and Density are formally tested here; the other two remain
% in the exported descriptive table.
metricNames = {'GlobalStrength','Density'};
globalResults = table();

for m = 1:length(metricNames)
    metricName = metricNames{m};
    metric = metricsTable.(metricName);

    [~,pWelch,~,statsWelch] = ttest2(metric(isPatient), metric(isCtrl), ...
        'Vartype', 'unequal');

    % Log-transform positive strength-like measures for the regression.
    if strcmp(metricName, 'Density')
        yModel = metric;
    else
        yModel = log1p(metric);
    end

    Tmodel = table(yModel, groupCat, age, gender, ...
        'VariableNames', {'Metric','Group','Age','Gender'});
    mdl = fitlm(Tmodel, 'Metric ~ Group + Age + Gender');

    [betaGroup, tGroup, pAdjusted] = getGroupEffect(mdl);

    fprintf('%s:\n', metricName);
    fprintf('  Control mean = %.6g; schizophrenia mean = %.6g\n', ...
        mean(metric(isCtrl)), mean(metric(isPatient)));
    fprintf('  Welch t = %.3f, p = %.5f\n', statsWelch.tstat, pWelch);
    fprintf('  Adjusted group beta = %.6g, t = %.3f, p = %.5f\n', ...
        betaGroup, tGroup, pAdjusted);

    newRow = table(string(metricName), mean(metric(isCtrl)), ...
        mean(metric(isPatient)), statsWelch.tstat, pWelch, ...
        betaGroup, tGroup, pAdjusted, ...
        'VariableNames', {'Metric','MeanControl','MeanSchizophrenia', ...
        'WelchT','WelchP','AdjustedGroupBeta','AdjustedT','AdjustedP'});
    globalResults = [globalResults; newRow]; %#ok<AGROW>
end

writetable(globalResults, fullfile(outputDir, 'global_metric_statistics.csv'));

%.............................................................................................

%% ========================================================================
% MODIFIED 8: Descriptive group mean matrices
% =========================================================================
% WHY:
% Raw differences, log-transformed differences, and total-normalized
% relative differences answer slightly different questions, so they are
% written separately rather than using one generic "difference" matrix.

meanCtrlRaw = mean(matsRaw(:,:,isCtrl), 3);
meanPatientRaw = mean(matsRaw(:,:,isPatient), 3);
diffRaw = meanPatientRaw - meanCtrlRaw;

meanCtrlLog = mean(matsLog(:,:,isCtrl), 3);
meanPatientLog = mean(matsLog(:,:,isPatient), 3);
diffLog = meanPatientLog - meanCtrlLog;

meanCtrlNorm = mean(matsNorm(:,:,isCtrl), 3);
meanPatientNorm = mean(matsNorm(:,:,isPatient), 3);
diffNorm = meanPatientNorm - meanCtrlNorm;

writematrix(meanCtrlRaw, fullfile(outputDir, 'mean_connectome_controls_RAW.csv'));
writematrix(meanPatientRaw, fullfile(outputDir, 'mean_connectome_schizophrenia_RAW.csv'));
writematrix(diffRaw, fullfile(outputDir, 'difference_scz_minus_ctrl_RAW.csv'));

writematrix(meanCtrlLog, fullfile(outputDir, 'mean_connectome_controls_LOG1P.csv'));
writematrix(meanPatientLog, fullfile(outputDir, 'mean_connectome_schizophrenia_LOG1P.csv'));
writematrix(diffLog, fullfile(outputDir, 'difference_scz_minus_ctrl_LOG1P.csv'));

writematrix(meanCtrlNorm, fullfile(outputDir, 'mean_connectome_controls_NORMALIZED.csv'));
writematrix(meanPatientNorm, fullfile(outputDir, 'mean_connectome_schizophrenia_NORMALIZED.csv'));
writematrix(diffNorm, fullfile(outputDir, 'difference_scz_minus_ctrl_NORMALIZED.csv')); 

%.............................................................................................

%% ========================================================================
% NEW 9: Primary edge-wise statistical analysis
% =========================================================================
% PRIMARY DATA: log1p(raw SIFT2 weights)
%
% For each undirected edge:
%   A) calculate prevalence (fraction of subjects with non-zero raw weight)
%   B) retain edges with prevalence >= minEdgePresence
%   C) calculate an unadjusted Welch group comparison
%   D) calculate Hedges' g effect size on log1p weights
%   E) fit: Connectivity ~ Group + Age + Gender
%   F) FDR-correct the adjusted GROUP p-values across tested edges
%
% IMPORTANT INTERPRETATION:
% Positive group beta / Hedges g = schizophrenia > controls.
% Negative group beta / Hedges g = schizophrenia < controls.

[rowIdx, colIdx] = find(upperIdx);
nEdges = numel(rowIdx);

edgePresence = zeros(nEdges,1);
meanCtrlEdge = nan(nEdges,1);
meanPatientEdge = nan(nEdges,1);
welchT = nan(nEdges,1);
welchP = nan(nEdges,1);
hedgesG = nan(nEdges,1);
adjustedBeta = nan(nEdges,1);
adjustedT = nan(nEdges,1);
adjustedP = nan(nEdges,1);
edgeTested = false(nEdges,1);

for e = 1:nEdges
    i = rowIdx(e);
    j = colIdx(e);

    rawAll = squeeze(matsRaw(i,j,:));
    edgePresence(e) = mean(rawAll > 0);

    if edgePresence(e) < minEdgePresence
        continue;
    end

    edgeTested(e) = true;

    y = squeeze(matsLog(i,j,:));
    ctrl = y(isCtrl);
    patient = y(isPatient);

    meanCtrlEdge(e) = mean(ctrl);
    meanPatientEdge(e) = mean(patient);

    % A) Unadjusted Welch test (schizophrenia minus control direction).
    [~,welchP(e),~,statsWelch] = ttest2(patient, ctrl, 'Vartype','unequal');
    welchT(e) = statsWelch.tstat;

    % B) Hedges' g effect size on log1p edge weights.
    hedgesG(e) = hedges_g(patient, ctrl);

    % C) Primary adjusted model.
    Tedge = table(y, groupCat, age, gender, ...
        'VariableNames', {'Connectivity','Group','Age','Gender'});
    mdl = fitlm(Tedge, 'Connectivity ~ Group + Age + Gender');

    [adjustedBeta(e), adjustedT(e), adjustedP(e)] = getGroupEffect(mdl);
end

% FDR correction ONLY across edges that were actually tested.
adjustedQ = bh_fdr(adjustedP);
welchQ = bh_fdr(welchP);   % secondary/unadjusted FDR result, also exported

sigAdjusted = adjustedQ < fdrAlpha;

region1 = nodeLabels(rowIdx);
region2 = nodeLabels(colIdx);
region1FreeSurfer = nodeFreeSurferNames(rowIdx);
region2FreeSurfer = nodeFreeSurferNames(colIdx);

edgeResults = table(region1, region2, region1FreeSurfer, region2FreeSurfer, rowIdx, colIdx, ...
    edgePresence, edgeTested, ...
    meanCtrlEdge, meanPatientEdge, ...
    welchT, welchP, welchQ, hedgesG, ...
    adjustedBeta, adjustedT, adjustedP, adjustedQ, sigAdjusted, ...
    'VariableNames', {'Region1','Region2','Region1_FreeSurfer','Region2_FreeSurfer','Node1','Node2', ...
    'PresenceFraction','Tested','MeanControl_LOG1P','MeanSchizophrenia_LOG1P', ...
    'WelchT','WelchP','WelchFDR_Q','HedgesG', ...
    'AdjustedGroupBeta','AdjustedT','AdjustedP','AdjustedFDR_Q','FDR_Significant'});

% Sort with the most negative adjusted group effects first.
edgeResults = sortrows(edgeResults, 'AdjustedGroupBeta', 'ascend');
writetable(edgeResults, fullfile(outputDir, 'edge_statistics_PRIMARY_LOG1P.csv'));

fprintf('\nPRIMARY EDGE-WISE ANALYSIS\n');
fprintf('--------------------------\n');
fprintf('Possible undirected edges: %d\n', nEdges);
fprintf('Edges tested (presence >= %.0f%%): %d\n', ...
    100*minEdgePresence, sum(edgeTested));
fprintf('Edges significant after adjusted-model FDR q < %.3f: %d\n', ...
    fdrAlpha, sum(sigAdjusted));

%.............................................................................................

%% ========================================================================
% NEW 10: Edge-wise effect-size and significance matrices
% =========================================================================
% WHY:
% Raw difference maps can be dominated by the absolute scale of different
% connections. Hedges' g provides a standardized effect size. A second map
% shows only edges surviving FDR correction in the adjusted model.

gMat = zeros(nNodes);
betaMat = zeros(nNodes);
qMat = nan(nNodes);
sigGMat = zeros(nNodes);

for e = 1:nEdges
    i = rowIdx(e);
    j = colIdx(e);

    if edgeTested(e)
        gMat(i,j) = hedgesG(e);
        gMat(j,i) = hedgesG(e);

        betaMat(i,j) = adjustedBeta(e);
        betaMat(j,i) = adjustedBeta(e);

        qMat(i,j) = adjustedQ(e);
        qMat(j,i) = adjustedQ(e);

        if sigAdjusted(e)
            sigGMat(i,j) = hedgesG(e);
            sigGMat(j,i) = hedgesG(e);
        end
    end
end

writematrix(gMat, fullfile(outputDir, 'edge_HedgesG_matrix.csv'));
writematrix(betaMat, fullfile(outputDir, 'edge_adjusted_group_beta_matrix.csv'));
writematrix(qMat, fullfile(outputDir, 'edge_adjusted_FDR_q_matrix.csv'));
writematrix(sigGMat, fullfile(outputDir, 'edge_HedgesG_FDRsignificant_only.csv'));


%.............................................................................................

%% ========================================================================
% NEW 11: Node-strength analysis
% =========================================================================
% WHY:
% Edge-wise analysis asks which specific connections differ. Node strength
% asks whether a REGION has altered total connectivity across all of its
% connections. This is a simple and useful complement, particularly when
% several abnormal edges converge on the same region.
%
% Raw node strength is summed first, then log1p-transformed for statistics.

nodeStrengthRaw = zeros(nSub,nNodes);
nodeStrengthLog = zeros(nSub,nNodes);

for s = 1:nSub
    nodeStrengthRaw(s,:) = sum(matsRaw(:,:,s), 2)';
end
nodeStrengthLog = log1p(nodeStrengthRaw);

nodeMeanCtrl = zeros(nNodes,1);
nodeMeanPatient = zeros(nNodes,1);
nodeHedgesG = zeros(nNodes,1);
nodeWelchT = zeros(nNodes,1);
nodeWelchP = zeros(nNodes,1);
nodeAdjustedBeta = zeros(nNodes,1);
nodeAdjustedT = zeros(nNodes,1);
nodeAdjustedP = zeros(nNodes,1);

for n = 1:nNodes
    y = nodeStrengthLog(:,n);
    ctrl = y(isCtrl);
    patient = y(isPatient);

    nodeMeanCtrl(n) = mean(ctrl);
    nodeMeanPatient(n) = mean(patient);

    [~,nodeWelchP(n),~,statsWelch] = ttest2(patient, ctrl, 'Vartype','unequal');
    nodeWelchT(n) = statsWelch.tstat;
    nodeHedgesG(n) = hedges_g(patient, ctrl);

    Tnode = table(y, groupCat, age, gender, ...
        'VariableNames', {'NodeStrength','Group','Age','Gender'});
    mdl = fitlm(Tnode, 'NodeStrength ~ Group + Age + Gender');

    [nodeAdjustedBeta(n), nodeAdjustedT(n), nodeAdjustedP(n)] = getGroupEffect(mdl);
end

nodeWelchQ = bh_fdr(nodeWelchP);
nodeAdjustedQ = bh_fdr(nodeAdjustedP);
nodeSigAdjusted = nodeAdjustedQ < fdrAlpha;

nodeResults = table(nodeLabels, nodeFreeSurferNames, (1:nNodes)', ...
    nodeMeanCtrl, nodeMeanPatient, nodeHedgesG, ...
    nodeWelchT, nodeWelchP, nodeWelchQ, ...
    nodeAdjustedBeta, nodeAdjustedT, nodeAdjustedP, nodeAdjustedQ, ...
    nodeSigAdjusted, ...
    'VariableNames', {'Region','FreeSurferName','Node','MeanControl_LOG1PStrength', ...
    'MeanSchizophrenia_LOG1PStrength','HedgesG','WelchT','WelchP','WelchFDR_Q', ...
    'AdjustedGroupBeta','AdjustedT','AdjustedP','AdjustedFDR_Q','FDR_Significant'});

nodeResults = sortrows(nodeResults, 'AdjustedGroupBeta', 'ascend');
writetable(nodeResults, fullfile(outputDir, 'node_strength_statistics.csv'));

fprintf('\nNODE-STRENGTH ANALYSIS\n');
fprintf('----------------------\n');
fprintf('Regions significant after adjusted-model FDR q < %.3f: %d / %d\n', ...
    fdrAlpha, sum(nodeSigAdjusted), nNodes);


%.............................................................................................
%% ========================================================================
% NEW 12: Sensitivity analysis using total-normalized connectomes
% =========================================================================
% WHY:
% The original normalized analysis is retained as a useful SECONDARY test.
% Here, a positive effect means an edge occupies a larger FRACTION of total
% connectome weight in schizophrenia; it should not automatically be called
% an absolute increase in structural connectivity.
%
% We repeat the same adjusted edge model and FDR correction, but do not
% duplicate all figures/statistics. The output can be compared with the
% primary log1p(raw SIFT2) result for robustness.

normBeta = nan(nEdges,1);
normT = nan(nEdges,1);
normP = nan(nEdges,1);
normHedgesG = nan(nEdges,1);

for e = 1:nEdges
    if ~edgeTested(e)
        continue;
    end

    i = rowIdx(e);
    j = colIdx(e);

    y = squeeze(matsNorm(i,j,:));
    ctrl = y(isCtrl);
    patient = y(isPatient);

    normHedgesG(e) = hedges_g(patient, ctrl);

    Tedge = table(y, groupCat, age, gender, ...
        'VariableNames', {'Connectivity','Group','Age','Gender'});
    mdl = fitlm(Tedge, 'Connectivity ~ Group + Age + Gender');

    [normBeta(e), normT(e), normP(e)] = getGroupEffect(mdl);
end

normQ = bh_fdr(normP);
normSig = normQ < fdrAlpha;

sensitivityResults = table(region1, region2, region1FreeSurfer, region2FreeSurfer, rowIdx, colIdx, ...
    edgePresence, edgeTested, normHedgesG, normBeta, normT, normP, normQ, normSig, ...
    'VariableNames', {'Region1','Region2','Region1_FreeSurfer','Region2_FreeSurfer','Node1','Node2','PresenceFraction','Tested', ...
    'HedgesG_NORMALIZED','AdjustedGroupBeta_NORMALIZED','AdjustedT_NORMALIZED', ...
    'AdjustedP_NORMALIZED','AdjustedFDR_Q_NORMALIZED','FDR_Significant_NORMALIZED'});

sensitivityResults = sortrows(sensitivityResults, ...
    'AdjustedGroupBeta_NORMALIZED', 'ascend');
writetable(sensitivityResults, ...
    fullfile(outputDir, 'edge_statistics_SENSITIVITY_TOTAL_NORMALIZED.csv'));

% NEW: Compact comparison between primary and normalized results.
comparisonResults = table(region1, region2, region1FreeSurfer, region2FreeSurfer, rowIdx, colIdx, edgePresence, ...
    adjustedBeta, adjustedQ, sigAdjusted, hedgesG, ...
    normBeta, normQ, normSig, normHedgesG, ...
    'VariableNames', {'Region1','Region2','Region1_FreeSurfer','Region2_FreeSurfer','Node1','Node2','PresenceFraction', ...
    'PrimaryBeta_LOG1P','PrimaryFDR_Q','PrimarySignificant','PrimaryHedgesG', ...
    'NormalizedBeta','NormalizedFDR_Q','NormalizedSignificant','NormalizedHedgesG'});

% Flag effects that point in the same direction in both representations.
comparisonResults.SameDirection = ...
    sign(comparisonResults.PrimaryBeta_LOG1P) == sign(comparisonResults.NormalizedBeta);
comparisonResults.SignificantBoth = ...
    comparisonResults.PrimarySignificant & comparisonResults.NormalizedSignificant;

comparisonResults = sortrows(comparisonResults, 'PrimaryFDR_Q', 'ascend');
writetable(comparisonResults, ...
    fullfile(outputDir, 'edge_primary_vs_normalized_comparison.csv'));

fprintf('\nSENSITIVITY ANALYSIS\n');
fprintf('--------------------\n');
fprintf('FDR-significant edges using total-normalized matrices: %d\n', sum(normSig));
fprintf('FDR-significant in BOTH primary and normalized analyses: %d\n', ...
    sum(comparisonResults.SignificantBoth));

%............................................................................................. 


%.............................................................................................

%% ========================================================================
% MODIFIED 13: Figures
% =========================================================================
% Descriptive matrix figures are retained. Additional standardized effect
% and FDR-significant effect figures are added. Figures are also exported at
% 300 dpi for convenient inspection and reporting.

% ----- Raw mean controls -----
fig = figure;
imagesc(meanCtrlRaw);
axis image;
colorbar;
title('Mean connectome: controls (raw SIFT2)');
%exportgraphics(fig, fullfile(outputDir, 'mean_controls_RAW.png'), 'Resolution', 300);

labels_table = readtable('Copy_2_of_label_fs_default.txt'); 
area_names = labels_table.Var2; 

figure('Position', [100, 100, 800, 700]); % Hacemos la figura un poco más grande
imagesc(meanCtrlRaw); % O meanPatient
colorbar;
title('\fontsize{9}Mean connectome: controls (raw SIFT2)');

n = length(area_names);
xticks(1:n);
xticklabels(area_names);
xtickangle(45); % Rotar a 45 grados es mejor que 90 para la lectura

yticks(1:n);
yticklabels(area_names);

% 4. Ajuste fino de estilo
set(gca, 'FontSize', 5); % Un tamaño pequeño pero legible
axis image; % Mantiene la matriz cuadrada
grid on; % Opcional: ayuda a seguir las filas/columnas
set(gca, 'TickLength', [0 0]);

%............................................................................................. 

% ----- Raw mean schizophrenia -----
fig = figure;
imagesc(meanPatientRaw);
axis image;
colorbar;
title('Mean connectome: schizophrenia (raw SIFT2)');
%exportgraphics(fig, fullfile(outputDir, 'mean_schizophrenia_RAW.png'), 'Resolution', 300);

figure('Position', [100, 100, 800, 700]); % Hacemos la figura un poco más grande
imagesc(meanPatientRaw); % O meanPatient
colorbar;
title('\fontsize{9}Mean: schizophrenia (raw SIFT2).png');

n = length(area_names);
xticks(1:n);
xticklabels(area_names);
xtickangle(45); % Rotar a 45 grados es mejor que 90 para la lectura

yticks(1:n);
yticklabels(area_names);

% 4. Ajuste fino de estilo
set(gca, 'FontSize', 5); % Un tamaño pequeño pero legible
axis image; % Mantiene la matriz cuadrada
grid on; % Opcional: ayuda a seguir las filas/columnas
set(gca, 'TickLength', [0 0]);

%.............................................................................................  
% ----- Raw schizophrenia - control -----
fig = figure;
maxAbsRaw = max(abs(diffRaw(:)));
if maxAbsRaw == 0, maxAbsRaw = 1; end
imagesc(diffRaw, [-maxAbsRaw maxAbsRaw]);
axis image;
colorbar;
title('Difference: schizophrenia - controls (raw SIFT2)');
%exportgraphics(fig, fullfile(outputDir, 'difference_RAW.png'), 'Resolution', 300);

figure('Position', [100, 100, 800, 700]); % Hacemos la figura un poco más grande
imagesc(diffRaw, [-maxAbsRaw maxAbsRaw]); % O meanPatient
colorbar;
title('\fontsize{9}Difference: schizophrenia - controls (raw SIFT2');

n = length(area_names);
xticks(1:n);
xticklabels(area_names);
xtickangle(45); % Rotar a 45 grados es mejor que 90 para la lectura

yticks(1:n);
yticklabels(area_names);

% 4. Ajuste fino de estilo
set(gca, 'FontSize', 5); % Un tamaño pequeño pero legible
axis image; % Mantiene la matriz cuadrada
grid on; % Opcional: ayuda a seguir las filas/columnas
set(gca, 'TickLength', [0 0]);


%............................................................................................. 

% ----- Total-normalized relative difference -----
fig = figure;
maxAbsNorm = max(abs(diffNorm(:)));
if maxAbsNorm == 0, maxAbsNorm = 1; end
imagesc(diffNorm, [-maxAbsNorm maxAbsNorm]);
axis image;
colorbar;
title('Relative difference: schizophrenia - controls (total-normalized)');
%exportgraphics(fig, fullfile(outputDir, 'difference_TOTAL_NORMALIZED.png'), 'Resolution', 300);

figure('Position', [100, 100, 800, 700]); % Hacemos la figura un poco más grande
imagesc(diffNorm, [-maxAbsNorm maxAbsNorm]); % O meanPatient
colorbar;
title('\fontsize{9}Relative difference: schizophrenia - controls (total-normalized)');

n = length(area_names);
xticks(1:n);
xticklabels(area_names);
xtickangle(45); % Rotar a 45 grados es mejor que 90 para la lectura

yticks(1:n);
yticklabels(area_names);

% 4. Ajuste fino de estilo
set(gca, 'FontSize', 5); % Un tamaño pequeño pero legible
axis image; % Mantiene la matriz cuadrada
grid on; % Opcional: ayuda a seguir las filas/columnas
set(gca, 'TickLength', [0 0]);

%............................................................................................. 
% ----- Hedges' g for all tested edges -----
fig = figure;
maxAbsG = max(abs(gMat(:)));
if maxAbsG == 0, maxAbsG = 1; end
imagesc(gMat, [-maxAbsG maxAbsG]);
axis image;
colorbar;
title('Edge effect size: schizophrenia - controls (Hedges'' g)');
%exportgraphics(fig, fullfile(outputDir, 'edge_HedgesG.png'), 'Resolution', 300);


figure('Position', [100, 100, 800, 700]); % Hacemos la figura un poco más grande
imagesc(gMat, [-maxAbsG maxAbsG]); % O meanPatient
colorbar;
title('\fontsize{9}Edge effect size: schizophrenia - controls (Hedges'' g)');

targetRegions = {'L.TTG', 'L.PA', 'R.TTG', 'R.PA'}; 


targetIndices = [];
for k = 1:length(targetRegions)
    idx = find(strcmp(area_names, targetRegions{k}));
    if ~isempty(idx)
        targetIndices = [targetIndices, idx];
    end
end

hold on; % Para no sobrescribir la matriz
for idx = targetIndices
    % Líneas verticales (eje X)
    xline(idx, '-', 'Color', [0.85 0.15 0.15], 'LineWidth', 1.0); 

    % Líneas horizontales (eje Y)
    yline(idx, '-', 'Color', [0.85 0.15 0.15], 'LineWidth', 1.0); 
end
hold off;

% 2. Supongamos que 'regionNames' es la celda (cell array) con las 84 etiquetas
%    (Asegúrate de reemplazar 'regionNames' por el nombre exacto de tu variable de etiquetas)
labelsX = area_names; 
labelsY = area_names;


% 3. Recorrer la lista y reemplazar por texto vacío '' las que no sean de interés
for k = 1:length(area_names)
    if ~ismember(area_names{k}, targetRegions)
        labelsX{k} = '';
        labelsY{k} = '';
    end
end


% 4. Aplicar los nuevos nombres filtrados a la figura actual
ax = gca;
ax.XTick = 1:length(area_names);
ax.YTick = 1:length(area_names);
ax.XTickLabel = labelsX;
ax.YTickLabel = labelsY; 


% (Opcional) Ajustar el tamaño y rotación de las 4 etiquetas visibles para que destaquen
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
ax.XTickLabelRotation = 45; % Rotar etiquetas del eje X si es necesario 

%............................................................................................. 
% ----- Hedges' g only for FDR-significant adjusted edges -----
fig = figure;
maxAbsSigG = max(abs(sigGMat(:)));
if maxAbsSigG == 0, maxAbsSigG = 1; end
imagesc(sigGMat, [-maxAbsSigG maxAbsSigG]);
axis image;
colorbar;
title(sprintf('FDR-significant edges only (q < %.2f): Hedges'' g', fdrAlpha));
%exportgraphics(fig, fullfile(outputDir, 'edge_HedgesG_FDRsignificant.png'), 'Resolution', 300);

figure('Position', [100, 100, 800, 700]); % Hacemos la figura un poco más grande
imagesc(sigGMat, [-maxAbsSigG maxAbsSigG]); % O meanPatient
colorbar;
title('\fontsize{9}Edge effect size: schizophrenia - controls (Hedges'' g)');

n = length(area_names);
xticks(1:n);
xticklabels(area_names);
xtickangle(45); % Rotar a 45 grados es mejor que 90 para la lectura

yticks(1:n);
yticklabels(area_names);

% 4. Ajuste fino de estilo
set(gca, 'FontSize', 5); % Un tamaño pequeño pero legible
axis image; % Mantiene la matriz cuadrada
grid on; % Opcional: ayuda a seguir las filas/columnas
set(gca, 'TickLength', [0 0]);


%............................................................................................. 

%% ========================================================================
% NEW 14: Print strongest effects for quick inspection
% =========================================================================
% The CSV files remain the authoritative results. This simply prints the
% strongest tested negative and positive adjusted effects to the console.

nShow = 10;
testedResults = edgeResults(edgeResults.Tested, :);

fprintf('\nSTRONGEST NEGATIVE ADJUSTED EDGE EFFECTS\n');
fprintf('----------------------------------------\n');
disp(testedResults(1:min(nShow,height(testedResults)), ...
    {'Region1','Region2','AdjustedGroupBeta','AdjustedP','AdjustedFDR_Q','HedgesG'}));

fprintf('\nSTRONGEST POSITIVE ADJUSTED EDGE EFFECTS\n');
fprintf('----------------------------------------\n');
positiveOrder = sortrows(testedResults, 'AdjustedGroupBeta', 'descend');
disp(positiveOrder(1:min(nShow,height(positiveOrder)), ...
    {'Region1','Region2','AdjustedGroupBeta','AdjustedP','AdjustedFDR_Q','HedgesG'}));

fprintf('\nAnalysis complete. Outputs written to:\n%s\n', outputDir);

%% ========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================

function g = hedges_g(group1, group2)
% HEDGES_G Standardized mean difference with small-sample correction.
% Direction is group1 - group2. In this script:
%   group1 = schizophrenia
%   group2 = controls

    group1 = group1(isfinite(group1));
    group2 = group2(isfinite(group2));

    n1 = numel(group1);
    n2 = numel(group2);

    if n1 < 2 || n2 < 2
        g = NaN;
        return;
    end

    pooledVar = ((n1-1)*var(group1) + (n2-1)*var(group2)) / (n1+n2-2);

    if pooledVar <= 0 || ~isfinite(pooledVar)
        % If both groups are constant and identical, effect is zero.
        if mean(group1) == mean(group2)
            g = 0;
        else
            g = NaN;
        end
        return;
    end

    d = (mean(group1) - mean(group2)) / sqrt(pooledVar);
    J = 1 - 3 / (4*(n1+n2) - 9);
    g = J * d;
end

function q = bh_fdr(p)
% BH_FDR Benjamini-Hochberg FDR-adjusted q-values.
% NaN p-values are left as NaN and are not included in the correction.

    q = nan(size(p));
    valid = isfinite(p);

    if ~any(valid)
        return;
    end

    pValid = p(valid);
    [pSorted, order] = sort(pValid(:), 'ascend');
    m = numel(pSorted);

    qSorted = pSorted .* m ./ (1:m)';

    % Enforce the required monotonic ordering from largest rank backwards.
    qSorted = flipud(cummin(flipud(qSorted)));
    qSorted(qSorted > 1) = 1;

    qValid = nan(m,1);
    qValid(order) = qSorted;
    q(valid) = qValid;
end

function fsDefaultFile = find_mrtrix_fs_default(manualPath)
% FIND_MRTRIX_FS_DEFAULT Locate MRtrix's default FreeSurfer DK lookup table.
%
% If manualPath is non-empty it is used directly. Otherwise common MRtrix
% installation layouts are searched, including a path inferred from the
% installed labelconvert executable.

    candidates = strings(0,1);

    if nargin >= 1 && strlength(strtrim(string(manualPath))) > 0
        candidates(end+1,1) = string(manualPath);
    end

    % Infer MRtrix prefix from labelconvert on the current PATH.
    [status, labelconvertExe] = system('command -v labelconvert');
    if status == 0
        labelconvertExe = strtrim(labelconvertExe);
        if ~isempty(labelconvertExe)
            binDir = fileparts(labelconvertExe);
            candidates(end+1,1) = string(fullfile(binDir, '..', 'share', ...
                'mrtrix3', 'labelconvert', 'fs_default.txt'));
        end
    end

    % Optional MRtrix installation root if defined by the local system.
    mrtrixHome = getenv('MRTRIX3_HOME');
    if ~isempty(mrtrixHome)
        candidates(end+1,1) = string(fullfile(mrtrixHome, 'share', ...
            'mrtrix3', 'labelconvert', 'fs_default.txt'));
        candidates(end+1,1) = string(fullfile(mrtrixHome, 'labelconvert', ...
            'fs_default.txt'));
    end

    % Common Linux / macOS installation locations.
    candidates = [candidates; ...
        "/usr/share/mrtrix3/labelconvert/fs_default.txt"; ...
        "/usr/local/share/mrtrix3/labelconvert/fs_default.txt"; ...
        "/opt/homebrew/share/mrtrix3/labelconvert/fs_default.txt"];

    % Remove duplicates while preserving order.
    candidates = unique(candidates, 'stable');

    for k = 1:numel(candidates)
        candidate = char(candidates(k));
        if isfile(candidate)
            fsDefaultFile = candidate;
            return;
        end
    end

    error(['Could not locate MRtrix fs_default.txt automatically. ' ...
        'Set mrtrixFsDefaultFile manually near the top of the script.\n' ...
        'Expected file: share/mrtrix3/labelconvert/fs_default.txt']);
end

function [beta, tValue, pValue] = getGroupEffect(mdl)
% GETGROUPEFFECT Extract the schizophrenia-vs-control coefficient from
% fitlm when Group is categorical with Control as the reference category.

    rowNames = string(mdl.Coefficients.Properties.RowNames);

    idx = contains(rowNames, 'Group') & contains(rowNames, 'Schizophrenia');

    if sum(idx) ~= 1
        error(['Could not uniquely identify the schizophrenia group coefficient ' ...
            'in the fitted model. Coefficient names were: %s'], ...
            strjoin(rowNames, ', '));
    end

    beta = mdl.Coefficients.Estimate(idx);
    tValue = mdl.Coefficients.tStat(idx);
    pValue = mdl.Coefficients.pValue(idx);
end



