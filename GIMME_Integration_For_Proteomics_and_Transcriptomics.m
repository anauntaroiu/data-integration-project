
%% Initialize COBRA
initCobraToolbox % this only needs to be run once

%% import proteomics data 

% Set working directory to where data is saved
cd 'C:\Users\Ana\Documents\MATLAB\Capstone\Combined Proteomics'

% Load proteomics data
proteomics_I = readtable('Proteomics_Data_For_GIMME_Sample_I.csv','Delimiter',',');
proteomics_II = readtable('Proteomics_Data_For_GIMME_Sample_II.csv','Delimiter',',');
proteomics_III = readtable('Proteomics_Data_For_GIMME_Sample_III.csv','Delimiter',',');
proteomics_IV = readtable('Proteomics_Data_For_GIMME_Sample_IV.csv','Delimiter',',');
proteomics_V = readtable('Proteomics_Data_For_GIMME_Sample_V.csv','Delimiter',',');
medians_proteomics = readtable('Medians_for_Complete_Proteomics.csv','Delimiter',',');

%% import transcriptomics data

% Set working directory to where data is saved
cd 'C:\Users\Ana\Documents\MATLAB\Capstone\Transcriptomics'

% Load transcriptomics data
transcriptomics_I = readtable('Transcriptomics_Data_For_GIMME_Sample_I.csv','Delimiter',',');
transcriptomics_II = readtable('Transcriptomics_Data_For_GIMME_Sample_II.csv','Delimiter',',');
transcriptomics_III = readtable('Transcriptomics_Data_For_GIMME_Sample_III.csv','Delimiter',',');
transcriptomics_IV = readtable('Transcriptomics_Data_For_GIMME_Sample_IV.csv','Delimiter',',');
transcriptomics_V = readtable('Transcriptomics_Data_For_GIMME_Sample_V.csv','Delimiter',',');
medians_transcriptomics = readtable('Medians_for_Complete_Transcriptomics.csv','Delimiter',',');

%% import model 

% Set working directory to find model
cd 'C:\Users\Ana\Documents\MATLAB\Capstone'

model = readCbModel('SA_USA300_FPR3757.xml');
model = creategrRulesField(model); % Need this field to delete genes
optimizeCbModel(model) % make sure the model is working

%% Remove Genes with No Corresponding Reactions
reaction_list = struct2cell(findRxnsFromGenes(model,...
    model.genes)); % Find corresponding reactions for all genes in model

gene_names = model.genes; % Corresponding Genes for the reaction list
genes_to_remove = cell(1, 1); % Create variable for genes to be removed
j = 1; % Index for genes_to_remove variable

% Find genes w/o reactions
for i = 1:length(reaction_list) % Loop through each rxn
    
    if (isempty(reaction_list{i,1})) % If a rxn does not exist
        genes_to_remove(j) = gene_names(i); % record gene name
        j = j + 1;
    end
    
end  

model = deleteModelGenes(model,...
    genes_to_remove); % Delete genes with no corresponding reactions

%% GIMME model generation
% Proteomics Data Integration

% Sample I
gimme_proteomics_I_I = GIMME(model, proteomics_I.Replicate_1,...
    table2array(medians_proteomics(1,1)));
gimme_proteomics_I_II = GIMME(model, proteomics_I.Replicate_2,...
    table2array(medians_proteomics(1,2)));
gimme_proteomics_I_III = GIMME(model, proteomics_I.Replicate_3,...
    table2array(medians_proteomics(1,3)));

% Sample II
gimme_proteomics_II_I = GIMME(model, proteomics_II.Replicate_1,...
    table2array(medians_proteomics(2,1)));
gimme_proteomics_II_II = GIMME(model, proteomics_II.Replicate_2,...
    table2array(medians_proteomics(2,2)));
gimme_proteomics_II_III = GIMME(model, proteomics_II.Replicate_3,...
    table2array(medians_proteomics(2,3)));

% Sample III
gimme_proteomics_III_I = GIMME(model, proteomics_III.Replicate_1,...
    table2array(medians_proteomics(3,1)));
gimme_proteomics_III_II = GIMME(model, proteomics_III.Replicate_2,...
    table2array(medians_proteomics(3,2)));
gimme_proteomics_III_III = GIMME(model, proteomics_III.Replicate_3,...
    table2array(medians_proteomics(3,3)));

% Sample IV
gimme_proteomics_IV_I = GIMME(model, proteomics_IV.Replicate_1,...
    table2array(medians_proteomics(4,1)));
gimme_proteomics_IV_II = GIMME(model, proteomics_IV.Replicate_2,...
    table2array(medians_proteomics(4,2)));
gimme_proteomics_IV_III = GIMME(model, proteomics_IV.Replicate_3,...
    table2array(medians_proteomics(4,3)));

% Sample V
gimme_proteomics_V_I = GIMME(model, proteomics_V.Replicate_1,...
    table2array(medians_proteomics(5,1)));
gimme_proteomics_V_II = GIMME(model, proteomics_V.Replicate_2,...
    table2array(medians_proteomics(5,2)));
gimme_proteomics_V_III = GIMME(model, proteomics_V.Replicate_3,...
    table2array(medians_proteomics(5,3)));

%% GIMME model generation
% Transcriptomics Data Integration

% Sample I
gimme_transcriptomics_I_I = GIMME(model, transcriptomics_I.Replicate_1,...
    table2array(medians_transcriptomics(1,1)));
gimme_transcriptomics_I_II = GIMME(model, transcriptomics_I.Replicate_2,...
    table2array(medians_transcriptomics(1,2)));
gimme_transcriptomics_I_III = GIMME(model, transcriptomics_I.Replicate_3,...
    table2array(medians_transcriptomics(1,3)));

% Sample II
gimme_transcriptomics_II_I = GIMME(model, transcriptomics_II.Replicate_1,...
    table2array(medians_transcriptomics(2,1)));
gimme_transcriptomics_II_II = GIMME(model, transcriptomics_II.Replicate_2,...
    table2array(medians_transcriptomics(2,2)));
gimme_transcriptomics_II_III = GIMME(model, transcriptomics_II.Replicate_3,...
    table2array(medians_transcriptomics(2,3)));

% Sample III
gimme_transcriptomics_III_I = GIMME(model, transcriptomics_III.Replicate_1,...
    table2array(medians_transcriptomics(3,1)));
gimme_transcriptomics_III_II = GIMME(model, transcriptomics_III.Replicate_2,...
    table2array(medians_transcriptomics(3,2)));
gimme_transcriptomics_III_III = GIMME(model, transcriptomics_III.Replicate_3,...
    table2array(medians_transcriptomics(3,3)));

% Sample IV
gimme_transcriptomics_IV_I = GIMME(model, transcriptomics_IV.Replicate_1,...
    table2array(medians_transcriptomics(4,1)));
gimme_transcriptomics_IV_II = GIMME(model, transcriptomics_IV.Replicate_2,...
    table2array(medians_transcriptomics(4,2)));
gimme_transcriptomics_IV_III = GIMME(model, transcriptomics_IV.Replicate_3,...
    table2array(medians_transcriptomics(4,3)));

% Sample V
gimme_transcriptomics_V_I = GIMME(model, transcriptomics_V.Replicate_1,...
    table2array(medians_transcriptomics(5,1)));
gimme_transcriptomics_V_II = GIMME(model, transcriptomics_V.Replicate_2,...
    table2array(medians_transcriptomics(5,2)));
gimme_transcriptomics_V_III = GIMME(model, transcriptomics_V.Replicate_3,...
    table2array(medians_transcriptomics(5,3)));

 %% essentiality predictions 
% singleRxnDeletion - deletes a reaction, sees how the model grows, and
% compares to previous growth

% Store the proteomics and transcriptomics models
proteomics_models = {gimme_proteomics_I_I;gimme_proteomics_I_II;gimme_proteomics_I_III;...
    gimme_proteomics_II_I;gimme_proteomics_II_II;gimme_proteomics_II_III;...
    gimme_proteomics_III_I;gimme_proteomics_III_II;gimme_proteomics_III_III;...
    gimme_proteomics_IV_I;gimme_proteomics_IV_II;gimme_proteomics_IV_III;...
    gimme_proteomics_V_I;gimme_proteomics_V_II;gimme_proteomics_V_III};

transcriptomics_models = {gimme_transcriptomics_I_I;gimme_transcriptomics_I_II;gimme_transcriptomics_I_III;...
    gimme_transcriptomics_II_I;gimme_transcriptomics_II_II;gimme_transcriptomics_II_III;...
    gimme_transcriptomics_III_I;gimme_transcriptomics_III_II;gimme_transcriptomics_III_III;...
    gimme_transcriptomics_IV_I;gimme_transcriptomics_IV_II;gimme_transcriptomics_IV_III;...
    gimme_transcriptomics_V_I;gimme_transcriptomics_V_II;gimme_transcriptomics_V_III};

% Conduct Single Reaction Deletions on all Transcriptomics Models
lethal_reactions_transcriptomics = cell(1,15); % each column = different model
for n = 1:15 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(transcriptomics_models{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(transcriptomics_models{n},'FBA',...
     transcriptomics_models{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio)) % If growth ratio is NA, set it to zero
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_transcriptomics(j,n) = transcriptomics_models{n}.rxns(i);
       j = j + 1;
    end
end
end

clearvars n j i grRatio 

% Conduct Single Reaction Deletions on all Proteomics Models
lethal_reactions_proteomics = cell(1,15); % each column = different model
for n = 1:15 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(proteomics_models{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(proteomics_models{n},'FBA',...
     proteomics_models{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_proteomics(j,n) = proteomics_models{n}.rxns(i);
       j = j + 1;
    end
end
end

clearvars n j i grRatio 

% Conduct Single Reaction Deletions for Unconstrained Model
lethal_reactions_unconstrained = cell(1,1);
j = 1;
for i = 1:length(model.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(model,'FBA',...
     model.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_unconstrained(j) = model.rxns(i);
       j = j + 1;
    end
end

%% Save Results

% Set working directory to where data should be saved
cd 'C:\Users\Ana\Documents\Capstone\R Code\Essential Reaction Overlap\Complete Data'

% Save data
P = table(lethal_reactions_proteomics);
writetable(P,'Essential_Rxns_Proteomics_Models.csv','Delimiter',',');

T = table(lethal_reactions_transcriptomics);
writetable(T,'Essential_Rxns_Transcriptomics_Models.csv','Delimiter',',');

U = table(lethal_reactions_unconstrained');
writetable(U,'Essential_Rxns_Unconstrained_Model.csv','Delimiter',',');

