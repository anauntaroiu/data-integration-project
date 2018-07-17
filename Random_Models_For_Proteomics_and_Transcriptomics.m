
%% Initialize COBRA
initCobraToolbox % this only needs to be run once

%% import data for transcriptomics models

% Set working directory to where data is saved
cd 'C:\Users\Ana\Documents\MATLAB\Capstone\Transcriptomics'

% Load random data for transcriptomics
Random_Transcriptomics = readtable('Random_Data_For_GIMME_Transcriptomics.csv','Delimiter',',');
Medians_Random_Transcriptomics = readtable('Medians_for_Transcriptomics_Random.csv','Delimiter',',');

%% import data for proteomics models

% Set working directory to where data is saved
cd 'C:\Users\Ana\Documents\MATLAB\Capstone\Combined Proteomics'

% Load random data for proteomics 
Random_Proteomics = readtable('Random_Data_For_GIMME_Proteomics.csv','Delimiter',',');
Medians_Random_Proteomics = readtable('Medians_for_Proteomics_Random.csv','Delimiter',',');

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
rxns_to_remove = cell(1, 1); % Create variable for genes to be removed
j = 1; % Index for genes_to_remove variable

% Find genes w/o reactions
for i = 1:length(reaction_list) % Loop through each rxn
    
    if (isempty(reaction_list{i,1})) % If a rxn does not exist
        rxns_to_remove(j) = gene_names(i); % record gene name
        j = j + 1;
    end
    
end  

model = deleteModelGenes(model,...
    rxns_to_remove); % Delete genes with no corresponding reactions

%% GIMME Integration
% Random models for transcriptomics
GIMME_Models_Transcriptomics = Generate_GIMME_Models(Random_Transcriptomics,...
    Medians_Random_Transcriptomics, model);

% Random models for proteomics
GIMME_Models_Proteomics = Generate_GIMME_Models(Random_Proteomics,...
    Medians_Random_Proteomics, model);

 %% essentiality predictions 
%%
 % singleRxnDeletion - deletes a reaction, sees how the model grows, and
% compares to previous growth
%% Essentiality simulations for transcriptomics

% Conduct Single Reaction Deletions on all stored models
lethal_reactions_T = cell(1,150); % each column = different model
for n = 1:150 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_Transcriptomics{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_Transcriptomics{n},'FBA',...
     GIMME_Models_Transcriptomics{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_T(j,n) = GIMME_Models_Transcriptomics{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_T);
writetable(P,'Essential_Rxns_Transcriptomics_Random.csv','Delimiter',',');

clearvars n j i grRatio P

%% Essentiality simulations for proteomics

% Conduct Single Reaction Deletions on all stored models
lethal_reactions_P = cell(1,150); % each column = different model
for n = 1:150 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_Proteomics{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_Proteomics{n},'FBA',...
     GIMME_Models_Proteomics{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_P(j,n) = GIMME_Models_Proteomics{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_P);
writetable(P,'Essential_Rxns_Proteomics_Random.csv','Delimiter',',');

clearvars n j i grRatio P
