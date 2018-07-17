
%% Initialize COBRA
initCobraToolbox % this only needs to be run once

%% import data for absolute proteomics data

% Set working directory to where data is saved
cd 'C:\Users\Ana\Documents\MATLAB\Capstone\Absolute Proteomics'

% Load data
Absolute = readtable('Absolute Proteomics For GIMME.csv','Delimiter',',');
Medians_Absolute = readtable('Medians For Absolute Proteomics.csv','Delimiter',',');

% Load random data
Random_Absolute = readtable('Random Samples For Absolute Proteomics.csv','Delimiter',',');
Medians_Random_Absolute = readtable('Random Sample Medians For Absolute Proteomics.csv','Delimiter',',');

%% import data for relative proteomics data

% Set working directory to where data is saved
cd 'C:\Users\Ana\Documents\MATLAB\Capstone\Relative Proteomics'

% Load data
Relative = readtable('Relative Proteomics For GIMME.csv','Delimiter',',');
Medians_Relative = readtable('Medians For Relative Proteomics.csv','Delimiter',',');

% Load random data
Random_Relative = readtable('Random Samples For Relative Proteomics.csv','Delimiter',',');
Medians_Random_Relative = readtable('Random Sample Medians For Relative Proteomics.csv','Delimiter',',');

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

% Models for absolute proteomics
GIMME_Models_absolute = Generate_GIMME_Models(Absolute,...
    Medians_Absolute, model);

% Random models for absolute proteomics
GIMME_Models_absolute_random = Generate_GIMME_Models(Random_Absolute,...
    Medians_Random_Absolute, model);

% Models for relative proteomics
GIMME_Models_relative = Generate_GIMME_Models(Relative,...
    Medians_Relative, model);

% Random models for relative proteomics
GIMME_Models_relative_random = Generate_GIMME_Models(Random_Relative,...
    Medians_Random_Relative, model);

 %% essentiality predictions 
%%
 % singleRxnDeletion - deletes a reaction, sees how the model grows, and
% compares to previous growth
%% Essentiality simulations for absolute proteomics

% Conduct Single Reaction Deletions on all stored models
lethal_reactions_A = cell(1,15); % each column = different model
for n = 1:15 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_absolute{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_absolute{n},'FBA',...
     GIMME_Models_absolute{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_A(j,n) = GIMME_Models_absolute{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_A);
writetable(P,'Essential_Rxns_Absolute_Proteomics.csv','Delimiter',',');

clearvars n j i grRatio P

%% Essentiality simulations for absolute proteomics - random

% Conduct Single Reaction Deletions on all stored models
lethal_reactions_A_R = cell(1,150); % each column = different model
for n = 1:150 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_absolute_random{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_absolute_random{n},'FBA',...
     GIMME_Models_absolute_random{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_A_R(j,n) = GIMME_Models_absolute_random{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_A_R);
writetable(P,'Essential_Rxns_Absolute_Proteomics_Random.csv','Delimiter',',');

clearvars n j i grRatio P

%% Essentiality simulations for relative proteomics

% Conduct Single Reaction Deletions on all stored models
lethal_reactions_R = cell(1,15); % each column = different model
for n = 1:15 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_relative{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_relative{n},'FBA',...
     GIMME_Models_relative{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_R(j,n) = GIMME_Models_relative{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_R);
writetable(P,'Essential_Rxns_Relative_Proteomics.csv','Delimiter',',');

clearvars n j i grRatio P

%% Essentiality simulations for relative proteomics - random

% Conduct Single Reaction Deletions on all stored models
lethal_reactions_R_R = cell(1,150); % each column = different model
for n = 1:150 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_relative_random{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_relative_random{n},'FBA',...
     GIMME_Models_relative_random{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_R_R(j,n) = GIMME_Models_relative_random{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_R_R);
writetable(P,'Essential_Rxns_Relative_Proteomics_Random.csv','Delimiter',',');

clearvars n j i grRatio P
