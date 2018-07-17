
%% Initialize COBRA
initCobraToolbox % this only needs to be run once

%% import data for subsetted proteomics data

% Set working directory to where data is saved
cd 'C:\Users\Ana\Documents\MATLAB\Capstone\Matching Data Points\Proteomics'

% Load proteomics data (10-90%)
% 10 Percent
Proteomics_10 = readtable('Proteomics_Data_10_percent.csv','Delimiter',',');
Medians_Proteomics_10 = readtable('Medians_for_Proteomics_10_percent.csv','Delimiter',',');
% 20 Percent
Proteomics_20 = readtable('Proteomics_Data_20_percent.csv','Delimiter',',');
Medians_Proteomics_20 = readtable('Medians_for_Proteomics_20_percent.csv','Delimiter',',');
% 30 Percent
Proteomics_30 = readtable('Proteomics_Data_30_percent.csv','Delimiter',',');
Medians_Proteomics_30 = readtable('Medians_for_Proteomics_30_percent.csv','Delimiter',',');
% 40 Percent
Proteomics_40 = readtable('Proteomics_Data_40_percent.csv','Delimiter',',');
Medians_Proteomics_40 = readtable('Medians_for_Proteomics_40_percent.csv','Delimiter',',');
% 50 Percent
Proteomics_50 = readtable('Proteomics_Data_50_percent.csv','Delimiter',',');
Medians_Proteomics_50 = readtable('Medians_for_Proteomics_50_percent.csv','Delimiter',',');
% 60 Percent
Proteomics_60 = readtable('Proteomics_Data_60_percent.csv','Delimiter',',');
Medians_Proteomics_60 = readtable('Medians_for_Proteomics_60_percent.csv','Delimiter',',');
% 70 Percent
Proteomics_70 = readtable('Proteomics_Data_70_percent.csv','Delimiter',',');
Medians_Proteomics_70 = readtable('Medians_for_Proteomics_70_percent.csv','Delimiter',',');
% 80 Percent
Proteomics_80 = readtable('Proteomics_Data_80_percent.csv','Delimiter',',');
Medians_Proteomics_80 = readtable('Medians_for_Proteomics_80_percent.csv','Delimiter',',');
% 90 Percent
Proteomics_90 = readtable('Proteomics_Data_90_percent.csv','Delimiter',',');
Medians_Proteomics_90 = readtable('Medians_for_Proteomics_90_percent.csv','Delimiter',',');

%% import data for subsetted transcriptomics data

% Set working directory to where data is saved
cd 'C:\Users\Ana\Documents\MATLAB\Capstone\Matching Data Points\Transcriptomics'

% Load proteomics data (10-90%)
% 10 Percent
Transcriptomics_10 = readtable('Matching_Transcriptomics_Data_10_percent.csv','Delimiter',',');
Medians_Transcriptomics_10 = readtable('Medians_for_Matching_Transcriptomics_10_percent.csv','Delimiter',',');
% 20 Percent
Transcriptomics_20 = readtable('Matching_Transcriptomics_Data_20_percent.csv','Delimiter',',');
Medians_Transcriptomics_20 = readtable('Medians_for_Matching_Transcriptomics_20_percent.csv','Delimiter',',');
% 30 Percent
Transcriptomics_30 = readtable('Matching_Transcriptomics_Data_30_percent.csv','Delimiter',',');
Medians_Transcriptomics_30 = readtable('Medians_for_Matching_Transcriptomics_30_percent.csv','Delimiter',',');
% 40 Percent
Transcriptomics_40 = readtable('Matching_Transcriptomics_Data_40_percent.csv','Delimiter',',');
Medians_Transcriptomics_40 = readtable('Medians_for_Matching_Transcriptomics_40_percent.csv','Delimiter',',');
% 50 Percent
Transcriptomics_50 = readtable('Matching_Transcriptomics_Data_50_percent.csv','Delimiter',',');
Medians_Transcriptomics_50 = readtable('Medians_for_Matching_Transcriptomics_50_percent.csv','Delimiter',',');
% 60 Percent
Transcriptomics_60 = readtable('Matching_Transcriptomics_Data_60_percent.csv','Delimiter',',');
Medians_Transcriptomics_60 = readtable('Medians_for_Matching_Transcriptomics_60_percent.csv','Delimiter',',');
% 70 Percent
Transcriptomics_70 = readtable('Matching_Transcriptomics_Data_70_percent.csv','Delimiter',',');
Medians_Transcriptomics_70 = readtable('Medians_for_Matching_Transcriptomics_70_percent.csv','Delimiter',',');
% 80 Percent
Transcriptomics_80 = readtable('Matching_Transcriptomics_Data_80_percent.csv','Delimiter',',');
Medians_Transcriptomics_80 = readtable('Medians_for_Matching_Transcriptomics_80_percent.csv','Delimiter',',');
% 90 Percent
Transcriptomics_90 = readtable('Matching_Transcriptomics_Data_90_percent.csv','Delimiter',',');
Medians_Transcriptomics_90 = readtable('Medians_for_Matching_Transcriptomics_90_percent.csv','Delimiter',',');
% 100 Percent
Transcriptomics_100 = readtable('Matching_Transcriptomics_Data_100_percent.csv','Delimiter',',');
Medians_Transcriptomics_100 = readtable('Medians_for_Matching_Transcriptomics_100_percent.csv','Delimiter',',');

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

%% GIMME Integration - Proteomics

% 10% Proteomics Models
GIMME_Models_proteomics_10 = Generate_GIMME_Models(Proteomics_10,...
    Medians_Proteomics_10, model);

clearvars Proteomics_10 Medians_Proteomics_10

% 10 Percent
% Conduct Single Reaction Deletions on all stored models
lethal_reactions_P_10 = cell(1,75); % each column = different model
for n = 1:75 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_proteomics_10{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_proteomics_10{n},'FBA',...
     GIMME_Models_proteomics_10{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_P_10(j,n) = GIMME_Models_proteomics_10{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_P_10);
writetable(P,'Essential_Rxns_Proteomics_10.csv','Delimiter',',');

clearvars n j i grRatio P lethal_reactions_P_10

% 20% Proteomics Models
GIMME_Models_proteomics_20 = Generate_GIMME_Models(Proteomics_20,...
    Medians_Proteomics_20, model);

clearvars Proteomics_20 Medians_Proteomics_20

% 20 Percent
% Conduct Single Reaction Deletions on all stored models
lethal_reactions_P_20 = cell(1,75); % each column = different model
for n = 1:75 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_proteomics_20{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_proteomics_20{n},'FBA',...
     GIMME_Models_proteomics_20{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_P_20(j,n) = GIMME_Models_proteomics_20{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_P_20);
writetable(P,'Essential_Rxns_Proteomics_20.csv','Delimiter',',');

clearvars n j i grRatio P lethal_reactions_P_20

% 30% Proteomics Models
GIMME_Models_proteomics_30 = Generate_GIMME_Models(Proteomics_30,...
    Medians_Proteomics_30, model);

clearvars Proteomics_30 Medians_Proteomics_30

% 30 Percent
% Conduct Single Reaction Deletions on all stored models
lethal_reactions_P_30 = cell(1,75); % each column = different model
for n = 1:75 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_proteomics_30{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_proteomics_30{n},'FBA',...
     GIMME_Models_proteomics_30{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_P_30(j,n) = GIMME_Models_proteomics_30{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_P_30);
writetable(P,'Essential_Rxns_Proteomics_30.csv','Delimiter',',');

clearvars n j i grRatio P lethal_reactions_P_30

% 40% Proteomics Models
GIMME_Models_proteomics_40 = Generate_GIMME_Models(Proteomics_40,...
    Medians_Proteomics_40, model);

clearvars Proteomics_40 Medians_Proteomics_40

% 40 Percent
% Conduct Single Reaction Deletions on all stored models
lethal_reactions_P_40 = cell(1,75); % each column = different model
for n = 1:75 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_proteomics_40{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_proteomics_40{n},'FBA',...
     GIMME_Models_proteomics_40{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_P_40(j,n) = GIMME_Models_proteomics_40{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_P_40);
writetable(P,'Essential_Rxns_Proteomics_40.csv','Delimiter',',');

clearvars n j i grRatio P lethal_reactions_P_40

% 50% Proteomics Models
GIMME_Models_proteomics_50 = Generate_GIMME_Models(Proteomics_50,...
    Medians_Proteomics_50, model);

clearvars Proteomics_50 Medians_Proteomics_50

% 50 Percent
% Conduct Single Reaction Deletions on all stored models
lethal_reactions_P_50 = cell(1,75); % each column = different model
for n = 1:75 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_proteomics_50{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_proteomics_50{n},'FBA',...
     GIMME_Models_proteomics_50{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_P_50(j,n) = GIMME_Models_proteomics_50{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_P_50);
writetable(P,'Essential_Rxns_Proteomics_50.csv','Delimiter',',');

clearvars n j i grRatio P lethal_reactions_P_50

% 60% Proteomics Models
GIMME_Models_proteomics_60 = Generate_GIMME_Models(Proteomics_60,...
    Medians_Proteomics_60, model);

clearvars Proteomics_60 Medians_Proteomics_60

% 60 Percent
% Conduct Single Reaction Deletions on all stored models
lethal_reactions_P_60 = cell(1,75); % each column = different model
for n = 1:75 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_proteomics_60{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_proteomics_60{n},'FBA',...
     GIMME_Models_proteomics_60{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_P_60(j,n) = GIMME_Models_proteomics_60{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_P_60);
writetable(P,'Essential_Rxns_Proteomics_60.csv','Delimiter',',');

clearvars n j i grRatio P lethal_reactions_60

% 70% Proteomics Models
GIMME_Models_proteomics_70 = Generate_GIMME_Models(Proteomics_70,...
    Medians_Proteomics_70, model);

clearvars Proteomics_70 Medians_Proteomics_70

% 70 Percent
% Conduct Single Reaction Deletions on all stored models
lethal_reactions_P_70 = cell(1,75); % each column = different model
for n = 1:75 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_proteomics_70{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_proteomics_70{n},'FBA',...
     GIMME_Models_proteomics_70{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_P_70(j,n) = GIMME_Models_proteomics_70{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_P_70);
writetable(P,'Essential_Rxns_Proteomics_70.csv','Delimiter',',');

clearvars n j i grRatio P lethal_reactions_70

% 80% Proteomics Models
GIMME_Models_proteomics_80 = Generate_GIMME_Models(Proteomics_80,...
    Medians_Proteomics_80, model);

clearvars Proteomics_80 Medians_Proteomics_80

% 80 Percent
% Conduct Single Reaction Deletions on all stored models
lethal_reactions_P_80 = cell(1,75); % each column = different model
for n = 1:75 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_proteomics_80{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_proteomics_80{n},'FBA',...
     GIMME_Models_proteomics_80{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_P_80(j,n) = GIMME_Models_proteomics_80{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_P_80);
writetable(P,'Essential_Rxns_Proteomics_80.csv','Delimiter',',');

clearvars n j i grRatio P lethal_reactions_P_80

% 90% Proteomics Models
GIMME_Models_proteomics_90 = Generate_GIMME_Models(Proteomics_90,...
    Medians_Proteomics_90, model);

clearvars Proteomics_90 Medians_Proteomics_90

% 90 Percent
% Conduct Single Reaction Deletions on all stored models
lethal_reactions_P_90 = cell(1,75); % each column = different model
for n = 1:75 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_proteomics_90{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_proteomics_90{n},'FBA',...
     GIMME_Models_proteomics_90{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_P_90(j,n) = GIMME_Models_proteomics_90{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_P_90);
writetable(P,'Essential_Rxns_Proteomics_90.csv','Delimiter',',');

clearvars n j i grRatio P lethal_reactions_90

%% GIMME Integration - Transcriptomics

% 10% Transcriptomics Models
GIMME_Models_transcriptomics_10 = Generate_GIMME_Models(Transcriptomics_10,...
    Medians_Transcriptomics_10, model);

clearvars Transcriptomics_10 Medians_Transcriptomics_10

% 10 Percent
% Conduct Single Reaction Deletions on all stored models
lethal_reactions_T_10 = cell(1,75); % each column = different model
for n = 1:75 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_transcriptomics_10{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_transcriptomics_10{n},'FBA',...
     GIMME_Models_transcriptomics_10{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_T_10(j,n) = GIMME_Models_transcriptomics_10{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_T_10);
writetable(P,'Essential_Rxns_Transcriptomics_10.csv','Delimiter',',');

clearvars n j i grRatio P lethal_reactions_T_10

% 20% Transcriptomics Models
GIMME_Models_transcriptomics_20 = Generate_GIMME_Models(Transcriptomics_20,...
    Medians_Transcriptomics_20, model);

clearvars Transcriptomics_20 Medians_Transcriptomics_20

% 20 Percent
% Conduct Single Reaction Deletions on all stored models
lethal_reactions_T_20 = cell(1,75); % each column = different model
for n = 1:75 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_transcriptomics_20{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_transcriptomics_20{n},'FBA',...
     GIMME_Models_transcriptomics_20{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_T_20(j,n) = GIMME_Models_transcriptomics_20{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_T_20);
writetable(P,'Essential_Rxns_Transcriptomics_20.csv','Delimiter',',');

clearvars n j i grRatio P lethal_reactions_T_20

% 30% Transcriptomics Models
GIMME_Models_transcriptomics_30 = Generate_GIMME_Models(Transcriptomics_30,...
    Medians_Transcriptomics_30, model);

clearvars Transcriptomics_30 Medians_Transcriptomics_30

% 30 Percent
% Conduct Single Reaction Deletions on all stored models
lethal_reactions_T_30 = cell(1,75); % each column = different model
for n = 1:75 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_transcriptomics_30{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_transcriptomics_30{n},'FBA',...
     GIMME_Models_transcriptomics_30{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_T_30(j,n) = GIMME_Models_transcriptomics_30{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_T_30);
writetable(P,'Essential_Rxns_Transcriptomics_30.csv','Delimiter',',');

clearvars n j i grRatio P lethal_reactions_T_30

% 40% Transcriptomics Models
GIMME_Models_transcriptomics_40 = Generate_GIMME_Models(Transcriptomics_40,...
    Medians_Transcriptomics_40, model);

clearvars Transcriptomics_40 Medians_Transcriptomics_40

% 40 Percent
% Conduct Single Reaction Deletions on all stored models
lethal_reactions_T_40 = cell(1,75); % each column = different model
for n = 1:75 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_transcriptomics_40{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_transcriptomics_40{n},'FBA',...
     GIMME_Models_transcriptomics_40{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_T_40(j,n) = GIMME_Models_transcriptomics_40{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_T_40);
writetable(P,'Essential_Rxns_Transcriptomics_40.csv','Delimiter',',');

clearvars n j i grRatio P lethal_reactions_T_40

% 50% Transcriptomics Models
GIMME_Models_transcriptomics_50 = Generate_GIMME_Models(Transcriptomics_50,...
    Medians_Transcriptomics_50, model);

clearvars Transcriptomics_50 Medians_Transcriptomics_50

% 50 Percent
% Conduct Single Reaction Deletions on all stored models
lethal_reactions_T_50 = cell(1,75); % each column = different model
for n = 1:75 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_transcriptomics_50{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_transcriptomics_50{n},'FBA',...
     GIMME_Models_transcriptomics_50{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_T_50(j,n) = GIMME_Models_transcriptomics_50{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_T_50);
writetable(P,'Essential_Rxns_Transcriptomics_50.csv','Delimiter',',');

clearvars n j i grRatio P lethal_reactions_T_50

% 60% Transcriptomics Models
GIMME_Models_transcriptomics_60 = Generate_GIMME_Models(Transcriptomics_60,...
    Medians_Transcriptomics_60, model);

clearvars Transcriptomics_60 Medians_Transcriptomics_60

% 60 Percent
% Conduct Single Reaction Deletions on all stored models
lethal_reactions_T_60 = cell(1,75); % each column = different model
for n = 1:75 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_transcriptomics_60{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_transcriptomics_60{n},'FBA',...
     GIMME_Models_transcriptomics_60{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_T_60(j,n) = GIMME_Models_transcriptomics_60{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_T_60);
writetable(P,'Essential_Rxns_Transcriptomics_60.csv','Delimiter',',');

clearvars n j i grRatio P lethal_reactions_T_60

% 70% Transcriptomics Models
GIMME_Models_transcriptomics_70 = Generate_GIMME_Models(Transcriptomics_70,...
    Medians_Transcriptomics_70, model);

clearvars Transcriptomics_70 Medians_Transcriptomics_70

% 70 Percent
% Conduct Single Reaction Deletions on all stored models
lethal_reactions_T_70 = cell(1,75); % each column = different model
for n = 1:75 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_transcriptomics_70{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_transcriptomics_70{n},'FBA',...
     GIMME_Models_transcriptomics_70{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_T_70(j,n) = GIMME_Models_transcriptomics_70{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_T_70);
writetable(P,'Essential_Rxns_Transcriptomics_70.csv','Delimiter',',');

clearvars n j i grRatio P lethal_reactions_T_70

% 80% Transcriptomics Models
GIMME_Models_transcriptomics_80 = Generate_GIMME_Models(Transcriptomics_80,...
    Medians_Transcriptomics_80, model);

clearvars Transcriptomics_80 Medians_Transcriptomics_80

% 80 Percent
% Conduct Single Reaction Deletions on all stored models
lethal_reactions_T_80 = cell(1,75); % each column = different model
for n = 1:75 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_transcriptomics_80{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_transcriptomics_80{n},'FBA',...
     GIMME_Models_transcriptomics_80{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_T_80(j,n) = GIMME_Models_transcriptomics_80{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_T_80);
writetable(P,'Essential_Rxns_Transcriptomics_80.csv','Delimiter',',');

clearvars n j i grRatio P lethal_reactions_T_80

% 90% Transcriptomics Models
GIMME_Models_transcriptomics_90 = Generate_GIMME_Models(Transcriptomics_90,...
    Medians_Transcriptomics_90, model);

clearvars Transcriptomics_90 Medians_Transcriptomics_90

% 90 Percent
% Conduct Single Reaction Deletions on all stored models
lethal_reactions_T_90 = cell(1,75); % each column = different model
for n = 1:75 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_transcriptomics_90{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_transcriptomics_90{n},'FBA',...
     GIMME_Models_transcriptomics_90{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_T_90(j,n) = GIMME_Models_transcriptomics_90{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_T_90);
writetable(P,'Essential_Rxns_Transcriptomics_90.csv','Delimiter',',');

clearvars n j i grRatio P lethal_reactions_T_90

% 100% Transcriptomics Models
GIMME_Models_transcriptomics_100 = Generate_GIMME_Models(Transcriptomics_100,...
    Medians_Transcriptomics_100, model);

clearvars Transcriptomics_100 Medians_Transcriptomics_100

% 100 Percent
% Conduct Single Reaction Deletions on all stored models
lethal_reactions_T_100 = cell(1,15); % each column = different model
for n = 1:15 % loop through the models
  j = 1; % reset the index for every model
for i = 1:length(GIMME_Models_transcriptomics_100{n}.rxns) % loop through the rxns
    [grRatio,~,~,~,~,~] = singleRxnDeletion(GIMME_Models_transcriptomics_100{n},'FBA',...
     GIMME_Models_transcriptomics_100{n}.rxns(i)); % predict growth rate after deletion
    if (isnan(grRatio))
       grRatio = 0; 
    end
    if (grRatio <= 0.4) % if growth ratio is under 0.4, consider rxn lethal
       lethal_reactions_T_100(j,n) = GIMME_Models_transcriptomics_100{n}.rxns(i);
       j = j + 1;
    end
end
end

P = table(lethal_reactions_T_100);
writetable(P,'Essential_Rxns_Transcriptomics_100.csv','Delimiter',',');

clearvars n j i grRatio P lethal_reactions_T_100
