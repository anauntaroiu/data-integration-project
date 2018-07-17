function gimme_models = Generate_GIMME_Models(data, medians, model)
% Used to run a set of expression data through GIMME algorithim to 
% generate conditional models
%
% INPUT: set of expression data, medians of each set of data, COBRA model
% OUTPUT: Array of conditional models for each set of data

gimme_models = {}; % Generate array to store conditional models

for i = 1:1:width(data) % Loop through sets of data
    % Generate conditional model by calling GIMME
    new_model = GIMME(model, table2array(data(:,i)), table2array(medians(i,1)));
    gimme_models(end+1) = {new_model}; % Save conditional model in array
end

end
