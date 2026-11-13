%[text] This code is used to generate an enzyme constrained Model - ecModel using GECKO.
%[text] 
%[text] To begin, first let's clean the console and the workspace.
clear all
clc
%%
%[text] ## Define directory to modelAdapter class
%[text] Now, we can defined the working directory as well use this live format to store the console output, which can be very useful to follow any error during the execution.
workDir = '~/Documents/Training/04_bioinformatics/10-metabolic-models/exercises/';
cd([workDir 'code'])
%%
%[text] GECKO uses a class adapter to get all the relevant info to reconstruct the ecModel. So, let's load it
% Set default model adapter
adapterLocation = fullfile(workDir,'data','s_cerevisiae', 's_cerevisiaeAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation, 'true');
%%
%[text] ## Load the conventional model to be expanded
%[text] Now, we need a conventional model to expand it and get the ecModel. Please see GECKO v3 manuscript for more details.
model = loadConventionalGEM();
%%
%[text] ## Get ecModel
%[text] Perfect, we can now get an initial ecModel to fill further with enzyme data. So, build the ecModel
ecModel = makeEcModel(model);

% Apply complex data
[ecModel, foundComplex, proposedComplex] = applyComplexData(ecModel); %[output:92e14e5d]

% Get ec numbers from uniprot and KEGG
ecModel = getECfromDatabase(ecModel);

% Use gecko database for determine kcat values
kcatList.fuzzy  = fuzzyKcatMatching(ecModel);

% There is not kcat for all the proteins reported in yeast. Then let's use to 
% strategies: (i) DLKcat approach, (ii) Add an standard kcat value

% DLKcat was performed before. So load the predicted values
kcatList.DLKcat = readDLKcatOutput(ecModel);

kcatList = mergeDLKcatAndFuzzyKcats(kcatList.DLKcat, kcatList.fuzzy);

ecModel  = selectKcatValue(ecModel, kcatList);

% Second strategy will use the values stored in ecModel.ec.kcat, and calculate
% median and mean base on this. Then first let's curate the kcats with some 
% manual data
[ecModel, rxnUpdated, notMatch] = applyCustomKcats(ecModel); %[output:998cf77c]

% Let's make sure all the isoenzymes have a kcat value assign
ecModel = getKcatAcrossIsozymes(ecModel);

% Get the standard kcat
[ecModel, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(ecModel);
disp(['Standard values assigned. MW: ' num2str(standardMW) '. Kcat: ' num2str(standardKcat)]) %[output:9d6112b9]

% Finally, apply the kcat constraints, and set the protein pool size
ecModel = applyKcatConstraints(ecModel);
ecModel = setProtPoolSize(ecModel);
%%
%[text] ## Testing the model
%[text] We have the model now. But first let's check if the model can grow at experimental growth rate. Otherwise, we need to automatically relax the kcat values.
% Set glucose to experimental values
ecModel = setParam(ecModel, 'lb', 'r_1714', -16);
sol = solveLP(ecModel);

disp(['Growth rate before sensitivity tunning: ' num2str(abs(sol.f))]) %[output:13ec2980]
%%
%[text] ## Flexibilize some limiting kcats
%[text] Since model cannot reach the experimental growth rate in WT strain, let's do a sensitivity tunning.
% Set to minimal medium before adjust 
ecModel = changeMediaBatchYeast(ecModel, 'D-glucose', 'MIN');

[ecModel, tunedKcats] = sensitivityTuning(ecModel, '', '', 2); %[output:465443dd]

% Save tunning results
writetable(struct2table(tunedKcats), ...
    fullfile([workDir 'results/s_cerevisiae/'], 'tuned_kcats.txt'), ...
    "Delimiter", '\t');
%%
% Test again
ecModel = setParam(ecModel, 'lb', 'r_1714', -16);

sol = solveLP(ecModel);
disp(['Growth rate after sensitivity tunning: ' num2str(abs(sol.f))]) %[output:56d43d51]

% Save simulation results
printFluxes(ecModel, sol.x, true, 10^-8, fullfile([workDir 'results/s_cerevisiae/'], 'exc_flux_ecYeast.txt')); %[output:93d41d93]

% Remove glucose constraint to save 
ecModel = setParam(ecModel, 'lb', 'r_1714', -1000);

% At this point we have an ecModel that works for further proteomics 
% integration. Let's save this ecModel
saveEcModel(ecModel);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:92e14e5d]
%   data: {"dataType":"text","outputData":{"text":"A total of 208 complex have full match, and 17 proposed.\n","truncated":false}}
%---
%[output:998cf77c]
%   data: {"dataType":"text","outputData":{"text":"WARNING: Protein(s) D4YIP5 were not found in the model.WARNING: Protein(s) Q8ZND6 were not found in the model.","truncated":false}}
%---
%[output:9d6112b9]
%   data: {"dataType":"text","outputData":{"text":"Standard values assigned. MW: 44898.5. Kcat: 12.3\n","truncated":false}}
%---
%[output:13ec2980]
%   data: {"dataType":"text","outputData":{"text":"Growth rate before sensitivity tunning: 0.35839\n","truncated":false}}
%---
%[output:465443dd]
%   data: {"dataType":"text","outputData":{"text":"Iteration 1: Growth: 0.36109\nIteration 2: Growth: 0.41359\nIteration 3: Growth: 0.44601\n","truncated":false}}
%---
%[output:56d43d51]
%   data: {"dataType":"text","outputData":{"text":"Growth rate after sensitivity tunning: 0.43224\n","truncated":false}}
%---
%[output:93d41d93]
%   data: {"dataType":"text","outputData":{"text":"File successfully saved.\n","truncated":false}}
%---
