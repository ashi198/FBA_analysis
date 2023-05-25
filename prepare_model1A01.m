model_1=importdata('model_1A01.mat');
model_old=convertOldStyleModel(model_1);
%% 
model_test=model_old; 
model_test.genes=cellArray;
cellArray=cellstr(model_old.genes);

%% edit rxnumbers 
% Create a sample cell array with elements in double quotes
rxNumbers= model_old.rxnECNumbers;
% Use cellfun to extract each string into a separate cell
extractedCellArray = cellfun(@(x) cellstr(x), rxNumbers, 'UniformOutput', false);

% Loop through the cell array to extract information from 5x1 strings
for i = 1:numel(extractedCellArray)
    % Check if the element is a 5x1 string
    if numel(extractedCellArray{i}) > 1
        % Concatenate the string elements with commas
        concatenatedString = strjoin(extractedCellArray{i}', ',');
        % Store the concatenated string in a single cell
        extractedCellArray{i} = {concatenatedString};
        if numel(extractedCellArray{i}) == 0
        extractedCellArray{i} = "";
        end 
    end
end
%% extract extractCellArray in a separate excel file and continue this 

filename = 'input_rxns.xlsx';
rxNumbers = readtable(filename);

% Convert numerical data to strings
rx_transformed=table2array(rxNumbers);
model_test.rxnECNumbers=rx_transformed;

%% 
model_test.rules{245}='(( x(1390) & x(1391) & x(1392) & x(1293) & x(1294) & x(1395) & x(1396) & x(1397) ) | ( x(1398) & x(1399) & x(1400) & x(1401) & x(1402) & x(1403) & x(1404) & x(1405) ))';
model_1A01=writeCbModel(model_test, 'model_1A01_gralka.xml', 'format', 'sbml');

Kbase=readCbModel('1A01_DEMETER.mat');