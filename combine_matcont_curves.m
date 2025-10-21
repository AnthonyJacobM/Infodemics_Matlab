function combined_data = combine_matcont_curves(file1, file2, new_filename)
% COMBINE_MATCONT_CURVES Merges two MATCONT continuation structures (e.g., FWD and BKWD runs).
%
% This function is essential for creating a complete bifurcation curve from 
% two runs started from the same initial point, continuing in opposite 
% directions. It concatenates the state vector data, tangent vectors, and 
% bifurcation labels.
%
% INPUTS:
%   file1      - String, filename of the first MATCONT run (e.g., 'EP_FWD.mat').
%   file2      - String, filename of the second MATCONT run (e.g., 'EP_BKWD.mat').
%   new_filename - String, filename to save the combined result (e.g., 'EP_FULL.mat').
%
% OUTPUT:
%   combined_data - The merged structure.

fprintf('Loading curve 1: %s\n', file1);
data1_struct = load(file1);
% Find the single primary structure variable name (e.g., 'EP_EP_Rej_eta_default_FWD')
varName1 = fieldnames(data1_struct);
data1 = data1_struct.(varName1{1});

fprintf('Loading curve 2: %s\n', file2);
data2_struct = load(file2);
varName2 = fieldnames(data2_struct);
data2 = data2_struct.(varName2{1});

% Ensure both structures have the necessary fields
if ~all(isfield(data1, {'x', 'v', 's', 'f'})) || ~all(isfield(data2, {'x', 'v', 's', 'f'}))
    error('One or both MATCONT structures are missing required fields (x, v, s, f).');
end

% --- Data Concatenation ---
% The primary data arrays ('x', 'v', 'f') are concatenated horizontally.
% We typically keep the first point from the 'data1' run, and exclude the 
% first point from the 'data2' run to avoid duplication, as both FWD and BKWD 
% runs share the starting point.

% 1. State vectors (x): [states x points]
% Concatenate all points from data1 with points 2-end from data2.
combined_data.x = [data1.x, data2.x(:, 2:end)]; 

% 2. Tangent vectors (v): [states x points]
combined_data.v = [data1.v, data2.v(:, 2:end)];

% 3. Function/FT values (f, h): [values x points]
% FT (ftest) is crucial for stability and is often the first element of f.
combined_data.f = [data1.f, data2.f(:, 2:end)]; 

% 4. Bifurcation points (s): Array of structures
% The 's' array stores the bifurcation labels. These must be concatenated 
% directly. We can't skip the first label here because the label index 
% might be needed for consistency, but often MATCONT places the starting 
% label ('00') first in both, which is ok.
combined_data.s = [data1.s, data2.s];

% --- Retaining Initialization/Global Data ---
% These fields should be identical and are usually small. Retain them from data1.
combined_data.h = data1.h; % History structure
combined_data.globals = data1.globals;
combined_data.initUsed = data1.initUsed;
if isfield(data1, 'gui')
    combined_data.gui = data1.gui; % GUI settings
end

fprintf('Successfully merged %d points from curve 1 and %d points from curve 2.\n', ...
    size(data1.x, 2), size(data2.x, 2));
fprintf('Total points in combined curve: %d\n', size(combined_data.x, 2));

% --- Save the Merged Data ---
% The variable name saved inside the file must match the new filename for 
% consistency when loading later.
save_variable_name = new_filename(1:end-4); % Remove '.mat' extension
eval(sprintf('%s = combined_data;', save_variable_name));

save(new_filename, save_variable_name);
fprintf('Combined data saved to: %s\n', new_filename);

end
