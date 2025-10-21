function matcont_data_extractor(filename, cont_param_name, plot_variable_name)
% MATCONT_DATA_EXTRACTOR
% Loads a MATCONT continuation file and extracts the state variables 
% and the continuation parameter into named variables, plotting a 
% bifurcation diagram with stability information.
%
% Usage:
% 1. Ensure the MATCONT .mat file (e.g., 'EQ1.mat') is in your MATLAB path.
% 2. Run the function: matcont_data_extractor('EQ1.mat', 'gam', 'ib');
%
% Inputs:
%   filename              - The path to the MATCONT .mat file (e.g., 'EQ1.mat').
%   cont_param_name       - The name of the parameter varied (e.g., 'gam', 'mu').
%   plot_variable_name    - The name of the state variable to plot 
%                           ('sg', 'sb', 'ib', 'v', or 'phi').

% Check input arguments
if nargin < 3
    fprintf('Error: Three input arguments required: filename, cont_param_name, plot_variable_name.\n');
    return;
end

% Check if the file exists
if ~isfile(filename)
    fprintf('Error: File not found at %s. Please check the path.\n', filename);
    return;
end

% --- 1. Load the Data ---
s = load(filename); 
data_field_name = fieldnames(s);

if isempty(data_field_name)
    fprintf('Error: Could not find data structure in %s.\n', filename);
    return;
end
data = s.(data_field_name{1});

% --- 2. Extract the Core State Array and Main Structure ---
if iscell(data)
    % Find the main results structure within the cell array (usually the last one)
    main_struct = data{end};
else
    main_struct = data;
end

if ~isfield(main_struct, 'x') || ~isfield(main_struct, 'label') || ~isfield(main_struct, 'ftest')
    fprintf('Error: The loaded structure is missing key fields (x, label, or ftest).\n');
    return;
end

X = main_struct.x;
% FT is the greatest real part of the eigenvalue (stability indicator)
FT = main_struct.ftest(1, :); 

% --- 3. Map Indices to Symbolic Variables (Based on ode_system.m) ---
% Index mapping: X(1, :)->sg, X(2, :)->sb, X(3, :)->ib, X(4, :)->v, X(5, :)->phi
state_vars = {'sg', 'sb', 'ib', 'v', 'phi'};
state_data = cell(1, 5);

if size(X, 1) < 5
    fprintf('Error: Data array X has fewer than 5 rows (expected 5 state variables).\n');
    return;
end

for i = 1:5
    state_data{i} = X(i, :);
end

% Assign the chosen plot variable
plot_idx = find(strcmpi(plot_variable_name, state_vars));
if isempty(plot_idx)
    fprintf('Error: Invalid plot variable name "%s". Must be one of: sg, sb, ib, v, phi.\n', plot_variable_name);
    return;
end
Y_plot = state_data{plot_idx};

% Extract Continuation Parameter (assumed to be the 6th row for a 5D system)
if size(X, 1) >= 6
    cont_param = X(6, :);
else
    fprintf('Warning: Continuation parameter not found in the 6th row.\n');
    cont_param = []; 
end

if isempty(cont_param)
    fprintf('Cannot plot without the continuation parameter data.\n');
    return;
end

% --- 4. Plot Bifurcation Diagram with Stability and Labels ---

fprintf('\nSuccessfully loaded data from: %s\n', filename);
fprintf('Total data points on curve: %d\n', length(Y_plot));
fprintf('Plotting %s vs. %s with stability analysis.\n', plot_variable_name, cont_param_name);

figure;
hold on;

% --- Apply Publication Formatting ---
h_axes = gca;
set(h_axes, 'FontName', 'Times New Roman', 'FontSize', 12); % Base font size for axes
set(h_axes, 'Box', 'off', 'TickDir', 'out');
set(gcf, 'Color', 'w'); % Set figure background to white

% Set Title and Labels with larger font
h_title = title(sprintf('Bifurcation Diagram: %s vs. %s', plot_variable_name, cont_param_name));
h_xlabel = xlabel(sprintf('Continuation Parameter: %s', cont_param_name));
h_ylabel = ylabel(sprintf('State Variable: %s', plot_variable_name));
set([h_title, h_xlabel, h_ylabel], 'FontName', 'Times New Roman', 'FontSize', 14);

grid on;

% --- 4a. Stability-based plotting ---
% Stable: FT < 0 -> Solid line (-)
% Unstable: FT > 0 -> Dashed line (--)

% Find indices where stability changes sign (around zero)
% We use a small tolerance (1e-6) to account for numerical zeros/noise
stab_change_indices = find(sign(FT(1:end-1) < -1e-6) ~= sign(FT(2:end) < -1e-6)) + 1;
all_indices = [1, stab_change_indices, length(FT)];

% Iterate through segments
for i = 1:(length(all_indices) - 1)
    start_idx = all_indices(i);
    end_idx = all_indices(i+1);
    
    % Determine stability for the segment (check the sign of the start point)
    is_stable = FT(start_idx) < -1e-6; 
    
    if is_stable
        line_style = '-'; % Solid for stable
    else
        line_style = '--'; % Dashed for unstable
    end
    
    % Plot the segment in blue
    plot(cont_param(start_idx:end_idx), Y_plot(start_idx:end_idx), ...
         'Color', 'b', ...
         'LineStyle', line_style, ...
         'LineWidth', 2.5);
end


% --- 4b. Plot and Label Bifurcation Points ---
bifurcation_points = main_struct.label;

for i = 1:length(bifurcation_points)
    bif_label = bifurcation_points(i).label;
    bif_index = bifurcation_points(i).index;
    
    % Extract coordinates
    bif_x = cont_param(bif_index);
    bif_y = Y_plot(bif_index);
    
    % Plot the bifurcation point as a red asterisk
    plot(bif_x, bif_y, 'r*', 'MarkerSize', 10, 'LineWidth', 1.5);
    
    % Add text label slightly offset
    text(bif_x, bif_y, [' \leftarrow ', bif_label], ...
         'Color', 'k', ...
         'FontSize', 10, ...
         'FontName', 'Times New Roman', ...
         'VerticalAlignment', 'bottom');
end

hold off;

end
