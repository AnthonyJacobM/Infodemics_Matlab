function filtered_array = filter_specific_labels(input_array)
% FILTER_SPECIFIC_LABELS Filters a string array to include only specific
% MATCONT bifurcation labels, allowing for optional trailing spaces.
%
% INPUT:
%   input_array - A MATLAB string array (e.g., ["00", "H ", "LP", "99"])
%
% OUTPUT:
%   filtered_array - A string array containing only the matching labels.

% --- 1. Define the specific target labels for inclusion ---
% H, LP, BT, BP, NS are common bifurcation labels in MATCONT (Hopf, Limit Point, etc.)
target_labels = {'H', 'LP', 'BT', 'BP', 'NS'};

% --- 2. Construct the Regular Expression Pattern ---
% This pattern ensures a match to EXACTLY one of the target labels, 
% followed by optional spaces, and nothing else.

% (H|LP|BT|BP|NS) : Groups the target labels using OR '|'
% \s* : Matches zero or more whitespace characters (optional trailing spaces)
% $               : Matches the end of the string

regex_pattern = ['^(' strjoin(target_labels, '|') ')\s*$']; 
% Example of the final pattern: ^(H|LP|BT|BP|NS)\s*$

fprintf('Regex Pattern Used: %s\n', regex_pattern);

% --- 3. Apply the Regex Filter ---

% The regexp function is used to find matches based on the pattern.
% 'once' stops after the first match attempt.
% The output is a cell array where matches are non-empty.
match_indices = ~cellfun('isempty', regexp(input_array, regex_pattern, 'once'));

% --- 4. Apply the Logical Index to Filter the Array ---
filtered_array = input_array(match_indices);

fprintf('\nFiltered Array (Only specific bifurcation labels): \n');
disp(filtered_array);

end

function filtered_array = filter_letters()
% Original function to filter elements containing ONLY letters.
% (Kept for completeness, though not used in the primary function)
input_array = ["00", "H", "LP", "99", "A1", "Z", "Test", "9B", "C2D"];
pattern = '^[a-zA-Z]+$';
match_indices = ~cellfun('isempty', regexp(input_array, pattern, 'once'));
filtered_array = input_array(match_indices);
end
