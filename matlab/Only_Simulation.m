clear
close all

start = 1;
i = start; % Initialize index
file_exists = true; % Flag to check if the file exists
mac = 1;

while file_exists
    if mac
        input_filename = sprintf('/Users/stefanomarchetti/Library/CloudStorage/OneDrive-PolitecnicodiMilano/Documenti/MATLAB/PhD/UMD/Inputs/Input_%d.mat', i);
        output_filename = sprintf('/Users/stefanomarchetti/Library/CloudStorage/OneDrive-PolitecnicodiMilano/Documenti/MATLAB/PhD/UMD/Outputs/Output_%d.mat', i);
    else
        input_filename = sprintf('C:/Users/ste08/OneDrive - Politecnico di Milano/Documenti/MATLAB/PhD/UMD/Inputs/Input_%d.mat', i);
        output_filename = sprintf('C:/Users/ste08/OneDrive - Politecnico di Milano/Documenti/MATLAB/PhD/UMD/Outputs/Output_%d.mat', i);
    end
    if exist(input_filename, 'file')
        load(input_filename);

        out = parsim(in_filtered);

        if exist('file_division','var') == 0
            file_division=1;
        end

        save(output_filename, 'out','file_division','-v7.3');

        clear in_filtered
        clear out

        i = i + 1; % Increment index
    else
        file_exists = false; % Stop if file doesn't exist
    end
end
