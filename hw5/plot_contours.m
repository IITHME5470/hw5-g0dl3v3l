function plot_comparison()

    folder = 'output_file';
    plotdir = 'plots';
    if ~exist(plotdir, 'dir')
        mkdir(plotdir);
    end

    procConfigs = [4, 8, 16];
    allConfigs = [1, 4, 8, 16];
    configColors = containers.Map({1, 4, 8, 16}, {'k', 'r', 'b', 'g'});

    % --- Detect time steps from parallel output ---
    files = dir(fullfile(folder, 'T_x_y_*.dat'));
    allSteps = [];

    for k = 1:length(files)
        name = files(k).name;
        parts = regexp(name, 'T_x_y_(\d+)_\d+_par_\d+\.dat', 'tokens');
        if ~isempty(parts)
            step = str2double(parts{1}{1});
            allSteps = [allSteps, step];
        end
    end

    allSteps = unique(allSteps);
    allSteps = sort(allSteps);

    if length(allSteps) < 3
        error('Not enough time steps available.');
    end

    % --- Randomly select time steps from start, middle, and end ---
    N = length(allSteps);
    rng('shuffle');
    start_sample  = allSteps(randi([1, floor(N/3)]));
    middle_sample = allSteps(randi([floor(N/3)+1, floor(2*N/3)]));
    end_sample    = allSteps(randi([floor(2*N/3)+1, N]));

    timeSteps = [start_sample, middle_sample, end_sample];
    fprintf("📌 Randomly chosen time steps: %d, %d, %d\n", timeSteps);

    for t = timeSteps
        fprintf("📊 Processing time step %d\n", t);
        solutions = containers.Map('KeyType', 'double', 'ValueType', 'any');

        % --- Serial File ---
        serial_filename = fullfile(folder, sprintf('T_x_y_%06d.dat', t));
        if exist(serial_filename, 'file')
            serial_data = dlmread(serial_filename);
            [x, y, Tgrid] = reassembleField(serial_data);
            solutions(1) = {x, y, Tgrid};
        else
            fprintf("⚠️ Serial file missing for time step %d\n", t);
        end

        % --- Parallel Files ---
        for p = procConfigs
            pattern = fullfile(folder, sprintf('T_x_y_%06d_*_par_%d.dat', t, p));
            files = dir(pattern);
            if isempty(files)
                fprintf('⚠️ No files for step %d with %d procs.\n', t, p);
                continue;
            end

            fullData = [];
            for f = 1:length(files)
                fileData = dlmread(fullfile(files(f).folder, files(f).name));
                fullData = [fullData; fileData];
            end
            [x, y, Tgrid] = reassembleField(fullData);
            solutions(p) = {x, y, Tgrid};
        end

        % --- Contour Plots (Includes Serial) ---
        keys = solutions.keys;
        for k = 1:length(keys)
            proc = keys{k};
            val = solutions(proc);
            x = val{1}; y = val{2}; Tgrid = val{3};

            figure;
            [X, Y] = meshgrid(y, x);
            contourf(X, Y, Tgrid, 20);
            colorbar;
            title(sprintf('Contour at t = %d, Procs = %d', t, proc));
            xlabel('x'); ylabel('y');
            saveas(gcf, fullfile(plotdir, sprintf('contour_ts%06d_par_%d.png', t, proc)));
            close;
        end

        % --- Line Plot at Centerline for All ---
        figure; hold on;
        for k = 1:length(keys)
            proc = keys{k};
            val = solutions(proc);
            x = val{1}; y = val{2}; Tgrid = val{3};
            center_j = round(length(y)/2);
            plot(x, Tgrid(:, center_j), ...
                'DisplayName', sprintf('Procs = %d', proc), ...
                'Color', configColors(proc), ...
                'LineWidth', 1.5);
        end
        title(sprintf('Center Line at t = %d', t));
        xlabel('x'); ylabel('T(x)');
        legend;
        grid on;
        saveas(gcf, fullfile(plotdir, sprintf('line_ts%06d.png', t)));
        close;
    end
end

function [x_unique, y_unique, Tgrid] = reassembleField(data)
    x_vals = data(:,1);
    y_vals = data(:,2);
    T_vals = data(:,3);

    x_unique = unique(x_vals);
    y_unique = unique(y_vals);
    nx = length(x_unique);
    ny = length(y_unique);
    Tgrid = zeros(nx, ny);

    [~, x_idx] = ismember(x_vals, x_unique);
    [~, y_idx] = ismember(y_vals, y_unique);

    for i = 1:length(T_vals)
        Tgrid(x_idx(i), y_idx(i)) = T_vals(i);
    end
end
