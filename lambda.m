clear; close all; clc;

dataDir = './data';

files = struct();
files.f30     = fullfile(dataDir, 'fig_lambda0_30.dat');
files.f30_t0  = fullfile(dataDir, 'fig_lambda0_30_t0.dat');
files.f100    = fullfile(dataDir, 'fig_lambda0_100.dat');
files.f100_t0 = fullfile(dataDir, 'fig_lambda0_100_t0.dat');

have_fig = all(cellfun(@(p) exist(p,'file')==2, struct2cell(files)));

if ~have_fig
    iterFiles = dir(fullfile(dataDir, 'iter_*.dat'));
    if isempty(iterFiles)
        error(['No fig_* nor iter_*.dat files found in "%s".\n' ...
               'Run your C++ parameter study (for fig_*) or run solve() to produce iter_*.dat.'], dataDir);
    end

    itNums = nan(numel(iterFiles),1);
    for k = 1:numel(iterFiles)
        tok = regexp(iterFiles(k).name, 'iter_(\d+)\.dat$', 'tokens','once');
        if ~isempty(tok), itNums(k) = str2double(tok{1}); end
    end
    [~, order] = sort(itNums);
    iterFiles = iterFiles(order(~isnan(itNums)));
    if numel(iterFiles) < 2
        error('Need at least two iter_*.dat files to show initial/final panels.');
    end
    files_fallback.t0 = fullfile(iterFiles(1).folder, iterFiles(1).name);
    files_fallback.fin = fullfile(iterFiles(end).folder, iterFiles(end).name);
end


fprintf('Reading data...\n');
if have_fig
    D30     = read_phase_field_dat(files.f30);
    D30_t0  = read_phase_field_dat(files.f30_t0);
    D100    = read_phase_field_dat(files.f100);
    D100_t0 = read_phase_field_dat(files.f100_t0);


    assert(isequal(size(D30.Phi), size(D30_t0.Phi)),  'Grid size mismatch for λ0=30.');
    assert(isequal(size(D100.Phi), size(D100_t0.Phi)),'Grid size mismatch for λ0=100.');
else
    Dfin = read_phase_field_dat(files_fallback.fin);
    Dt0  = read_phase_field_dat(files_fallback.t0);
    assert(isequal(size(Dfin.Phi), size(Dt0.Phi)), 'Grid size mismatch between earliest and latest iter files.');
end


fprintf('Generating figure...\n');
figure('Color','w','Position',[100 100 800 650],'Name',' Lambda-Effect');
threshold = 0.5;

if have_fig
    subplot(2,1,1);
    draw_panel(D30.X, D30.Y, D30_t0.Phi, D30.Phi, threshold);
    title('(a) \lambda_0 = 30', 'FontSize', 12);
    ylabel('y', 'FontSize', 11); set(gca,'FontSize',10);

    subplot(2,1,2);
    draw_panel(D100.X, D100.Y, D100_t0.Phi, D100.Phi, threshold);
    title('(b) \lambda_0 = 100', 'FontSize', 12);
    xlabel('x','FontSize',11); ylabel('y','FontSize',11); set(gca,'FontSize',10);

    sgtitle('Effect of \lambda_0', ...
        'FontSize',14,'FontWeight','bold');

else
    subplot(2,1,1);
    draw_panel(Dfin.X, Dfin.Y, Dt0.Phi, Dfin.Phi, threshold);
    ttlA = sprintf('(a) Earliest vs Latest: %s → %s', ...
        strip_filename(files_fallback.t0), strip_filename(files_fallback.fin));
    title(ttlA,'FontSize',12);
    ylabel('y','FontSize',11); set(gca,'FontSize',10);

    subplot(2,1,2);
    draw_panel(Dfin.X, Dfin.Y, Dt0.Phi, Dfin.Phi, threshold);
    title('(b) Fallback panel (same run)', 'FontSize',12);
    xlabel('x','FontSize',11); ylabel('y','FontSize',11); set(gca,'FontSize',10);

    sgtitle('Phase-field evolution (fallback to iter_*.dat)', ...
        'FontSize',14,'FontWeight','bold');

end

fprintf('Visualization generated scuuessfully\n');

%% ----------------------------- Helpers --------------------------------
function S = read_phase_field_dat(filepath)
    fprintf('  Reading: %s\n', filepath);
    opts = detectImportOptions(filepath, 'FileType','text', 'CommentStyle','#');
    A = readmatrix(filepath, 'FileType','text', 'CommentStyle','#');
    if size(A,2) < 3
        error('Data in %s has fewer than 3 columns.', filepath);
    end
    x   = A(:,1);
    y   = A(:,2);
    phi = A(:,3);

    [xu, ~, ix] = unique(x, 'stable');
    [yu, ~, iy] = unique(y, 'stable');

    S.Nx  = numel(xu);
    S.Ny  = numel(yu);
    S.X   = xu(:)';
    S.Y   = yu(:)';
    S.Phi = accumarray([iy ix], phi, [S.Ny S.Nx]);
end

function draw_panel(X, Y, Phi_t0, Phi_final, threshold)
    mask_initial = (Phi_t0   >= threshold);
    mask_final   = (Phi_final>= threshold);
    mask_growth  = (mask_final & ~mask_initial);

    img = 2 * ones(size(Phi_final));
    img(mask_final) = 0;
    img(mask_growth) = 1;

    imagesc(X, Y, img);
    set(gca,'YDir','normal'); axis tight; box on;

    paper_cmap = [0 0 0; 0 1 0; 1 1 1];
    colormap(gca, paper_cmap);
    caxis([0 2]);

    hold on;
    contour(X, Y, Phi_final, [threshold threshold], 'Color',[0.1 0.1 0.1], 'LineWidth',0.5);
    hold off;
end

function s = strip_filename(p)
    [~, s, ext] = fileparts(p); s = [s ext];
end
