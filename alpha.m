clear; close all; clc;

dataDir = './data';

targetAlphas = [0.05, 0.5, 1.0];

canon = cell(numel(targetAlphas),1);
for i = 1:numel(targetAlphas)
    a = targetAlphas(i);
    finP = fullfile(dataDir, sprintf('fig_alpha_%0.6f.dat', a));
    t0P  = fullfile(dataDir, sprintf('fig_alpha_%0.6f_t0.dat', a));
    canon{i} = struct('a', a, 'fin', finP, 't0', t0P, ...
                      'ok', (exist(finP,'file')==2) && (exist(t0P,'file')==2));
end

okFlags = cellfun(@(s) s.ok, canon);
if all(okFlags)
    pairs = canon;
else
    finFiles = dir(fullfile(dataDir, 'fig_alpha_*.dat'));
    tmpPairs = struct('a',{},'fin',{},'t0',{},'ok',{});
    for k = 1:numel(finFiles)
        name = finFiles(k).name;
        tok = regexp(name, '^fig_alpha_([0-9]*\.?[0-9]+)\.dat$', 'tokens','once');
        if isempty(tok), continue; end
        aval = str2double(tok{1});
        finP = fullfile(finFiles(k).folder, name);
        t0P  = fullfile(finFiles(k).folder, sprintf('fig_alpha_%0.6f_t0.dat', aval));
        if exist(t0P,'file')==2
            tmpPairs(end+1) = struct('a', aval, 'fin', finP, 't0', t0P, 'ok', true);
        end
    end

    if isempty(tmpPairs)
        error(['No valid fig_alpha_*.dat + *_t0.dat pairs found in "%s".\n' ...
               'Please run your α-parameter study so that BOTH files exist for each α:\n' ...
               '  fig_alpha_<value>.dat  and  fig_alpha_<value>_t0.dat\n' ...
               '(At least the canonical α = 0.05, 0.5, 1.0).'], dataDir);
    end

    [~, idx] = sort([tmpPairs.a]);
    tmpPairs = tmpPairs(idx);

    sel = [];
    for want = targetAlphas
        hit = find(abs([tmpPairs.a] - want) < 1e-12, 1);
        if ~isempty(hit), sel(end+1) = hit; end 
    end
    rem = setdiff(1:numel(tmpPairs), sel, 'stable');
    sel = unique([sel, rem], 'stable');
    sel = sel(1:min(3, numel(sel)));

    pairs = arrayfun(@(k) tmpPairs(k), sel, 'UniformOutput', false);
end

alphaVals = cellfun(@(p) p.a, pairs);
nPanels   = numel(alphaVals);
assert(nPanels >= 1, 'No α-pairs found to plot.');

fprintf('Reading data files...\n');
Datas = cell(nPanels,1);
for i = 1:nPanels
    p = pairs{i};
    fprintf('  α = %.6f\n', p.a);
    Datas{i}.fin = read_phase_field_dat(p.fin);
    Datas{i}.t0  = read_phase_field_dat(p.t0);

    assert(isequal(size(Datas{i}.fin.Phi), size(Datas{i}.t0.Phi)), ...
        'Grid mismatch for α=%.6f between final and t0.', p.a);
end

fprintf('Generating figure...\n');
figure('Color','w','Position',[100 100 800 200 + 300*nPanels], 'Name','Alpha-Effect');

threshold = 0.5;
for i = 1:nPanels
    subplot(nPanels, 1, i);
    draw_panel_like_paper(Datas{i}.fin.X, Datas{i}.fin.Y, Datas{i}.t0.Phi, Datas{i}.fin.Phi, threshold);

    panelLetter = char('a' + (i-1));
    if i <= 3
        ttl = sprintf('(%s) \\alpha = %g', panelLetter, alphaVals(i));
    else
        ttl = sprintf('\\alpha = %g', alphaVals(i));
    end
    title(ttl, 'FontSize', 12);

    if i == nPanels
        xlabel('x','FontSize',11);
    end
    ylabel('y','FontSize',11);
    set(gca,'FontSize',10);
end

sgtitle('Effect of \alpha', ...
    'FontSize',14,'FontWeight','bold');

fprintf('Visualization generated successfully\n');

function S = read_phase_field_dat(filepath)
    fprintf('    Reading: %s\n', filepath);
    A = readmatrix(filepath, 'FileType','text', 'CommentStyle','#');
    if isempty(A) || size(A,2) < 3
        % Fallback parser
        fid = fopen(filepath,'r');
        if fid == -1, error('Cannot open %s', filepath); end
        C = textscan(fid, '%f %f %f %*[^\n]', 'CommentStyle', '#', 'CollectOutput', true);
        fclose(fid);
        A = C{1};
        if isempty(A) || size(A,2) < 3
            error('Data in %s has fewer than 3 numeric columns.', filepath);
        end
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

function draw_panel_like_paper(X, Y, Phi_t0, Phi_final, threshold)
    mask_initial = (Phi_t0   >= threshold);
    mask_final   = (Phi_final>= threshold);
    mask_growth  = (mask_final & ~mask_initial);

    img = 2 * ones(size(Phi_final));
    img(mask_final)  = 0;
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
