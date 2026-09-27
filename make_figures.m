function make_figures(results_dir, fig_dir)
% Figures from the CS-SOM run: true vs predicted field over the room (cf. paper
% Fig. 6), relative error map (cf. Fig. 8a), and CS-SOM vs truncated SVD.

if nargin < 1, results_dir = 'results'; end
if nargin < 2, fig_dir = 'figures'; end
if ~isfolder(fig_dir), mkdir(fig_dir); end
R = load(fullfile(results_dir, 'error_on_grid_nm1_25_0.55.mat'));
G = load(fullfile(results_dir, 'Est_mtx_grid_4obj_temp.mat'), 'index');

n = size(R.E_true, 1);
ax = linspace(-4.75, 4.75, n);                 % grid in wavelengths
mask = false(n); mask(G.index) = true;         % inside / on the objects
Et = R.E_true; Ep = R.E_pred;
Et(mask) = NaN; Ep(mask) = NaN;
rel = abs(Ep - Et) ./ abs(Et);

fig = figure('Color', 'w', 'Position', [100 100 1350 390]);
cl = [0 prctile(abs(Et(~mask)), 99)];
panels = {abs(Et), abs(Ep), min(rel, 1)};
titles = {'True field |E|', sprintf('CS-SOM prediction (error %.1f%%)', 100*R.gerror(end)), ...
          'Relative error'};
for p = 1:3
    subplot(1, 3, p);
    h = imagesc(ax, fliplr(ax), panels{p}); set(h, 'AlphaData', ~isnan(panels{p}));
    axis xy image; colorbar; set(gca, 'Color', [0.85 0.85 0.85]);
    if p < 3, caxis(cl); colormap(gca, parula); else, caxis([0 1]); colormap(gca, hot); end
    xlabel('x [\lambda]'); ylabel('y [\lambda]'); title(titles{p}); set(gca, 'FontSize', 12);
end
exportgraphics(fig, fullfile(fig_dir, 'field_prediction.png'), 'Resolution', 150);

fig2 = figure('Color', 'w', 'Position', [100 100 700 420]);
b = bar([mean(R.gerror1) mean(R.gerror)] * 100); b.FaceColor = 'flat';
b.CData = [0.6 0.6 0.6; 0.85 0.33 0.1];
set(gca, 'XTickLabel', {'Truncated SVD', 'CS-SOM'}, 'FontSize', 13);
ylabel('Field prediction error \Delta_G [%]'); grid on
title(sprintf('387 measurements (0.55 x unknowns), 25 dB SNR, %d trials', numel(R.gerror)));
exportgraphics(fig2, fullfile(fig_dir, 'error_comparison.png'), 'Resolution', 150);
fprintf('Figures saved to %s\n', fig_dir);
end
