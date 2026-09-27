function run_all
% Runs the full pipeline of "Compressive sensing approaches for the prediction of
% scattered electromagnetic fields" (C. Bhat, K. Sastry and U. K. Khankhoje,
% JOSA A 2020) for the paper's 4-object room, then makes the figures.
% Takes about 10 minutes on a laptop.
%
% Requires CVX (http://cvxr.com/cvx) on the MATLAB path for the CS-SOM step.
% All intermediate .mat files are written to ./results, figures to ./figures.

root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root, 'code'));
out = fullfile(root, 'results');
if ~exist(out, 'dir'), mkdir(out); end
here = pwd; cleanup = onCleanup(@() cd(here));
cd(out);

steps = {
    'Closed_oneObj1',                        'Forward solver (boundary integral, lambda/40): exact surface fields'
    'true_2d_field',                         'True field on the 10 lambda x 10 lambda grid'
    'generate_truetangfields_all_surfaces',  'True tangential fields on the bounding boxes'
    'tang_interpolate',                      'True tangential fields resampled at lambda/5'
    'generate_matrices',                     'Random measurement sets, noisy data and system matrices'
    'generate_state_eqn',                    'State matrix (extinction theorem)'
    'Prediction_matrix_one_time',            'Prediction matrix for the grid (Huygens'' principle)'
    'Inverse_solver_25dB_055',               'CS-SOM and truncated-SVD field prediction'
    };
for s = 1:size(steps, 1)
    fprintf('\n=== Step %d/%d: %s\n', s, size(steps, 1), steps{s, 2});
    run_script(steps{s, 1});
    close all
end
make_figures(out, fullfile(root, 'figures'));
end

function run_script(name)
% Runs a research script in this function's own workspace (the scripts start
% with "clear all", which would otherwise wipe the caller's variables).
eval(name);
end
