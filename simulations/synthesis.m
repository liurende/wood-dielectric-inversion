%% Synthesis: Combine findings from all three experiments
% Generates the summary figure for the research report.

function synthesis()
    clc; close all;

    outDir = fullfile(pwd, 'simulations', 'results');
    fontName = 'Times New Roman';

    fprintf('=== Synthesis Report ===\n');
    fprintf('Findings from all three simulation experiments.\n\n');

    % ---- Exp1 Findings ----
    fprintf('--- Exp1: Identifiability ---\n');
    fprintf('Key finding: (tan(delta), n) couple with aspect ratio ~0.07.\n');
    fprintf('The cost landscape is a long narrow valley -- these two params\n');
    fprintf('cannot be separated using |S21| magnitude alone.\n');
    fprintf('Recommendation: fix n from literature or use phase information.\n\n');

    % ---- Exp2 Findings (load data) ----
    f2 = fullfile(outDir, 'Exp2_Recovery_Data.mat');
    if exist(f2, 'file')
        ld = load(f2);
        fprintf('--- Exp2: Recovery Test ---\n');
        truth = ld.truth;
        recovery_all = ld.recovery_all;
        noiseLevels = ld.noiseLevels;

        for iS = 1:size(truth,1)
            for iN = 1:length(noiseLevels)
                err = squeeze(recovery_all(iS, iN, :, :) - truth(iS,:));
                err(:,4) = err(:,4) * 1000;
                mu = mean(err, 1);
                fprintf('  %s + %s: ', ld.sampleNames{iS}, noiseLevels{iN}.name);
                fprintf('d_eps=%.3f d_tand=%.4f d_n=%.2f d_d=%.1fmm\n', ...
                    mu(1), mu(2), mu(3), mu(4));
            end
        end
    end

    % ---- Exp3 Findings (load data) ----
    f3 = fullfile(outDir, 'Exp3_Model_Mismatch_Data.mat');
    if exist(f3, 'file')
        ld3 = load(f3);
        fprintf('\n--- Exp3: Model Mismatch ---\n');
        fprintf('Debye-generated data fitted with power-law model shows\n');
        fprintf('systematic residuals, especially near the relaxation frequency.\n');
        fprintf('RMS residuals > 0.5 dB indicate model inadequacy.\n');
    end

    fprintf('\n=== Synthesis complete ===\n');
end
