function plot_figure_1( res )

n_runs = length(res.runs);
n_models = length(res.models);

l = zeros(n_runs, n_models);
l_std = zeros(n_runs, n_models);

wald = zeros(n_runs, n_models);
wald_std = zeros(n_runs, n_models);

for i = 1:n_runs
    n_samples = length(res.runs{i});
    
    x = zeros(n_samples, n_models);
    w = zeros(n_samples, n_models);
    for j = 1:n_samples
        for k = 1:n_models
            x(j, k) = lambda_GC(res.runs{i}{j}{k}.wald_stat_null);
            w(j, k) = mean(res.runs{i}{j}{k}.wald_stat_causal);
        end
    end

    l(i, :) = mean(x);
    l_std(i, :) = std(x);
    wald(i, :) = mean(w);
    wald_std(i, :) = std(w);
end

run_names = cell(n_runs, 1);
model_names = cell(n_models, 1);

for i = 1:n_runs
    run_names{i} = res.data_generators{i}.name;
end

for i = 1:n_models
    model_names{i} = res.models{i}.name;
end

%%
font_size = 14;

subplot(2, 1, 1);
hold on

barwitherr(wald_std(1:2, :)/sqrt(n_samples), wald(1:2, :));

ylabel('mean Wald statistic', 'FontSize', font_size);
set(gca, 'XTick', 1:(2), 'XTickLabel', run_names(1:2), 'FontSize', font_size);
%legend(model_names(1:end), 'location', 'BestOutside', 'FontSize', font_size);

hold off;

subplot(2, 1, 2);
hold on

barwitherr(wald_std(3:end, :)/sqrt(n_samples), wald(3:end, :));

ylabel('mean Wald statistic', 'FontSize', font_size);
set(gca, 'XTick', 1:(2), 'XTickLabel', run_names(3:4), 'FontSize', font_size);
%legend(model_names(1:end), 'location', 'BestOutside', 'FontSize', font_size);

hold off;

table_from_lambda_GC(l, l_std/sqrt(n_samples));

end

function table_from_lambda_GC(l, l_se)

fprintf('\\begin{tabular}{ r || c | c | c | c }\n');
fprintf('mean $\\lambda_{GC}$ & pop. strat. & pop. strat. & & \\\\\n');
fprintf('(std. error) & $p = 0.05$ & $p = 0.005$ & $p = 0.05$ & $p = 0.005$ \\\\\n');
fprintf('\\hline\n');
fprintf('\\hline\n');

methods = {'Linear regression', 'Linear reg. with PCs', ...
    'LMM', 'FaST-LMM Select', 'PC-Select' };

for i = 1:5
    fprintf('%s ', methods{i});
    
    for j = 1:4
        fprintf('& %.2f (%.2f)', l(j, i), l_se(j, i));
    end
    fprintf('\\\\\n');
end

fprintf('\\end{tabular}\n');

end
