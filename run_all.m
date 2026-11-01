function run_all()
% Reproduce both experiments, tables, regression checks, and figures.
run_experiment_a();
run_experiment_b();
compute_paper_metrics();
verify_results();
plot_results();
fprintf('\nAll experiments, regression checks, and figures completed.\n');
end
