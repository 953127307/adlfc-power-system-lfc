function stats=summarize_run(t,x_cell,u_cell,ACE,governor_log,dt)
stats=struct('duration',t(end),'dt',dt);
stats.max_state_rate=zeros(1,4); stats.integration_residual=zeros(1,4);
stats.peak_frequency=zeros(1,4); stats.max_abs_control=zeros(1,4);
for ai=1:4
    rate=diff(x_cell{ai}(:,4))/dt;
    stats.max_state_rate(ai)=max(abs(rate));
    stats.integration_residual(ai)=max(abs(rate-governor_log.actual_rate(1:end-1,ai)));
    stats.peak_frequency(ai)=max(abs(x_cell{ai}(:,2)));
    stats.max_abs_control(ai)=max(abs(u_cell{ai}));
end
stats.peak_ACE=max(abs(ACE),[],1);
stats.max_applied_rate=max(abs(governor_log.actual_rate(1:end-1,:)),[],1);
stats.finite=all(cellfun(@(v)all(isfinite(v(:))),x_cell)) && ...
    all(cellfun(@(v)all(isfinite(v(:))),u_cell));
stats.rate_constraint_pass=stats.finite && max(stats.max_state_rate)<=2+1e-8 && ...
    max(stats.max_applied_rate)<=2 && max(stats.integration_residual)<1e-8;
assert(stats.rate_constraint_pass,'Governor rate or integration consistency failed.');
end
