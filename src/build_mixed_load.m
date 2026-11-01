function [t,P_d,delta_load]=build_mixed_load(config)
% Exact random stream and smooth onsets used for the submitted figures.
dt=config.integration_step_s;
t=(0:round(config.duration_s/dt))'*dt;
P_d=repmat(config.baseline_pu(:)',numel(t),1);
delta_load=zeros(size(P_d)); onsets=config.onset_s_by_area;
active=find(t>=onsets(1)); delta_load(active,1)=config.step_pu;
active=find(t>=onsets(2)); tau=t(active)-onsets(2);
ramp=0.5*(1-cos(pi*min(tau/config.smooth_start_s,1)));
frequencies=config.multisine_frequencies_hz(:)'; amplitudes=config.multisine_amplitudes_pu(:);
delta_load(active,2)=ramp.*(sin(2*pi*tau*frequencies)*amplitudes);
active=find(t>=onsets(3)); tau=t(active)-onsets(3);
ramp=0.5*(1-cos(pi*min(tau/config.smooth_start_s,1)));
stream=RandStream('mt19937ar','Seed',config.seed);
a=exp(-dt/config.ou_correlation_time_s); b=config.ou_stationary_std_pu*sqrt(1-a*a);
ou=zeros(numel(active),1); innovation=randn(stream,numel(active)-1,1);
for k=2:numel(active), ou(k)=a*ou(k-1)+b*innovation(k-1); end
delta_load(active,3)=ramp.*ou;
P_d=P_d+delta_load;
assert(all(P_d(:)>0));
end
