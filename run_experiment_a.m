function run_experiment_a()
% Section IV-A: three noise cases, 600 s, forward Euler at 0.001 s.
root=fileparts(mfilename('fullpath'));
oldpath=path; cleanup=onCleanup(@()path(oldpath)); %#ok<NASGU>
addpath(fullfile(root,'src'));
f=load(fullfile(root,'config','frozen_settings.mat'),'settings_A');
cfg=f.settings_A; pam=cfg.pam;
d=load(fullfile(root,'data','controller_design.mat'),'w'); w=d.w;
out=fullfile(root,'results'); if ~isfolder(out), mkdir(out); end
for ci=1:3
    bar_d=cfg.noise_bounds(ci);
    fprintf('A case %d: duration=%g s, dt=%g s, noise bound=%g\n',ci,cfg.duration,cfg.dt,bar_d);
    [t,x_cell,u_cell,P_f,P_f_d,~,~,~,~,~,governor_log]= ...
        A_mfile_sim(pam,cfg.duration,cfg.dt,bar_d,'verify',w,[]);
    gamma=cell2mat(cellfun(@(x)x(:,2),x_cell,'UniformOutput',false)');
    ACE=(P_f-P_f_d)*pam.B'+gamma.*pam.B_f';
    metrics=A_S5_metrics(t,x_cell,P_f,P_f_d,pam,cfg.dt);
    stats=summarize_run(t,x_cell,u_cell,ACE,governor_log,cfg.dt);
    save(fullfile(out,sprintf('A_case_%d.mat',ci)), ...
        't','x_cell','u_cell','P_f','P_f_d','ACE','pam','w','metrics','governor_log','stats','-v7.3');
    clear t x_cell u_cell P_f P_f_d ACE gamma governor_log metrics stats
end
end
