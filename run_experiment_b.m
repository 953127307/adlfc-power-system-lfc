function run_experiment_b()
% Section IV-B: the same mixed-load trajectory for all three controllers.
root=fileparts(mfilename('fullpath'));
oldpath=path; cleanup=onCleanup(@()path(oldpath)); %#ok<NASGU>
addpath(fullfile(root,'src'));
f=load(fullfile(root,'config','frozen_settings.mat'),'settings_B'); cfg=f.settings_B;
d=load(fullfile(root,'data','controller_design.mat'),'w'); w=d.w;
config=jsondecode(fileread(fullfile(root,'config','mixed_load.json')));
[t,P_d,delta_load]=build_mixed_load(config);
reference=load(fullfile(root,'data','mixed_load_profile.mat'),'t','P_d','delta_load');
assert(isequal(t,reference.t) && isequal(P_d,reference.P_d) && ...
    isequal(delta_load,reference.delta_load),'Mixed-load generation changed.');
clear reference
out=fullfile(root,'results'); if ~isfolder(out), mkdir(out); end
save(fullfile(out,'mixed_load_profile.mat'),'t','P_d','delta_load','config','-v7');
for ci=1:3
    pam=cfg.pam; pam.ctrl_method=ci; pam.P_d_trajectory=P_d;
    fprintf('B case %d: duration=%g s, dt=%g s\n',ci,cfg.duration,cfg.dt);
    [tt,x_cell,u_cell,P_f,P_f_d,governor_log]= ...
        B_mfile_sim(pam,cfg.duration,cfg.dt,cfg.noise_bound,'verify',w,[]);
    assert(isequal(t,tt));
    gamma=cell2mat(cellfun(@(x)x(:,2),x_cell,'UniformOutput',false)');
    ACE=(P_f-P_f_d)*pam.B'+gamma.*pam.B_f'; E_ACE=max(abs(ACE),[],2);
    stats=summarize_run(t,x_cell,u_cell,ACE,governor_log,cfg.dt);
    save(fullfile(out,sprintf('B_case_%d.mat',ci)), ...
        't','x_cell','u_cell','P_f','P_f_d','ACE','E_ACE','pam','w','governor_log','stats','-v7.3');
    clear tt x_cell u_cell P_f P_f_d ACE E_ACE gamma governor_log stats
end
end
