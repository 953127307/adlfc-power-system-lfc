function [w, certificate] = compute_adlfc_weights()
% Recollect offline data on the corrected plant and apply Algorithm 1.
% Plant parameters are used by the simulator, not by the weight equations.
source_dir = fileparts(mfilename('fullpath'));
previous_dir = pwd;
restore_dir = onCleanup(@() cd(previous_dir));
old_path = path;
restore_path = onCleanup(@() path(old_path));
addpath(fullfile(source_dir,'src'));
cd(fullfile(source_dir,'src','offline'));
output_dir = fullfile(source_dir,'results');
if ~isfolder(output_dir), mkdir(output_dir); end
run('P_sim.m');
bar_d = bar_d_case(1);
original_kappa = karpi;
run('A_S2_calculate_parameters.m');

certificate = struct;
certificate.model_tag = 'negative_frequency_damping_R2';
certificate.data_bound = bar_d;
certificate.sample_count = k;
certificate.kappa = karpi;
certificate.original_kappa = original_kappa;
certificate.requested_feedback_gain = pam.feedback_gain;
certificate.rank = zeros(4,4);
certificate.lambda = lambda;
certificate.robust_margin = zeros(4,4);
certificate.Q_residual = zeros(4,4);
certificate.p_residual = zeros(4,4);
certificate.data_identity_residual = zeros(4,4);
for ai=1:4
    for sj=1:4
        gg=G{ai,sj};
        certificate.rank(ai,sj)=rank(gg./vecnorm(gg,2,2));
        certificate.robust_margin(ai,sj)=Z{ai,sj}*p{ai,sj}+sqrt(k)*bar_d*norm(p{ai,sj});
        certificate.Q_residual(ai,sj)=norm(gg*Q{ai,sj}-Y{ai,sj},'fro');
        certificate.p_residual(ai,sj)=norm(Psi{ai,sj}*p{ai,sj}-C{ai,sj});
        if sj<4
            next_signal=X{ai,sj+1};
        else
            next_signal=U{ai};
        end
        % Independent identity check; does not enter the synthesis.
        eta_check=eta(ai,sj);
        nuisance=(Z{ai,sj}-D{ai,sj}-eta_check*next_signal)/Psi{ai,sj};
        certificate.data_identity_residual(ai,sj)=norm(Z{ai,sj}-D{ai,sj}-eta_check*next_signal-nuisance*Psi{ai,sj});
        assert(certificate.rank(ai,sj)==size(gg,1),'Data rank condition failed');
        assert(certificate.robust_margin(ai,sj)<0,'Robust decay condition failed');
        assert(certificate.Q_residual(ai,sj)<1e-6,'Data reconstruction residual too large');
    end
end
save(fullfile(output_dir,'controller_regenerated.mat'),'w','certificate','X','U','Z','Psi','D','p','Q');
fprintf('DATA_DRIVEN_SYNTHESIS_CERTIFIED: max margin %.4g, max Q residual %.3e\n',...
    max(certificate.robust_margin(:)),max(certificate.Q_residual(:)));
end
