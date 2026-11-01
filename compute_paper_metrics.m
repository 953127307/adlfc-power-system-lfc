function compute_paper_metrics()
% Export the three manuscript tables and Fig. 13 peak comparisons.
root=fileparts(mfilename('fullpath')); out=fullfile(root,'results');
R=zeros(12,6); N=zeros(12,5); S=zeros(12,6);
for ci=1:3
    d=load(fullfile(out,sprintf('A_case_%d.mat',ci)),'metrics'); m=d.metrics;
    for ai=1:4
        row=3*(ai-1)+ci;
        R(row,:)=[ai,ci,m.RoCoF(ai,:)];
        N(row,:)=[ai,ci,m.Nadir(ai,1:3)];
        S(row,:)=[ai,ci,m.StepSettling(ai,:)];
    end
end
writetable(array2table(R,'VariableNames',{'Area','Case','t60','t200','t360','t500'}),fullfile(out,'Table_III_RoCoF.csv'));
writetable(array2table(N,'VariableNames',{'Area','Case','t60','t200','t360'}),fullfile(out,'Table_IV_Nadir.csv'));
writetable(array2table(S,'VariableNames',{'Area','Case','t60','t200','t360','t500'}),fullfile(out,'Table_V_ACE_settling.csv'));
peaks=zeros(3,3); starts=[60,200,360];
for ci=1:3
    d=load(fullfile(out,sprintf('B_case_%d.mat',ci)),'t','E_ACE');
    for j=1:3
        mask=d.t>=starts(j) & d.t<=starts(j)+30;
        peaks(ci,j)=max(d.E_ACE(mask));
    end
end
reduction_percent=100*(1-peaks(1,:)./peaks(3,:));
writetable(array2table([starts',peaks',reduction_percent'], ...
    'VariableNames',{'Event_s','ADLFC','ACE_SMC','Integral_ACE','Reduction_percent'}), ...
    fullfile(out,'System_ACE_peaks.csv'));
report=struct('RoCoF',R,'Nadir',N,'Settling',S,'system_ACE_peak',peaks, ...
    'reduction_percent',reduction_percent);
fid=fopen(fullfile(out,'paper_metrics.json'),'w'); assert(fid>=0);
c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));
fprintf('Paper tables written. System ACE peak reduction: %s %%\n',mat2str(reduction_percent,6));
end
