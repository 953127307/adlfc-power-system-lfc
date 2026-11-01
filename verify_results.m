function verify_results()
% Compare regenerated results to compact samples of the submitted trajectories.
root=fileparts(mfilename('fullpath')); out=fullfile(root,'results');
f=load(fullfile(root,'reference','reference_traces.mat'),'reference'); refs=f.reference;
max_error=0;
for j=1:numel(refs)
    ref=refs(j); d=load(fullfile(out,sprintf('%s_case_%d.mat',ref.part,ref.case_id)), ...
        't','x_cell','u_cell','P_f','ACE','stats');
    values={d.t(ref.indices),cell2mat(d.x_cell'),cell2mat(d.u_cell'),d.P_f,d.ACE};
    expected={ref.t,ref.x,ref.u,ref.P_f,ref.ACE};
    for k=1:numel(values)
        actual=values{k}; if k>1, actual=actual(ref.indices,:); end
        err=abs(actual-expected{k});
        assert(all(err(:)<=1e-9+1e-10*abs(expected{k}(:))), ...
            'Reference mismatch: %s case %d variable %d',ref.part,ref.case_id,k);
        max_error=max(max_error,max(err(:)));
    end
    assert(d.stats.rate_constraint_pass);
end
expected=jsondecode(fileread(fullfile(root,'reference','paper_tables.json')));
names={'rocof','nadir','settling'}; filenames={'Table_III_RoCoF.csv','Table_IV_Nadir.csv','Table_V_ACE_settling.csv'};
for j=1:3
    table=readtable(fullfile(out,filenames{j})); actual=table{:,3:end}; rows=expected.(names{j});
    for row=1:numel(rows)
        assert(table.Area(row)==rows(row).area && table.Case(row)==rows(row).case);
        for col=1:size(actual,2)
            if j==3, rounded=str2double(sprintf('%.2f',actual(row,col)));
            else, rounded=str2double(sprintf('%.2e',actual(row,col))); end
            assert(rounded==rows(row).values(col),'Paper table mismatch: %s row %d col %d',names{j},row,col);
        end
    end
end
report=struct('six_cases_match_reference',true,'paper_tables_match',true, ...
    'governor_rate_checks_pass',true,'max_sample_absolute_error',max_error, ...
    'matlab_version',version,'platform',computer);
fid=fopen(fullfile(out,'verification.json'),'w'); assert(fid>=0);
c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));
fprintf('VERIFIED: six cases, all three paper tables, rate limit. Max sampled error %.3g\n',max_error);
end
