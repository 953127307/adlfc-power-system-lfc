function plot_results()
% Plot the reproduced trajectories. Numeric data match the paper; layout is simplified.
root=fileparts(mfilename('fullpath')); out=fullfile(root,'results');
figdir=fullfile(out,'figures'); if ~isfolder(figdir), mkdir(figdir); end
colors=lines(4); styles={'-','--','-.',':'};
for partcell={'A','B'}
    part=partcell{1}; cases=cell(3,1);
    for ci=1:3, cases{ci}=load(fullfile(out,sprintf('%s_case_%d.mat',part,ci)),'t','x_cell','u_cell','P_f','P_f_d','ACE'); end
    if part=='A', count=6; labels={'Frequency deviation (rad/s)','Tie-line power (p.u.)','Turbine power (p.u.)','Governor output (p.u.)','Control input (p.u.)','ACE (p.u.)'};
    else, count=4; labels={'Frequency deviation (rad/s)','Tie-line power (p.u.)','ACE (p.u.)','Control input (p.u.)'}; end
    for number=1:count
        fig=figure('Visible','off','Color','w','Position',[100,100,900,850]);
        for ci=1:3
            d=cases{ci}; ax=subplot(3,1,ci,'Parent',fig); hold(ax,'on');
            if number==1, values=cell2mat(cellfun(@(v)v(:,2),d.x_cell,'UniformOutput',false)');
            elseif number==2, values=d.P_f;
            elseif part=='A' && number<=4, values=cell2mat(cellfun(@(v)v(:,number),d.x_cell,'UniformOutput',false)');
            elseif (part=='A' && number==5) || (part=='B' && number==4), values=cell2mat(d.u_cell');
            else, values=d.ACE; end
            for ai=1:size(values,2)
                idx=envelope_indices(values(:,ai));
                plot(ax,d.t(idx),values(idx,ai),'Color',colors(ai,:),'LineStyle',styles{ai},'LineWidth',0.9);
            end
            if number==2
                for ai=1:3, yline(ax,d.P_f_d(1,ai),':','Color',[.5,.5,.5],'HandleVisibility','off'); end
                legend(ax,{'P12','P23','P34'},'Location','best');
            else, legend(ax,{'Area 1','Area 2','Area 3','Area 4'},'Location','best'); end
            xlabel(ax,'Time (s)'); ylabel(ax,labels{number}); grid(ax,'on'); xlim(ax,[0,d.t(end)]);
            if part=='A', title(ax,sprintf('Noise case %d',ci));
            else, methods={'ADLFC','ACE-SMC','Integral ACE'}; title(ax,methods{ci}); end
        end
        exportgraphics(fig,fullfile(figdir,sprintf('Fig_%s_%d.png',part,number)),'Resolution',180); close(fig);
    end
end
fig=figure('Visible','off','Color','w','Position',[100,100,900,850]);
for ci=1:3
    d=load(fullfile(out,sprintf('B_case_%d.mat',ci)),'t','E_ACE'); ax=subplot(3,1,ci,'Parent',fig);
    idx=envelope_indices(d.E_ACE); plot(ax,d.t(idx),d.E_ACE(idx),'LineWidth',1); grid(ax,'on');
    methods={'ADLFC','ACE-SMC','Integral ACE'}; title(ax,methods{ci}); xlabel(ax,'Time (s)'); ylabel(ax,'System-wide ACE (p.u.)');
end
exportgraphics(fig,fullfile(figdir,'Fig_System_ACE.png'),'Resolution',180); close(fig);
fprintf('Eleven trajectory figures written to results/figures.\n');
end

function idx=envelope_indices(values)
% Keep original minima/maxima in each block; do not smooth the response.
n=numel(values); block=max(1,ceil(n/4000)); idx=[1;n];
for first=1:block:n
    last=min(first+block-1,n); [~,lo]=min(values(first:last)); [~,hi]=max(values(first:last));
    idx=[idx;first+lo-1;first+hi-1]; %#ok<AGROW>
end
idx=unique(idx);
end
