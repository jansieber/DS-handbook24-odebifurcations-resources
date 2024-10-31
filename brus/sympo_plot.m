%% plot one-parameter results from symmetric branch
%% load runs, define variables
clear
ty=@(r,s)strcmp(r{:,'TYPE'},s);
alw={'linewidth',2};
lw={'linewidth',3};
txt={'FontSize',20,'FontName','courier','FontWeight','bold'};
ltx={'Interpreter','latex'};
srun=@(i)sprintf('run%d',i);
for i=4:-1:1
    run{i}=coco_bd_table(srun(i));
    specsel{i}=ty(run{i},'UST')|ty(run{i},'FP')|ty(run{i},'BP');
end
names={'$P_1$: $\Pi\sim(1\,2\,4\,3)_4$',...
    '$P_4$: $\Pi\sim(2\,3)_2(1\,4)_2$',...
    '$P_5$: $\Pi\sim(2\,3)_1(1\,4)_2$ ($\sim P_2$)',...
    '$P_6$: $\Pi\sim(2\,3)_1$'};
shortnames=cellfun(@(n){n(1:5)},names);
proflab={{'BP',1},{'FP',1},{'FP',5},{'EP',2}};
clri=[1,2,4,5];
nruns=length(run);
%% plot of all runs
figure(2);clf;tl=tiledlayout(nruns,2);
clr=lines();
comp_order={'ColorOrder',clr};
tl.TileSpacing='compact';
nexttile([nruns,1]);ax1=gca;hold(ax1,'on');
for i=1:nruns
    pl(i)=plot(ax1,run{i}.lambda,run{i}.amplitude,lw{:},'DisplayName',names{i},...
        'color',clr(clri(i),:));
    uz=plot(ax1,run{i}.lambda(specsel{i}),run{i}.amplitude(specsel{i}),'ks',...
        'markerfacecolor','k','MarkerSize',8);
end
uz.DisplayName='UST/BP/FP';
xlim(ax1,[0,26])
legend(ax1,[pl,uz],ltx{:},txt{:},'Location','North');
xlabel(ax1,'$\epsilon\times10^{3}$',ltx{:},txt{:})
ylabel(ax1,'amplitude',ltx{:},txt{:})
set(ax1,txt{:},alw{:},'box','on');
ax1.LabelFontSizeMultiplier=1.2;
for i=1:nruns
    pl_labs{i}=coco_bd_labs(srun(i),proflab{i}{1});
    tprof{i}=po_read_solution(srun(i),pl_labs{i}(proflab{i}{2}));
    nexttile(2*i);ax(i)=gca;hold(ax(i),'on');
    t=tprof{i}.tbp/tprof{i}.T;
    tpl{i}=plot(ax(i),t,tprof{i}.xbp(:,1:2:end),lw{:});
    for k=1:4
        set(tpl{i}(k),'DisplayName',sprintf('$x_%d$',k));
    end
    pl_points{i}=run{i}(ty(run{i},proflab{i}{1}),{'lambda','amplitude'});
    pl_points{i}=pl_points{i}(proflab{i}{2},:);
    marker={'o','MarkerEdgeColor',clr(clri(i),:),...
        'MarkerFaceColor',clr(clri(i),:),'MarkerSize',10};
    plot(ax1,pl_points{i}.lambda,pl_points{i}.amplitude,marker{:},'HandleVisibility','off')
    plnamelab{i}=plot(ax(i),NaN,NaN,marker{:},'DisplayName',shortnames{i});
    legend(ax(i),[plnamelab{i};tpl{i}],txt{:},ltx{:},'Location','EastOutside');
    set(ax(i),txt{:},alw{:},'box','on',comp_order{:});
    ax(i).XTick=linspace(0,1,5);
    if i<nruns
        ax(i).XTickLabel={};
    end
    ax(i).YTickLabel={};
    ax(i).YLimitMethod='padded';
    ax(i).XLimitMethod='padded';
end
set(tpl{3}(3),'LineStyle','--');
set(tpl{4}(3),'LineStyle','--');
xlabel(ax(end),'time $t/T$',ltx{:});
ax(end).LabelFontSizeMultiplier=1.2;
%%
exportgraphics(figure(2),'../figures/brus_sympo.pdf')
