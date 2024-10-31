clear
lw={'linewidth',2};
ltx={'interpreter','latex','FontSize',22};
%% create data and plots from cstr
bdpc=coco_bd_read('po_run');
bdp=cell2table(bdpc(2:end,:),'VariableNames',bdpc(1,:));
po_ep=coco_bd_labs('po_run','EP');
clr=lines();
figure(2);clf;ax2=gca;
%plot(bdp.('sigma'),bdp.('po.period'),'o-')
yyaxis(ax2,'left');
plot(ax2,bdp{:,'PT'}+1,bdp.('po.orb.NTST'),'o',lw{:},'DisplayName','mesh size $L$');
hold(ax2,'on');
plot(ax2,bdp{:,'PT'}+1,-10*log10(bdp.('po.orb.coll.err')),'+','color',clr(1,:),lw{:},...
    'DisplayName','-10\,log10(error)');
ax2.YLim(1)=0;
yyaxis(ax2,'right');
plot(ax2,bdp{:,'PT'}+1,bdp.('po.period'),'o',lw{:},'DisplayName','period $T$');
legend(ax2,ltx{:},'Location','southeast');
grid(ax2,'on');
set(ax2,'FontSize',18,'FontName','courier','FontWeight','bold','LineWidth',1);
xlabel(ax2,'point number along branch');
%%
exportgraphics(figure(2),'../figures/po_discretization.pdf',...
    'BackgroundColor','none','ContentType','vector')

