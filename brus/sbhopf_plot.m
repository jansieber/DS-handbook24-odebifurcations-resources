%% plot SB continuation
clear
ty=@(r,s)strcmp(r{:,'TYPE'},s);
alw={'linewidth',2};
lw={'linewidth',3};
txt={'FontSize',20,'FontName','courier','FontWeight','bold'};
ltx={'Interpreter','latex'};
runSB=coco_bd_table('runSB');
runHB=coco_bd_table('runHB');
EPstart=runSB(ty(runSB,'EP')&runSB.PT==0,:);
startlab=EPstart.LAB{1};
prob=ode_SN2SN(coco_prob(),'','runSB',startlab);
[data,uidx,u]=coco_get_func_data(prob,'po.orb.coll.var','data','uidx','u0');
[cs,cv]=deal(data.pr.coll_seg,data.pr.coll_var);
t=cs.mesh.tbp;
x=reshape(u(cv.xbp_idx),cv.xbp_shp);
v=reshape(u(cv.v_idx),cv.xbp_shp);
[x,v]=deal(x(1:2:end,:),v(1:2:end,:));
figure(3);tl=tiledlayout(1,2);
tl.TileSpacing='tight';
clr=lines();
comp_order={'ColorOrder',clr};
tl.TileSpacing='tight';
nexttile();ax1=gca;hold(ax1,'on');
plsb=plot(ax1,runSB.lambda,runSB.B,'Color',clr(1,:),'DisplayName',...
    sprintf('Symmetry breaking:\n$(1\\,2\\,4\\,3)_4\\mapsto(1\\,4)_2(2\\,3)_2$'),lw{:});
plhb=plot(ax1,runHB.lambda,runHB.B,'Color',clr(2,:),'DisplayName',...
    sprintf('Equivariant Hopf\nbifurcation'),lw{:});
GHB=runSB(abs(runSB.amplitude)<1e-4,:);
tmarker={'kd','MarkerFaceColor',clr(2,:),alw{:},'MarkerSize',10};
plep=plot(ax1,EPstart.lambda(1),EPstart.B(1),tmarker{:},'DisplayName',...
    sprintf('Starting point'),alw{:});
plghb=plot(ax1,GHB.lambda(1),GHB.B(1),'ks','MarkerFaceColor',clr(3,:),'DisplayName',...
    sprintf('Degenerate equivariant\nHopf bifurcation'),alw{:},'MarkerSize',10);
legend(ax1,[plsb,plhb,plep,plghb],txt{:},ltx{:},'Location','NorthEast');
ax1.XLim=[min(runSB.lambda),max(runSB.lambda)*1.05];
ylim(ax1,[5.5,7.5]);
set(ax1,txt{:},'box','on',alw{:});
xlabel(ax1,'$\epsilon\times10^{3}$',ltx{:})
ylabel(ax1,'$B$',ltx{:})
ax1.XTick=15:17;
ax1.YTick=5:7;
ax1.LabelFontSizeMultiplier=1.2;
nexttile();ax2=gca;hold(ax2,'on');
tpl=plot(ax2,t,x,lw{:});
for k=1:4
    set(tpl(k),'DisplayName',sprintf('$x_%d$',k));
end
pmark=plot(ax2,NaN,NaN,tmarker{:},'DisplayName','');
dummy=plot(ax2,NaN,NaN,'wo','DisplayName','');
ax2.ColorOrderIndex=1;
vpl=plot(ax2,t,v*0.95,'-.',lw{:});
for k=1:4
    set(vpl(k),'DisplayName',sprintf('$\\delta^x_%d$',k));
end
legend(ax2,[pmark;tpl;dummy;vpl],'Location','eastoutside',txt{:},'FontSize',24,ltx{:});
yline(ax2,1,'k-','LineWidth',1,'HandleVisibility','off');
set(ax2,txt{:},alw{:},'box','on');
xlabel(ax2,'time $t/T$',ltx{:});
ax2.LabelFontSizeMultiplier=1.2;
ax2.YTick=-1:4;
ax2.XTick=linspace(0,1,5);
ax2.YLim=[-1,2.5];
%%
exportgraphics(figure(3),'../figures/brus_2dbif.pdf')

