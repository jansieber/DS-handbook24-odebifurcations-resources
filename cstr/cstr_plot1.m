clear
doplot=false;
%% collect data from runs, define functions, labels
clr=lines();
ty=@(r,s)strcmp(r{:,'TYPE'},s);
ep=coco_bd_table('ep_run');
epeqfun=@(s,name)ep{ep.('ep.test.USTAB')==s,name};
epbiffun=@(s,name)ep{strcmp(ep.TYPE,s),name};
% bifurcations: saddle nodes, hopf, codim 2
sn1=coco_bd_table('sn_run');
hb1=coco_bd_table('hb_run1');
% sn formula (for comparison, not plotted)
xrg=linspace(0.01,0.99,1000);
sne=exp(1./(xrg-1));
[snd,sns]=deal(xrg.*sne./(1-xrg), sne./(1-xrg).^2);
hbe=exp(-xrg.*(1-xrg)./(xrg.*(1-xrg)-sn1.gamma(1)));
[hbd,hbs]=deal(xrg.*hbe./(1-xrg), xrg.*hbe./(xrg.*(1-xrg)-sn1.gamma(1)));
% find BTP & determine what is Hopf bifurcation
ibtp=find(strcmp(hb1.TYPE,'BTP'));
nhb=length(hb1.sigma);
nsa1sel=hb1.("ep.test.BTP")<1e-3&(1:nhb)'<=ibtp(1);
nsa2sel=hb1.("ep.test.BTP")<0&(1:nhb)'>ibtp(2);
nsa1=hb1(nsa1sel,:);
nsa2=hb1(nsa2sel,:);
hbsel=hb1.("ep.test.BTP")>-1e-3;
% Hopf surface - extract data
hb2=coco_bd_table('hb_run2');
dh1=coco_bd_table('dh_run1');
dh2=coco_bd_table('dh_run2');
Fsigma=scatteredInterpolant(hb2.x,hb2.beta,hb2.sigma,'linear','none');
Fdelta=scatteredInterpolant(hb2.x,hb2.beta,hb2.delta,'linear','none');
[xm,bm]=meshgrid(0.2:0.001:0.9,0:0.0001:0.06);
sd=cat(3,Fsigma(xm,bm),Fdelta(xm,bm));
% periodic orbits
pf=coco_bd_table('po_sn_run'); %pofold
pf=pf(find(strcmp(pf.TYPE,'FP')):-1:1,:); % po DH
hom=coco_bd_table('po_long_run'); % homoclinics
hom=hom(hom.sigma>0.4,:);
homvnames={'sigma','delta','xeq','yeq','det','tr'};
homfine=hom(:,homvnames);
hom=hom(hom.sigma>0.4,:);
% interpolate homoclinic line
thom=table2array(homfine);
thint=num2cell(interp1((1:size(thom,1))',thom,1:0.1:size(thom,1),'linear'),1);
homfine=table(thint{:},'VariableNames',homvnames);
is_snic=abs(homfine.det)<1e-1;
orig=[homfine.sigma(end);homfine.delta(end)];
bd=find(diff(is_snic));
% periodic orbit families (single parameter)
h2snic=coco_bd_table('po_run');
h2hom=coco_bd_table('po_run2');
%% line specs, marker specs
ltx={'interpreter','latex'};
txt={'Fontsize',18,'FontName','Courier','FontWeight','bold'};
lw={'linewidth',3};
lspec.sn={'-','LineWidth',4,'DisplayName','Saddle-Node Bifurcations', ...
    'color',clr(3,:)};
lspec.hb={'-','LineWidth',3,'DisplayName','Hopf bifurcation', ...
    'color',0.9*[1,0,0]};
lspec.nsa={':',  'LineWidth',2,'DisplayName','neutral saddle',...
    'color','r'};
lspec.snic={'-','color',[0,0,0.5],lw{:},'DisplayName','SNIC'};
lspec.hom={'-','color',clr(1,:),lw{:},'DisplayName','Homoclinic to saddle'};
lspec.snpo={':',lw{:},'color',clr(7,:),'LineWidth',6,'DisplayName','Saddle-node of periodic orbit'};
lspec.snpo_con={':',lw{:},'color',(clr(7,:)+1)/2,'LineWidth',6,'DisplayName','Saddle-node of p. o. (conjectured)'};
mspec=ep_plot_theme('ep');
[mspec.HB{5},mspec.SN{5}]=deal(12); % reset markersize
mspec.saddle={'+',lw{:},'MarkerSize',12,...
    'color',clr(2,:),'MarkerEdgeColor',clr(2,:),'MarkerFaceColor',clr(2,:),...
    'DisplayName','Saddle'};
mspec.ep={'o','MarkerSize',10,'MarkerEdgeColor',clr(1,:),'MarkerFaceColor',clr(1,:),'DisplayName','End Point'};
mspec.cusp={'o'  'MarkerFaceColor' clr(3,:)   'MarkerSize' 13  'DisplayName'  'Cusp' ...
    'Linewidth',0.1,'MarkerEdgeColor',clr(3,:),'DisplayName','Cusp'};
mspec.dh={'ko'  'MarkerFaceColor'  'b'  'MarkerSize'  13  'DisplayName'  'Degenerate Hopf'};
mspec.bp={'ks'  'MarkerFaceColor' clr(4,:)   'MarkerSize'  13  'DisplayName'  'Branch Point'};
mspec.btp={'kp','MarkerFaceColor','g','MarkerSize',15,'DisplayName',...
    'Bogdanov-Takens bifurcation','LineWidth',1.5};
mspec.hnsa={'o','MarkerEdgeColor',clr(4,:),'MarkerFaceColor',clr(4,:),'MarkerSize',12,'DisplayName','Homoclinic to neutral saddle'};
mspec.ncs={'ks','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',12,'DisplayName',sprintf('Noncentral connection to s.\\,n.')};
%%
figure(1);clf;tiledlayout(1,2,'TileSpacing','tight');
nexttile;ax1=gca;
hold(ax1,'on');
plot3(ax1,epeqfun(0,'x'),epeqfun(0,'sigma'),epeqfun(0,'delta'),'-','color',clr(1,:),...
    lw{:},'DisplayName',sprintf('equilibria, USTAB=%d',0));
plot3(ax1,epeqfun(1,'x'),epeqfun(1,'sigma'),epeqfun(1,'delta'),'--','color',clr(1,:),...
    lw{:},'DisplayName',sprintf('equilibria, USTAB=%d',1));
plot3(ax1,epeqfun(2,'x'),epeqfun(2,'sigma'),epeqfun(2,'delta'),':','color',clr(1,:),...
    lw{:},'DisplayName',sprintf('equilibria, USTAB=%d',2));
plot3(ax1,epbiffun('HB','x'),epbiffun('HB','sigma'),epbiffun('HB','delta'),mspec.HB{:})
plot3(ax1,epbiffun('SN','x'),epbiffun('SN','sigma'),epbiffun('SN','delta'),mspec.SN{:})
plot3(ax1,epbiffun('EP','x'),epbiffun('EP','sigma'),epbiffun('EP','delta'),'o','MarkerSize',10,...
    'MarkerEdgeColor',clr(1,:),'MarkerFaceColor',clr(1,:),'DisplayName','End Point')
grid(ax1,'on')
legend(ax1,'location','best',ltx{:},'FontSize',20);
set(ax1,txt{:},'LineWidth',2);
xlabel(ax1,'$x$',ltx{:},'Fontsize',30);
ylabel(ax1,'$\sigma$',ltx{:},'Fontsize',30);
zl=zlabel(ax1,'$\delta$',ltx{:},'Fontsize',30,'Rotation',0);
view(ax1,[-38,25]);
% new figure
nexttile(2);ax2=gca;
hold(ax2,'on');
% plot SN bifurcation
pl.sn=plot(ax2,sn1.sigma,sn1.delta,lspec.sn{:});
% plot Hopf bifurcation
pl.hb=plot(ax2,hb1.sigma(hbsel),hb1.delta(hbsel),lspec.hb{:});
pl.nsa=plot(ax2,hb1.sigma(nsa1sel),hb1.delta(nsa1sel),lspec.nsa{:});
plot(ax2,hb1.sigma(nsa2sel),hb1.delta(nsa2sel),lspec.nsa{:});
% plot special points
[icusp, idh, ibp]=deal(find(ty(sn1,'FP')),find(ty(hb1,'DH')),find(ty(hb1,'BP')));
pl.cusp=plot(ax2,sn1.sigma(icusp),sn1.delta(icusp),mspec.cusp{:});
pl.btp=plot(ax2,hb1.sigma(ibtp),hb1.delta(ibtp),mspec.btp{:});
pl.dh=plot(ax2,hb1.sigma(idh),hb1.delta(idh),mspec.dh{:});
pl.bp=plot(ax2,hb1.sigma(ibp),hb1.delta(ibp),mspec.bp{:});
legend(ax2,[pl.sn,pl.hb,pl.nsa,pl.cusp,pl.btp,pl.dh,pl.bp],...
    'location','best',ltx{:},'Fontsize',20)
set(ax2,txt{:},'LineWidth',2,'Box','on');
xlim(ax2,[0,1.05])
yl=ylabel(ax2,'$\delta$',ltx{:},'Fontsize',30,'Rotation',0);
xlabel(ax2,'$\sigma$',ltx{:},'Fontsize',30);
yl.Position(1:2)=[-0.08,0.23];
an1=annotation('textbox','String','(a)','Interpreter','latex',...
    'LineStyle','none','FontSize',30);
an1.Position(1:2)=[0.02,0.85];
ax2.YTick=0:0.1:0.3;
an2=annotation('textbox','String','(b)','Interpreter','latex',...
    'LineStyle','none','FontSize',30);
an2.Position(1:2)=[0.48,0.85];
%%
export(doplot,@svg2pdf,1,'cstr_ep');
%% plot surfaces
figure(3);clf;tiledlayout(1,2,'TileSpacing','compact');
atlas = coco_bd_read('hb_run2', 'atlas');
nexttile(1);ax8=gca;
plot_atlas_kd(atlas.charts, 1,3,5)
hold(ax8,'on');
thm =  struct('lspec', {{{'b', 'LineWidth', 3},{'b-.', 'LineWidth', 3}}},'ustab','');
h0pl=plot3(ax8,hb1.sigma,hb1.x,hb1.beta,'color','r',lw{:},'DisplayName','Hopf at $\beta=0$');
legend(ax8,h0pl,'Location','best','Interpreter','latex',txt{:})
grid(ax8,'on');
axis(ax8,[0 inf 0.2 0.9 0 inf]);
view(ax8,[-120,40])
hold(ax8,'off');
set(ax8,txt{:},'linewidth',2);
xlabel(ax8,'$\sigma$',ltx{:},'FontSize',24);
ylabel(ax8,'$x$',ltx{:},'FontSize',24);
zlabel(ax8,'$\beta$','Rotation',0,ltx{:},'FontSize',24);
%
%%
is_old=loc_isold()
%%
nexttile(2);ax9=gca;
colormap('default');
hold(ax9,'off');
dhpl=plot(ax9,dh1.sigma,dh1.beta,'color','r',lw{:},'DisplayName','degen. Hopf');
hold(ax9,'on');
[cs,hs]=contourf(ax9,sd(:,:,1),bm,sd(:,:,2),[0.08,0.12,0:0.05:0.3],'ShowText','on',...
    'EdgeColor','k','FaceAlpha',0.3,'EdgeAlpha',0.8,'LineWidth',2);
clabel(cs,hs, 'FontSize', 15, 'FontWeight','bold','FontName','Courier');
plot(ax9,dh2.sigma,dh2.beta,'color','r',lw{:});
cbs=colorbar(ax9);
cbs.Title.String='$\delta$';
cbs.Title.Interpreter='latex';
cbs.Title.FontSize=24;
set(ax9,txt{:},'LineWidth',2)
ax9.YTick=0:0.01:0.05;
xlim(ax9,[0.35,1]);
xlabel(ax9,'$\sigma$',ltx{:},'FontSize',24);
ylabel(ax9,'$\beta$',ltx{:},'Rotation',0,'FontSize',24);
legend(ax9,dhpl,'Location','northeast',ltx{:},txt{:})
an3=annotation('textbox','String','(a)','Interpreter','latex',...
    'LineStyle','none','FontSize',25);
an3.Position(1:2)=[0.048,0.88];
an4=annotation('textbox','String','(b)','Interpreter','latex',...
    'LineStyle','none','FontSize',25);
an4.Position(1:2)=[0.5,0.88];
%%
export(doplot,@exportgraphics,figure(3),'../Figures/cstr_hopf2.pdf','Resolution',600);
%% periodic orbits, connecting to SNIC and homoclinic
[sigbd_snic,tick_snic]=sigtrafo(h2snic,[2,3],0.01);
[sigbd_hom,tick_hom]=sigtrafo(h2hom,[0.1,1],0.002);
[xtickpos,xtickperm]=sort([tick_hom.pos,tick_snic.pos]);
xticklab=[tick_hom.label,tick_snic.label];
xticklab=arrayfun(@(l){num2str(l)},xticklab(xtickperm));
%
poend=h2snic(strcmp(h2snic.TYPE,'EP'),:);
figure(2);clf;tiledlayout(1,4,'TileSpacing','tight');
% plot amplitude
nexttile;ax3=gca;
hold(ax3,'on');
pothm=po_plot_theme('po');
plot(ax3,sigbd_snic,h2snic{:,'amplitude'},'-','Color',clr(1,:),...
    lw{:},'DisplayName','per. orb. to SNIC');
hold(ax3,'on');
plot(ax3,sigbd_hom,h2hom{:,'amplitude'},':','Color',clr(1,:),...
    lw{:},'DisplayName','per. orb. to homoclinic');
plot(ax3,sigbd_snic(1),h2snic{1,'amplitude'},mspec.ep{:})
plot(ax3,sigbd_snic(end),h2snic{end,'amplitude'},mspec.HB{:})
plot(ax3,sigbd_hom(1),h2hom{1,'amplitude'},mspec.ep{:},'HandleVisibility','off');
plot(ax3,sigbd_hom(end),h2hom{end,'amplitude'},mspec.HB{:},'HandleVisibility','off');
text(ax3,1.5,0,'//','VerticalAlignment','middle','HorizontalAlignment','center',txt{:});
set(ax3,txt{:},lw{:},'box','on');
xlabel(ax3,'$\sigma$','Fontsize',30,'Interpreter','latex');
ylabel(ax3,'amplitude ($x$)','Fontsize',24,'Interpreter','latex');
legend(ax3,'location','northeast','Interpreter','latex',txt{:})
ylim(ax3,[0,1.2])
ax3.YTick=0:0.2:0.8;
ax3.XTick=xtickpos;
ax3.XTickLabel=xticklab;
% plot sigma vs period
nexttile;ax4=gca;
plot(ax4,sigbd_snic,h2snic{:,'po.period'},'-','Color',clr(1,:),...
    lw{:},'DisplayName','per. orb. to SNIC');
hold(ax4,'on');
plot(ax4,sigbd_hom,h2hom{:,'po.period'},':','Color',clr(1,:),...
    lw{:},'DisplayName','per. orb. to homoclinic');
plot(ax4,sigbd_snic(1),h2snic{1,'po.period'},mspec.ep{:})
plot(ax4,sigbd_hom(1),h2hom{1,'po.period'},mspec.ep{:},'HandleVisibility','off')
plot(ax4,sigbd_snic(end),h2snic{end,'po.period'},mspec.HB{:})
plot(ax4,sigbd_hom(end),h2hom{end,'po.period'},mspec.HB{:},'HandleVisibility','off');
text(ax4,1.5,0,'//','VerticalAlignment','middle','HorizontalAlignment','center',txt{:});
set(ax4,txt{:},lw{:})
legend(ax4,'location','northeast','Interpreter','latex',txt{:})
ylim(ax4,[0,35])
xlabel(ax4,'$\sigma$',ltx{:},'Fontsize',30);
ylabel(ax4,'period ($T$)','Fontsize',24,'Interpreter','latex');
ax4.YTick=0:5:20;
ax4.XTick=xtickpos;
ax4.XTickLabel=xticklab;
% plot phase portrait of SNIC
nexttile;ax5=gca;
polab=coco_bd_labs('po_run','EP');
snic=po_read_solution('po_run',polab(1));
po2lab=coco_bd_labs('po_run2','EP');
homorb=po_read_solution('po_run2',po2lab(1));
dxdtnormsnic=diff(snic.xbp,[],1)./diff(snic.tbp);
[~,iminsnic]=min(max(abs(dxdtnormsnic)));
dxdtnormhom=diff(homorb.xbp,[],1)./diff(homorb.tbp);
[~,iminhom]=min(max(abs(dxdtnormhom)));
plot(ax5,snic.xbp(:,1),snic.xbp(:,2),'-','Linewidth',4,'color',clr(1,:),...
    'DisplayName','SNIC');
hold(ax5,'on')
plot(ax5,homorb.xbp(:,1),homorb.xbp(:,2),':',lw{:},'color',clr(6,:),...
    'DisplayName',sprintf('Homoclinic\nto saddle'));
plot(ax5,snic.xbp(iminsnic,1),snic.xbp(iminsnic,2),mspec.SN{:},'DisplayName','Saddle-node');
plot(ax5,homorb.xbp(iminhom,1),homorb.xbp(iminhom,2),mspec.saddle{:});
set(ax5,txt{:},'box','on',lw{:});
xlim(ax5,[0,1])
ylim(ax5,[0,7])
xlabel(ax5,'$x$',ltx{:},'Fontsize',30);
ylabel(ax5,'$y$',ltx{:},'Fontsize',30,'Rotation',0);
legend(ax5,'location','NorthWest','interpreter','latex',txt{:})
ax5.YTick(end)=[];
% plot time profile of SNIC
nexttile;
yyaxis left;ax6=gca;
plot(ax6,snic.tbp(:,1),snic.xbp(:,1),'-',lw{:},'color',clr(1,:),...
    'DisplayName','SNIC $x(t)$');
hold(ax6,'on');
[thom,ihom]=sort(mod(homorb.tbp(:,1)-1,homorb.T));
plot(ax6,thom,homorb.xbp(ihom,1),':',lw{:},'color',clr(1,:),...
    'DisplayName','Homoclinic $x(t)$');
ax6.YLim=[0,1];
ax6.YTick=0:0.2:0.8;
set(ax6,txt{:},'box','on',lw{:});
yyaxis right;ax7=gca;
plot(ax7,snic.tbp(:,1),snic.xbp(:,2),'-',lw{:},'color',clr(2,:),...
    'DisplayName','SNIC $y(t)$');
hold(ax7,'on');
plot(ax7,thom,homorb.xbp(ihom,2),':',lw{:},'color',clr(2,:),...
    'DisplayName','homoclinic $y(t)$');
set(ax7,txt{:},'box','on',lw{:});
xlabel(ax7,'$t$',ltx{:},'Fontsize',30);
%ylabel(ax6,'$x$',ltx{:},'Fontsize',30,'Rotation',0);
legend(ax6,'location','NorthWest','interpreter','latex',txt{:})
ax7.YLim=[0,7];
ax7.YTick=0:2:6;
%
an3=annotation('textbox','String','(a)','Interpreter','latex',...
    'LineStyle','none','FontSize',25);
an3.Position(1:2)=[0.04,0.88];
an4=annotation('textbox','String','(b)','Interpreter','latex',...
    'LineStyle','none','FontSize',25);
an4.Position(1:2)=[0.27,0.88];
an5=annotation('textbox','String','(c)','Interpreter','latex',...
    'LineStyle','none','FontSize',25);
an5.Position(1:2)=[0.5,0.88];
an6=annotation('textbox','String','(d)','Interpreter','latex',...
    'LineStyle','none','FontSize',25);
an6.Position(1:2)=[0.73,0.88];
%%
export(doplot,@exportgraphics,figure(2),'../Figures/cstr_po1.pdf','BackgroundColor','none','ContentType','vector')
%
%% zoomed-in 2D bifurcation diagram
fg=figure(4);clf;axpb=gca;
hold(axpb,'on');
rotang=150;
smat=[cosd(rotang),sind(rotang); -sind(rotang),cosd(rotang)];
trafo=@(sd)[sd(:,1)-orig(1),sd(:,2)-orig(2)]*smat;%/rotmat;
trhb=trafo([hb1.sigma,hb1.delta]);
trnsa=trhb(nsa1sel,:);
trhb=trhb(hbsel,:);
trsn=trafo([sn1.sigma,sn1.delta]);
trpf=trafo([pf.sigma,pf.delta]);
trhomfine=trafo([homfine.sigma,homfine.delta]);
nsahom=trafo(hom{ty(hom,'NSA')&hom.amplitude>1e-4,{'sigma','delta'}});
ncshom=trafo(hom{ty(hom,'NCS')&hom.amplitude>1e-4,{'sigma','delta'}});
[~,i_nsahom_close]=min(sum(nsahom.^2,2),[],1);
icusp=find(strcmp(sn1.TYPE,'FP'));
trh2snic=[trafo(h2snic{:,{'sigma','delta'}}),h2snic{:,'po.period'}];
% plot coordinate system
frc=0.001;
sigma_unit=trafo(orig'+repmat([0,frc,1],2,1)'*diag([1,0]));
plot(axpb,sigma_unit(:,1),sigma_unit(:,2),'+:','LineWidth',2,'color',0.5*[1,1,1],'DisplayName','$\sigma$ axis')
delta_unit=trafo(orig'+repmat([0,frc,1],2,1)'*diag([0,1]));
plot(axpb,delta_unit(:,1),delta_unit(:,2),'+--','LineWidth',2,'color',0.5*[1,1,1],'DisplayName','$\delta$ axis')
text(sigma_unit(2,1),sigma_unit(2,2),'\ $(\sigma_\mathrm{BT}+10^{-3},\delta_\mathrm{BT})$',...
    'VerticalAlignment','middle','HorizontalAlignment','left',txt{:},ltx{:});
text(sigma_unit(1,1),sigma_unit(1,2),'$(\sigma_\mathrm{BT},\delta_\mathrm{BT})$',...
    'VerticalAlignment','top','HorizontalAlignment','center',txt{:},ltx{:});
text(delta_unit(2,1),delta_unit(2,2),'$\ (\sigma_\mathrm{BT},\delta_\mathrm{BT}+10^{-3})$',...
    'VerticalAlignment','middle','HorizontalAlignment','left',txt{:},ltx{:});
% plot 1-par family HOpf 2 SNIC
plot(axpb,trh2snic(:,1),trh2snic(:,2),'k:','color',[0,0,0],lw{:},'DisplayName','p.o. family Hopf to SNIC');
% plot codimension-1 bifurcations of equilibria
plot(axpb,trhb(:,1),trhb(:,2),lspec.hb{:});
plot(axpb,trnsa(:,1),trnsa(:,2),lspec.nsa{:});
plot(axpb,trsn(:,1),trsn(:,2),lspec.sn{:});
% plot fold of periodic orbits
plot(axpb,trpf(:,1),trpf(:,2),lspec.snpo{:});
plot(axpb,[trpf(end,1);nsahom(i_nsahom_close,1)],[trpf(end,2);nsahom(i_nsahom_close,2)],lspec.snpo_con{:});
% plot homoclinics of both types
homclr=[[0,0,0.5];clr(1,:)];
hname={'SNIC','Homoclinic to saddle'};
[trsnic,trhsad]=deal(trhomfine);
trsnic(~is_snic,:)=NaN;
trhsad(is_snic,:)=NaN;
plot(axpb,trsnic(:,1),trsnic(:,2),lspec.snic{:});
plot(axpb,trhsad(:,1),trhsad(:,2),lspec.hom{:});
plot(axpb,nsahom(:,1),nsahom(:,2),mspec.hnsa{:});
plot(axpb,trhomfine(bd(1:end-1),1),trhomfine(bd(1:end-1),2),mspec.ncs{:});
plot(axpb,trsn(icusp,1),trsn(icusp,2),mspec.cusp{:})
plot(axpb,trnsa(end,1),trnsa(end,2),mspec.btp{:})
axpb.YLim=[min(trhomfine(bd(1:end-1),2))*1.45,max(trhomfine(bd(1:end-1),2))*1.1];
axpb.XLim=[1.1*trsn(icusp,1),0.05];%max(nsahom(:,1))*1.];
axpb.XTick=-0.01:0.01:0.05;
set(axpb,txt{:},'linewidth',2,'box','on')
lg=legend(axpb,'Location','eastoutside',ltx{:},txt{:});
xlabel(axpb,'$\tilde\sigma$',ltx{:},'FontSize',24)
ylabel(axpb,'$\tilde\delta$',ltx{:},'Fontsize',24,'Rotation',0)
%%
export(doplot,@exportgraphics,figure(4),'../Figures/cstr_pobif.pdf','ContentType','auto','BackgroundColor','none')
%% Plot Hopf bubbles

figure(5);clf;axhb=gca;hold(axhb,'on');
for i=1:4
    bubble{i}=coco_bd_table(sprintf('bubble%d',i));
end
hvtoggle={'off','on'};
for i=4:-1:1
    curve=table2array(bubble{i}(:,{'sigma','delta','amplitude'}));
    ustab=table2array(bubble{i}(:,'po.test.USTAB'));
    [cs,cu]=deal(curve);
    cs(ustab>0,:)=NaN;
    cu(ustab==0,:)=NaN;
    hv={'HandleVisibility',hvtoggle{double(i==1)+1}};
    cl=interp1([0,6]',[clr(1,:);[1,1,1]],i,'linear');
    plot3(axhb,cs(:,1),cs(:,2),cs(:,3),'-','color',cl,lw{:},'DisplayName',...
        sprintf('periodic orbit (stable), $\\delta=$%3.2f',bubble{i}.delta(1)));
    plot3(axhb,cu(:,1),cu(:,2),cu(:,3),':','color',cl,lw{:},'DisplayName','periodic orbit (unstable)',hv{:});
    ibhopf=find(abs(curve(:,3))<1e-2);
    plot3(axhb,curve(ibhopf,1),curve(ibhopf,2),curve(ibhopf,3),'o',...
        'markerfacecolor',cl,'MarkerSize',7,'HandleVisibility','off','MarkerEdgeColor','none');
    ibsn=find(strcmp(bubble{i}.TYPE,'SN'));
    plot3(axhb,curve(ibsn,1),curve(ibsn,2),curve(ibsn,3),'kd','DisplayName','Saddle-node of periodic orbit',...
        'markerfacecolor','g','MarkerSize',10,hv{:},'MarkerEdgeColor','k',lw{:});
end
plot3(axhb,hb1{:,'sigma'},hb1{:,'delta'},0*hb1{:,'delta'},'r',lw{:},'DisplayName','Hopf bifurcation');
idh=find(strcmp(hb1.TYPE,'DH'));
plot3(axhb,hb1{idh,'sigma'},hb1{idh,'delta'},0*hb1{idh,'delta'},mspec.dh{:});
ipfdh=find(strcmp(pf.TYPE,'FP'));
plot3(axhb,pf{1:ipfdh,'sigma'},pf{1:ipfdh,'delta'},pf{1:ipfdh,'amplitude'},':',...
    'color',clr(7,:),'LineWidth',6,'DisplayName','Saddle-node of periodic orbit');
grid(axhb,'on');
view(axhb,[-12,70]);
xlim(axhb,[0.5,1.03]);
ylim(axhb,[0.12,0.36]);
legend(axhb,'Location','EastOutside',ltx{:},'FontSize',20)
set(axhb,txt{:},'LineWidth',2);
xlabel(axhb,'$\sigma$',ltx{:},'FontSize',30)
ylabel(axhb,'$\delta$',ltx{:},'FontSize',30)
zlabel(axhb,'amplitude',ltx{:},'FontSize',24)
%%
export(doplot,@exportgraphics,figure(5),'../figures/cstr_hopfbubbles.pdf')
%% find crossing of neutral saddle with homoclinic
function nsahom=cross_nsa_hom(pfend,hom,nsa)
ppn=interp1(nsa(:,1),nsa(:,2),'linear','pp');
isgn=find(diff(sign(hom(:,2)-fnval(ppn,hom(:,1)))));
isgn=isgn(hom(isgn,1)>0);
[~,iminns]=min(sum((hom(isgn,:)-repmat(pfend,length(isgn),1)).^2,2));
nsahom=hom(isgn(iminns),:);
end
%%
function svg2pdf(fg,fname)
fullname=['../figures/',fname];
svg=[fullname,'.svg'];
pdf=[fullname,'.pdf'];
saveas(figure(fg),svg);
system(['inkscape --export-filename ',pdf,' ',svg])
end
%%
function [sig,ticks]=sigtrafo(bifd,base,scale)
bd=bifd{[1,end],'sigma'};
trafo=@(s)(s-bd(end))/(bd(end)-bd(1))*(base(end)-base(1))+base(end);
scale=-scale*sign(diff(bd));
ticks.label=round(bd(end),3):scale:round(bd(1),3);
ticks.pos=round(trafo(ticks.label),round(abs(log10(scale))));
sig=trafo(bifd{:,'sigma'});
end
%%
function export(doplot,exportfun,varargin)
if ~doplot
    return
end
feval(exportfun,varargin{:});
end
%%
function isold=loc_isold()
v=version;
Rloc=strfind(v,'R');
year=str2double(v(Rloc+(1:4)));
rel=v(Rloc+5);
isold=year<2022 || (year==2022 && rel=='b');
end