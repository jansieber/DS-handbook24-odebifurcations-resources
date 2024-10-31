%% test branching off from Hopf bifurcation
clear
%% define r.h.s.
if exist('sys_brus', 'file')~=2
  syms y1 y2 y3 y4 y5 y6 y7 y8 A B lambda
  F = sco_sym2funcs([
    A-(B+1)*y1+y1^2*y2+lambda/1000*(-4*y1+y1+y3+y5+y7); ...
    B*y1-y1^2*y2+lambda/100*(-4*y2+y2+y4+y6+y8); ...
    A-(B+1)*y3+y3^2*y4+lambda/1000*(-4*y3+y1+y3+y5+y7); ...
    B*y3-y3^2*y4+lambda/100*(-4*y4+y2+y4+y6+y8); ...
    A-(B+1)*y5+y5^2*y6+lambda/1000*(-4*y5+y1+y3+y5+y7); ...
    B*y5-y5^2*y6+lambda/100*(-4*y6+y2+y4+y6+y8); ...
    A-(B+1)*y7+y7^2*y8+lambda/1000*(-4*y7+y1+y3+y5+y7); ...
    B*y7-y7^2*y8+lambda/100*(-4*y8+y2+y4+y6+y8)], ...
    {[y1; y2; y3; y4; y5; y6; y7; y8], [A; B; lambda]}, ...
    {'x', 'p'}, 'filename', 'sys_brus');
  funcs = {F(''), F('x'), F('p')};
else
  F = sco_gen(@sys_brus);
  funcs = {F(''), F('x'), F('p')};
end
%% Set common parameters
basedim=2;
n_osc=4;
alldim=n_osc*basedim;
A = 2;
B = 5.9;
%% create initial periodic solution guess
delta = B-1-A^2;
omega = sqrt(11*A^2*(11-9*delta)-100*delta^2)/11;

period = 2*pi/omega;
radius = 1e-4;

mu = A^2;
nu = 1i*omega-A^2-10*delta/11;

vec = [mu; nu; (1i)^3*mu; (1i)^3*nu; 1i*mu; 1i*nu; (1i)^2*mu; (1i)^2*nu];
re = real(vec);
im = imag(vec);

t0 = 0:period/400:period;
xhopf=[A; B/A; A; B/A; A; B/A; A; B/A];
[pnames,phopf]=deal(...
    {'A', 'B',   'lambda'},...
     [A;   B; 1000*delta/44]);
parbd={{'lambda', 'po.period' 'amplitude'},...
        {[1 30],   [0,100],    [0,1]}};
%%
cycles1={[1,2,4,3],4};
prob=ode_isol2ep(coco_prob(),'',funcs{:},xhopf,pnames,phopf);
prob=hopf_equiv_po(prob,'',omega,n_osc,cycles1,radius);
prob = coco_set(prob, 'cont', 'PtMX', [100 0], 'NPR', 100,...
    'NAdapt', 1, 'norm', inf,'al_max',90,'h0',1e-3);
prob = coco_set(prob, 'corr','TOL',1e-5); % JS! initial correction failed
coco(prob, 'hopf_eq1', [], 1, parbd{:});
%% plot of hopf_eq1
lw={'linewidth',2};
txt={'FontSize',20,'FontName','courier','FontWeight','bold'};
ltx={'Interpreter','latex'};
ty=@(r,s)strcmp(r{:,'TYPE'},s);
hopf_eq1=coco_bd_table('hopf_eq1');
figure(1);clf;ax1=gca;hold(ax1,'on');
plot(ax1,hopf_eq1.lambda,hopf_eq1.amplitude,'o-',lw{:});
ust1labs=coco_bd_labs('hopf_eq1','UST');
ust1sel=ty(hopf_eq1,'UST');
plot(ax1,hopf_eq1.lambda(ust1sel),hopf_eq1.amplitude(ust1sel),'ks','DisplayName','BP/UST',lw{:});
hopf_eq1prof=po_read_solution('hopf_eq1',ust1labs(1));
figure(2);clf;ax2=gca;
t1profs=plot(ax2,hopf_eq1prof.tbp,hopf_eq1prof.xbp(:,1:2:end),lw{:});
for i=1:4
 set(t1profs(i),'DisplayName',sprintf('$x_%d$',i));
end
legend(ax2,t1profs,txt{:},ltx{:});
%%
cycles2={[2,1],2; [3,4],1};
prob=ode_isol2ep(coco_prob(),'',funcs{:},xhopf,pnames,phopf);
prob=hopf_equiv_po(prob,'',omega,n_osc,cycles2,radius);
prob = coco_set(prob, 'cont', 'PtMX', [100 0], 'NPR', 200, 'NAdapt', 1, 'norm', inf);
prob = coco_set(prob, 'corr','TOL',1e-5); % JS! initial correction failed
coco(prob, 'hopf_eq2', [], 1, parbd{:});
%%
hopf_eq2=coco_bd_table('hopf_eq2');
plot(ax1,hopf_eq2.lambda,hopf_eq2.amplitude,lw{:});
pl2labs=coco_bd_labs('hopf_eq2','FP');
ust2sel=ty(hopf_eq2,'UST')|ty(hopf_eq2,'BP');
plot(ax1,hopf_eq2.lambda(ust2sel),hopf_eq2.amplitude(ust2sel),'ks','DisplayName','BP/UST',lw{:});
pl2labs=coco_bd_labs('hopf_eq2','UST');
hopf_eq2prof=po_read_solution('hopf_eq2',pl2labs(1));
figure(3);clf;ax3=gca;
t2profs=plot(ax3,hopf_eq2prof.tbp,hopf_eq2prof.xbp(:,1:2:end),lw{:});
set(t2profs(4),'LineStyle','--');
for i=1:4
 set(t2profs(i),'DisplayName',sprintf('$x_%d$',i));
end
legend(ax3,t2profs,txt{:},ltx{:});
%%
cycles3={[2,3,4],1};
prob=ode_isol2ep(coco_prob(),'',funcs{:},xhopf,pnames,phopf);
prob=hopf_equiv_po(prob,'',omega,n_osc,cycles3,radius);
prob = coco_set(prob, 'cont', 'PtMX', [0,90], 'NPR', 200, 'NAdapt', 1, 'norm', inf);
prob = coco_set(prob, 'corr','TOL',1e-5); % JS! initial correction failed
coco(prob, 'hopf_eq3', [], 1, parbd{:});
%%
hopf_eq3=coco_bd_table('hopf_eq3');
plot(ax1,hopf_eq3.lambda,hopf_eq3.amplitude,lw{:});
ust3labs=coco_bd_labs('hopf_eq3','UST');
ust3sel=ty(hopf_eq3,'UST')|ty(hopf_eq3,'BP');
plot(ax1,hopf_eq3.lambda(ust3sel),hopf_eq3.amplitude(ust3sel),'ks','DisplayName','BP/UST',lw{:});
ust3labs=coco_bd_labs('hopf_eq3','UST');
hopf_eq3prof=po_read_solution('hopf_eq3',ust3labs(2));
figure(4);clf;ax4=gca;
t3profs=plot(ax4,hopf_eq3prof.tbp,hopf_eq3prof.xbp(:,1:2:end),lw{:});
set(t3profs(4),'LineStyle','--');
set(t3profs(2),'LineStyle',':','LineWidth',4);
for i=1:4
 set(t3profs(i),'DisplayName',sprintf('$x_%d$',i));
end
legend(ax4,t3profs,txt{:},ltx{:});