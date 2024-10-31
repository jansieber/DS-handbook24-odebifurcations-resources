%% runs for symmetric problem
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
A = 2;
B = 5.9;
[n_osc,basedim]=deal(4,2);
%% create initial periodic solution guess
delta = B-1-A^2;
omega = sqrt(11*A^2*(11-9*delta)-100*delta^2)/11;

period = 2*pi/omega;
epsilon = 1e-4;

mu = A^2;
nu = 1i*omega-A^2-10*delta/11;

vec = [mu; nu; (1i)^3*mu; (1i)^3*nu; 1i*mu; 1i*nu; (1i)^2*mu; (1i)^2*nu];
re = real(vec);
im = imag(vec);

t0 = 0:period/400:period;
x0 = [A; B/A; A; B/A; A; B/A; A; B/A] + epsilon*...
  (re.*cos(t0*2*pi/period)-im.*sin(t0*2*pi/period));
%% initial run (run1), branching off from equivariant Hopf bifurcation
prob = coco_prob;
prob = coco_set(prob, 'coll', 'MXCL', 'off');
prob = coco_set(prob, 'po', 'SN', 'off', 'PD', 'off', 'TR', 'off');
prob = ode_isol2po(prob, '', funcs{:}, t0', x0', {'A', 'B', 'lambda'}, [A; B; 1000*delta/44]); 
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
prob = coco_add_func(prob, 'amplitude', ...
  @(p,d,u) deal(d, max(u(1:2:end))-min(u(1:2:end))), ...
  struct('cid', 'po.orb.coll'), 'singular', ...
  'amplitude', 'uidx', uidx(data.coll_seg.maps.xbp_idx), ...
  'remesh', @ampremesh);
prob = po_construct_sbtest(prob, 'po', 5);
prob = coco_add_event(prob, 'UST', 'po.test.USTAB', (0:5-1)'+0.5);
prob = coco_set(prob, 'cont', 'PtMX', [100 0], 'NPR', 100, 'NAdapt', 1, 'norm', inf);
disp('Run1: branch off from Equiv. Hopf')
coco(prob, 'run1', [], 1, {'lambda' 'po.period' 'amplitude'}, [1 30]);
%% plot of run1
lw={'linewidth',2};
txt={'FontSize',20,'FontName','courier','FontWeight','bold'};
ltx={'Interpreter','latex'};
ty=@(r,s)strcmp(r{:,'TYPE'},s);
run1=coco_bd_table('run1');
figure(1);clf;ax1=gca;hold(ax1,'on');
plot(ax1,run1.lambda,run1.amplitude,lw{:});
bplabs=coco_bd_labs('run1','BP');
bpsel=ty(run1,'BP');
plot(ax1,run1.lambda(bpsel),run1.amplitude(bpsel),'ks','DisplayName','BP/UST',lw{:});
run1prof=po_read_solution('run1',bplabs(1));
figure(2);clf;ax2=gca;
tprofs=plot(ax2,run1prof.tbp,run1prof.xbp(:,1:2:end),lw{:});
for i=1:4
 set(tprofs(i),'DisplayName',sprintf('$x_%d$',i));
end
legend(ax2,tprofs,txt{:},ltx{:});
%% run2, branching off from (second) UST1 point near lambda=17 into 2x2 symmetry
UST1=coco_bd_labs('run1','UST');
prob = coco_prob;
prob = ode_BP2po(prob, '', 'run1', UST1(2));
prob = coco_set(prob, 'cont', 'branch', true);
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
prob = coco_add_func(prob, 'amplitude', ...
  @(p,d,u) deal(d, max(u(1:2:end))-min(u(1:2:end))), ...
  struct('cid', 'po.orb.coll'), 'singular', ...
  'amplitude', 'uidx', uidx(data.coll_seg.maps.xbp_idx), ...
  'remesh', @ampremesh);
prob = coco_add_event(prob, 'UST', 'po.test.USTAB', (0:5-1)'+0.5);
prob = coco_set(prob, 'cont', 'NPR', 0, 'PtMX', [0 201], 'NAdapt', 10, 'norm', inf, 'h_max', 10);
disp('Run2: branch off from secondary branch point, toward (14)_2 (23)_2 symmetry')
coco(prob, 'run2', [], 1, {'lambda', 'po.period' 'amplitude'}, [14 30]);
%% plot of run2
ty=@(r,s)strcmp(r{:,'TYPE'},s);
run2=coco_bd_table('run2');
figure(3);clf;ax3=gca;hold(ax3,'on');
plot(ax3,run2.lambda,run2.amplitude,lw{:});
specsel=ty(run2,'BP')|ty(run2,'UST');
plot(ax3,run2.lambda(specsel),run2.amplitude(specsel),'ks','DisplayName','BP/UST',lw{:});
fp2labs=coco_bd_labs('run2','FP');
run2prof=[po_read_solution('run2',fp2labs(1)),po_read_solution('run2',fp2labs(2))];
figure(4);clf;ax4=gca;hold(ax4,'on');
t2snprofs=plot(ax4,run2prof(1).tbp,run2prof(1).xbp(:,1:2:end),lw{:});
ax4.ColorOrderIndex=1;
t2saprofs1=plot(ax4,run2prof(2).tbp,run2prof(2).xbp(:,1),'o');
ax4.ColorOrderIndex=1;
t2saprofs=plot(ax4,run2prof(2).tbp,run2prof(2).xbp(:,1:2:end),':',lw{:});
for i=1:4
 set(t2snprofs(i),'DisplayName',sprintf('$x_%d$',i));
 set(t2saprofs(i),'DisplayName',sprintf('$x_%d$',i));
end
legend(ax4,[t2snprofs;t2saprofs],txt{:},ltx{:});
%% run 3 branch off from symmetry-adding bifurcation at lambda=14.5 into symmetry (1=2)+2
BP2=cell2mat(run2{ty(run2,'BP')&run2{:,'lambda'}<15,'LAB'});
prob = coco_prob;
prob = ode_BP2po(prob, '', 'run2', BP2(1));
prob = coco_set(prob, 'cont', 'branch', true);
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
prob = coco_add_func(prob, 'amplitude', ...
  @(p,d,u) deal(d, max(u(1:2:end))-min(u(1:2:end))), ...
  struct('cid', 'po.orb.coll'), 'singular', ...
  'amplitude', 'uidx', uidx(data.coll_seg.maps.xbp_idx), ...
  'remesh', @ampremesh);
prob = coco_add_event(prob, 'UST', 'po.test.USTAB', (0:5-1)'+0.5);
prob = coco_set(prob, 'cont', 'NPR', 0, 'PtMX', [-100 301], 'NAdapt', 10, 'norm', inf, 'h_max', 10);
disp('Run3: branch off from tertiary branch point adding symmetry, now (xy)_1 (uv)_2 symmetry')
coco(prob, 'run3', [], 1, {'lambda', 'po.period' 'amplitude'}, [1 30]);
%% plot run3
run3=coco_bd_table('run3');
figure(5);clf;ax5=gca;hold(ax5,'on');
plot(ax5,run3.lambda,run3.amplitude,lw{:});
specsel=ty(run3,'UST');
plot(ax5,run3.lambda(specsel),run3.amplitude(specsel),'ks','DisplayName','FP/BP/UST',lw{:});
%% run 4, branch off from symmetry-breaking bifurcation near lambda=19.6UST = coco_bd_labs('run3', 'UST');
UST3 = cell2mat(run3{ty(run3,'UST')&run3{:,'lambda'}>19&run3{:,'lambda'}<20,'LAB'});%coco_bd_labs('run3', 'UST');
prob = coco_prob;
prob = ode_BP2po(prob, '', 'run3', UST3);
prob = coco_set(prob, 'cont', 'branch', true);
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
prob = coco_add_func(prob, 'amplitude', ...
  @(p,d,u) deal(d, max(u(1:2:end))-min(u(1:2:end))), ...
  struct('cid', 'po.orb.coll'), 'singular', ...
  'amplitude', 'uidx', uidx(data.coll_seg.maps.xbp_idx), ...
  'remesh', @ampremesh);
prob = coco_add_event(prob, 'UST', 'po.test.USTAB', (0:5-1)'+0.5);
prob = coco_set(prob, 'cont', 'NPR', 0, 'PtMX', [0 51], 'NAdapt', 10, 'norm', inf, 'h_max', 10);
disp('Run4: branch off from quartiary branch point losing symmetry, now (xy)_1 symmetry')
coco(prob, 'run4', [], 1, {'lambda', 'po.period' 'amplitude'}, [1 30]);
%% Try and combine all periodic-orbit plots
ty=@(r,s)strcmp(r{:,'TYPE'},s);
run1=coco_bd_table('run1');
run2=coco_bd_table('run2');
run3=coco_bd_table('run3');
run4=coco_bd_table('run4');
figure(6);clf;axc=gca;hold(axc,'on');
clr=lines();
axc.ColorOrderIndex=1;
% run1 
plot(axc,run1.lambda,run1.amplitude,lw{:});
bp1sel=ty(run1,'BP');
plot(axc,run1.lambda(bp1sel),run1.amplitude(bp1sel),'ks','DisplayName','BP/UST',lw{:});
% run2
plot(axc,run2.lambda,run2.amplitude,lw{:});
spec2sel=ty(run2,'BP')|ty(run2,'UST');
plot(axc,run2.lambda(spec2sel),run2.amplitude(spec2sel),'ks','DisplayName','BP/UST',lw{:});
% run3
axc.ColorOrderIndex=4;
plot(axc,run3.lambda,run3.amplitude,lw{:});
spec3sel=ty(run3,'UST');
plot(axc,run3.lambda(spec3sel),run3.amplitude(spec3sel),'ks','DisplayName','FP/BP/UST',lw{:});
% run 4
plot(axc,run4.lambda,run4.amplitude,lw{:});
spec4sel=find(ty(run4,'EP'));
plot(axc,run4.lambda(spec4sel(2)),run4.amplitude(spec4sel(2)),'ro','DisplayName','EP',lw{:});
%% plot time profile of small-amplitude orbit of run3
r3sm=po_read_solution('run3',run3{find(run3.amplitude<1e-2&ty(run3,'FP')),'LAB'}{1});
figure(7);clf;axsm3=gca;hold(axsm3,'on');
tsm3profs=plot(axsm3,r3sm.tbp,r3sm.xbp(:,1:2:end),lw{:});

%% Track symmetry breaking
UST = coco_bd_labs('run1', 'UST1');
prob = coco_prob;
prob = coco_set(prob, 'coll', 'MXCL', 'off');
prob = coco_set(prob, 'po', 'SN', 'off', 'PD', 'off', 'TR', 'off');
prob = ode_SN2SN(prob, '', 'run1', UST(2));
pmat=cycle2perm([n_osc,basedim],{[1,2,4,3],4});
prob = po_add_func(prob, '', 'symmetry', @symmetry, @d_symmetry, ...
  struct('i', 1:2, 'tau', 1/4, 'gamma', ...
  inv(pmat)),{'sym1' 'sym2'}, 'inactive');
prob = coco_set_parival(prob, {'sym1' 'sym2'}, zeros(2,1));
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
prob = coco_add_func(prob, 'amplitude', ...
  @(p,d,u) deal(d, max(u(1:2:end))-min(u(1:2:end))), ...
  struct('cid', 'po.orb.coll'), 'singular', ...
  'amplitude', 'uidx', uidx(data.coll_seg.maps.xbp_idx), ...
  'remesh', @ampremesh);
%prob = po_construct_sbtest(prob, 'po', 5);
%prob = coco_add_event(prob, 'UST', 'po.test.USTAB', (0:5-1)'+0.5);
prob = coco_set(prob, 'cont', 'PtMX', [20 20], 'NPR', 0, 'NAdapt', 1, 'norm', inf, 'h_min', 1e-5);
disp('RunSB: track symmetry breaking in 2 parameters')
coco(prob, 'runSB', [], 1, {'lambda' 'B' 'po.period' 'amplitude'}, [1 30]);
