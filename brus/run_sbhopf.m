function run_symSB(figs)

if nargin<1
  figs = false;
end

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
    {'x', 'p'}, 'filename', 'sys_brus'); %#ok<NODEF>
  funcs = struct('f', F(''), 'dfdx', F('x'), 'dfdp', F('p'));
else
  F = sco_gen(@sys_brus);
  funcs = struct('f', F(''), 'dfdx', F('x'), 'dfdp', F('p'));
end

%% continuation along branch of equilibria

A = 2;
B = 5.9;

x0 = [A; B/A; A; B/A; A; B/A; A; B/A];
pnames = {'A', 'B', 'lambda'};
p0 = [A; B; 1];

prob = coco_prob;
prob = coco_set(prob, 'ep', 'SN', 'off', 'HB', 'off');
prob = ode_isol2ep(prob, '', funcs, x0, pnames, p0);

prob = coco_add_func(prob, 'ustab', @ep_USTAB, [], ...
  'discrete', 'ep.USTAB', 'passChart', 'requires', 'ep.test');
prob = coco_add_event(prob, 'UST', 'ep.USTAB');

prob = coco_set(prob, 'cont', 'NPR', inf, 'BP', 'off');
coco(prob, 'run_eq', [], 1, 'lambda', [1 30]);

%% run_po1, branch off from equivariant Hopf bifurcation

prob = coco_prob;
prob = osc_equivHB2po(prob, '', 'run_eq', 2, {[1,2,4,3],4});

prob = coco_add_func(prob, 'ustab', @po_USTAB, [], 'discrete', ...
  'po.USTAB', 'passChart', 'requires', {'po.orb.coll.test' 'po.test'});
tst.poid     = 'po';
tst.ustnames = strcat('UST', arrayfun(@(i){num2str(i)},1:8));
prob = coco_add_chart_data(prob, 'po.sb', [], []);
prob = coco_add_event(prob, @evhan_SB, tst, 'po.USTAB');

prob = coco_set(prob, 'cont', 'PtMX', [100 0], 'NPR', inf,...
  'NAdapt', 1, 'norm', inf);
disp('run_po1: branch off from Equiv. Hopf')
coco(prob, 'run_po1', [], 1, {'lambda' 'po.period' 'amplitude'}, [1 30]);

%% Track symmetry breaking

UST  = coco_bd_labs('run_po1', 'UST1');
prob = coco_prob;
prob = coco_set(prob, 'po', 'SN', 'off', 'PD', 'off', 'TR', 'off');

sol  = po_read_solution('', 'run_po1', UST(2));
chart = coco_read_solution('', 'run_po1', UST(2), 'chart');
cdata = coco_get_chart_data(chart, 'po.sb');
sol.var.v = cdata.sb.v;
sol.sn.u0 = -cdata.sb.b;
sol.sn.t0 = [];
prob = ode_coll2coll(prob, 'po.orb', 'run_po1', 'po.orb', UST(2), ...
  '-var', sol.var.v);

data = coco_read_solution('po', 'run_po1', UST(2), 'data');
data = po_init_data(prob, data, '', 'ode');
[prob, data] = po_add(prob, data, '-no-test');
prob = po_add_SN(prob, data, sol);
prob = ode_add_tb_info(prob, '', 'po', 'po', 'po', po_sol_info('SB'));

pmat = cycle2perm([4,2],{[1,2,4,3],4});
prob = po_add_func(prob, '', 'symmetry', @symmetry, @d_symmetry, ...
  struct('i', 1:2, 'tau', 1/4, 'gamma', inv(pmat)), ...
  {'sym1' 'sym2'}, 'inactive');
prob = coco_set_parival(prob, {'sym1' 'sym2'}, zeros(2,1));

[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
prob = coco_add_func(prob, 'amplitude', ...
  @(p,d,u) deal(d, max(u(1:2:end))-min(u(1:2:end))), ...
  struct('cid', 'po.orb.coll'), 'singular', ...
  'amplitude', 'uidx', uidx(data.coll_seg.maps.xbp_idx), ...
  'remesh', @ampremesh);

prob = coco_set(prob, 'cont', 'NPR', inf, 'PtMX', [20 20], ...
  'NAdapt', 1, 'norm', inf);
disp('run_SB: track symmetry breaking in 2 parameters')
coco(prob, 'run_SB', [], 1, {'lambda' 'B' 'po.period' 'amplitude'}, [15.5 18.5]);

%% Track equivariant Hopf bifurcation

prob = coco_prob;
prob = coco_set(prob, 'ep', 'SN', 'off', 'HB', 'off');

oid = '';
cycle = {[1], 1; [2,3,4],1};
run = 'run_eq';
lab = 2;

[sol, data] = ep_read_solution(oid, run, lab);
las = sol.ep_test.la;
omega = imag(las(1));

funcs = struct('f', data.fhan);
fnames = fieldnames(data);
for i=1:length(fnames)
  if strncmp(fnames(i), 'dfd', 3)
    funcs.(fnames{i}(1:end-3)) = data.(fnames{i});
  end
end
J = data.dfdxhan(sol.x,sol.p);

n_osc = 0;
for i=1:size(cycle,1)
  n_osc = n_osc + numel(cycle{i,1});
end
basedim = length(sol.x)/n_osc;
[pmat, irot] = cycle2perm([n_osc, basedim], cycle);
Jpl = [J-1j*omega*eye(8); pmat-diag(exp(2*pi*1j./irot))];
[a,~] = svds(Jpl',1,'smallest');
v = real(a);
w = imag(a);

k = omega^2;
al = [-w'*v v'*v];
al = al/norm(al);
nv = al(1)*v + al(2)*w;
nv = nv/norm(nv);

sol.var.v  = [ v/norm(v) -sqrt(k)*w/norm(v) ];
sol.var.w  = [ -sqrt(k)*w/norm(v) -k*v/norm(v) ];
sol.var.u0 = [ sol.var.v(:) ; sol.var.w(:) ];
sol.var.t0 = [];

sol.hb.k  = k;
sol.hb.nv = nv;
sol.hb.u0 = k;
sol.hb.t0 = [];

data = ode_init_data(prob, data, oid, 'ep');
[prob, data] = ep_add(prob, data, sol, '-no-test', '-cache-jac');
[prob, data] = ep_add_var(prob, data, sol);
prob = ep_add_HB(prob, data, sol);
prob = ode_add_tb_info(prob, oid, 'ep', 'ep', 'ep', ep_sol_info('equivHB'));
prob = ep_HB_add_func(prob,'','sym',@v_symmetry,...
  struct('gamma', pmat, 'tau', irot, 'dfdxhan', funcs.dfdx), 'inactive',...
  arrayfun(@(i){sprintf('sym%02d',i)},1:2*size(pmat,1)));

[data, uidx] = coco_get_func_data(prob, 'ep.HB', 'data', 'uidx');
prob = coco_add_pars(prob, 'kmon', uidx(data.pr.ep_hb.k_idx), 'k');

prob = coco_set(prob, 'cont', 'PtMX', [20 0], 'NPR', 100, 'norm', inf);
coco(prob, 'runHB', [], 1, {'lambda','B','k'}, [15 30]);

%% graphical representation

if figs
  clear
  ty  = @(r,s)strcmp(r{:,'TYPE'},s);
  alw = {'linewidth',2};
  lw  = {'linewidth',3};
  txt = {'FontSize',20,'FontName','courier','FontWeight','bold'};
  ltx = {'Interpreter','latex'};
  runSB = coco_bd_table('run_SB');
  runHB = coco_bd_table('runHB');
  EPstart = runSB(ty(runSB,'EP')&runSB.PT==0,:);
  startlab = EPstart.LAB{1};
  
  figure(3); tl = tiledlayout(1,2); tl.TileSpacing='tight'; clr=lines();
  nexttile(); ax1 = gca; hold(ax1,'on');
  plsb = plot(ax1,runSB.lambda,runSB.B,'Color',clr(1,:),'DisplayName',...
    sprintf('Symmetry breaking:\n$(1\\,2\\,4\\,3)_4\\mapsto(1\\,4)_2(2\\,3)_2$'),lw{:});
  plhb = plot(ax1,runHB.lambda,runHB.B,'Color',clr(2,:),'DisplayName',...
    sprintf('Equivariant Hopf\nbifurcation'),lw{:});
  GHB = runSB(abs(runSB.amplitude)<1e-4,:);
  tmarker = {'kd', 'MarkerFaceColor', clr(2,:), alw{:}, 'MarkerSize',10};
  plep = plot(ax1, EPstart.lambda(1), EPstart.B(1), tmarker{:}, ...
    'DisplayName', sprintf('Starting point'), alw{:});
  plghb = plot(ax1, GHB.lambda(1), GHB.B(1), 'ks', 'MarkerFaceColor', ...
    clr(3,:),'DisplayName', ...
    sprintf('Degenerate equivariant\nHopf bifurcation'), alw{:}, ...
    'MarkerSize',10);
  legend(ax1, [plsb, plhb, plep, plghb], txt{:}, ltx{:}, ...
    'Location', 'NorthEast');
  ax1.XLim = [min(runSB.lambda),max(runSB.lambda)*1.05];
  ylim(ax1,[5.5,7.5]);
  set(ax1, txt{:}, 'box', 'on', alw{:});
  xlabel(ax1,'$\epsilon\times10^{3}$', ltx{:})
  ylabel(ax1,'$B$', ltx{:})
  ax1.XTick = 15:18;
  ax1.YTick = 5:7;
  ax1.LabelFontSizeMultiplier = 1.2;
  
  [sol, data] = coll_read_solution('po.orb', 'run_SB', startlab);
  nexttile(); ax2 = gca; hold(ax2,'on');
  tpl = plot(ax2, sol.tbp/sol.T, sol.xbp(:,1:2:8), lw{:});
  for k=1:4
    set(tpl(k),'DisplayName',sprintf('$x_%d$',k));
  end
  pmark=plot(ax2,NaN,NaN,tmarker{:},'DisplayName','');
  dummy=plot(ax2,NaN,NaN,'wo','DisplayName','');
  ax2.ColorOrderIndex=1;
  vbp = reshape(sol.var.vbp,[8 numel(data.coll_seg.maps.tbp_idx)])';
  vpl=plot(ax2,sol.tbp/sol.T,vbp(:,1:2:8)*0.95,'-.',lw{:});
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
end

end

%% utilities

function [prob, status, xtr] = ampremesh(prob, data, chart, old_u, old_V) %#ok<INUSD>

[fdata, uidx] = coco_get_func_data(prob, data.cid, 'data', 'uidx');
seg  = fdata.coll_seg;
maps = seg.maps;
data.coll_seg = seg;

xtr       = [];
prob      = coco_change_func(prob, data, 'uidx', uidx(maps.xbp_idx));
status    = 'success';

end

function prob = osc_equivHB2po(prob, oid, run, lab, cycle)
% construct periodic orbit problem with given spatiotemporal symmetry from
% equivariant Hopf bifurcation, but do not impose the symmetry to eliminate
% branch points.

[sol, data] = ep_read_solution(oid, run, lab);
las = sol.ep_test.la;
omega = imag(las(1));

funcs = struct('f', data.fhan);
fnames = fieldnames(data);
for i=1:length(fnames)
  if strncmp(fnames(i), 'dfd', 3)
    funcs.(fnames{i}(1:end-3)) = data.(fnames{i});
  end
end
J = data.dfdxhan(sol.x,sol.p);

n_osc = 0;
for i=1:size(cycle,1)
  n_osc = n_osc + numel(cycle{i,1});
end
basedim = length(sol.x)/n_osc;
[pmat, irot] = cycle2perm([n_osc, basedim], cycle);
Jpl = [J-1j*omega*eye(8); pmat-diag(exp(2*pi*1j./irot))];
[a,~] = svds(Jpl',1,'smallest');
re = real(a);
im = imag(a);

period = 2*pi/omega;
t0 = 0:period/400:period;
x0 = sol.x + 1e-4*(re.*cos(t0*2*pi/period)-im.*sin(t0*2*pi/period));

prob = coco_set(prob, 'po', 'SN', 'off', 'PD', 'off', 'TR', 'off');
prob = coco_set(prob, 'cont', 'Valpha', 0);
prob = coco_set(prob, 'coll', 'NTST', 20);
prob = ode_isol2po(prob, '', funcs, t0', x0', data.pnames, sol.p);
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
prob = coco_add_func(prob, 'amplitude', ...
  @(p,d,u) deal(d, max(u(1:2:end))-min(u(1:2:end))), ...
  struct('cid', 'po.orb.coll'), 'singular', ...
  'amplitude', 'uidx', uidx(data.coll_seg.maps.xbp_idx), ...
  'remesh', @ampremesh);

% for i=1:size(cycle,1)
%   [c,irot] = deal(cycle{i,:});
%   name = ['cycle',sprintf('%d',c)];
%   pmat = cycle2perm([n_osc,basedim],{c,irot});
%   inds = sort(sub2ind([basedim,n_osc], repmat(1:basedim,1,length(c)), ...
%     reshape(repmat(c,basedim,1),1,[])));
%   prob = po_add_func(prob, '', name, @symmetry, @d_symmetry, ...
%     struct('i', inds, 'tau', 1/irot, 'gamma',inv(pmat)), ...
%     arrayfun(@(j){sprintf('sym%01d_%01d',i,j)},1:length(inds)), ...
%     'inactive');
% end

end

function [data, chart, y] = ep_USTAB(~, data, chart, ~)
%EP_USTAB   Monitor functions for number of unstable eigenvalues

cdata = coco_get_chart_data(chart, 'ep.test');
if ~isempty(cdata) && isfield(cdata, 'la')
  la = cdata.la;
end
y = sum(real(la)>=0);

end

function [data, chart, y] = po_USTAB(~, data, chart, ~)
%PO_USTAB   Monitor functions for number of unstable Floquet multipliers

cdata = coco_get_chart_data(chart, 'po.test');
if ~isempty(cdata) && isfield(cdata, 'la')
  la = cdata.la;
end
y = sum(abs(la)>=1);

end

function [data, cseg, msg] = evhan_SB(prob, data, cseg, cmd, msg)
%EVHAN_SB  symmetry breaking event handler: store nullvectors

pid = data.poid;
tid = coco_get_id(pid, 'test');
cid = coco_get_id(pid, 'orb.coll');
fid = coco_get_id(pid, 'sb');
switch cmd
  case 'init'
    if isfield(msg, 'finish') || strcmp(msg.action, 'warn')
      msg.action = 'finish';
    elseif strcmp(msg.action, 'locate')
      msg.action = 'warn';
    else
      cdata = coco_get_chart_data(cseg.ptlist{1}, tid);
      la0   = cdata.la;
      cdata = coco_get_chart_data(cseg.ptlist{end}, tid);
      la1   = cdata.la;
      data.degree_change = abs(sum(abs(la1)>1)-sum(abs(la0)>1));
      msg.point_type = data.ustnames{data.degree_change};
      msg.action = 'locate';
      msg.idx    = 1;
    end
  case 'check'
    cdata = coco_get_chart_data(cseg.curr_chart, fid);
    
    [fdata, uidx] = coco_get_func_data(prob, cid, 'data', 'uidx');
    pdata = coco_get_func_data(prob, pid, 'data');
    u     = cseg.curr_chart.x(uidx);
    maps  = fdata.coll_seg.maps;
    f0    = fdata.ode_F(fdata, 0, u(maps.x0_idx), u(maps.p_idx));
    M     = [ pdata.po_M-eye(numel(f0)) f0 ; f0' 0];
    [V, D]   = eig(M);
    [~, idx] = min(abs(diag(D)));
    V     = V(:,idx); % generalized eigenvector for double eigenvalue at 1
    V     = V/norm(V(1:end-1));
    sb.v = V(1:end-1);
    sb.b = V(end);
    sb.M = pdata.po_M;
    
    cdata.sb   = sb;
    cseg.curr_chart = coco_set_chart_data(cseg.curr_chart, fid, cdata);
    
    msg.action = 'add';
    msg.finish = true;
end

end
