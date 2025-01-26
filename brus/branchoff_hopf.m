function branchoff_hopf(figs)
%% test branching off from Hopf bifurcation

if nargin<1
  figs = true;
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
    {'x', 'p'}, 'filename', 'sys_brus');
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

%% common settings for continuation of periodic orbits

parbd = {{'lambda', 'po.period' 'amplitude'}, {[1 30], [0,100], [0,1]}};

%% continuation along branch P1

prob = coco_prob;
prob = osc_equivHB2po(prob, '', 'run_eq', 2, {[1,2,4,3],4});
prob = coco_add_func(prob, 'ustab', @po_USTAB, [], 'discrete', ...
  'po.USTAB', 'passChart', 'requires', {'po.orb.coll.test' 'po.test'});
prob = coco_add_event(prob, 'UST', 'po.USTAB');

prob = coco_set(prob, 'cont', 'PtMX', [100 0], 'NPR', inf,...
  'NAdapt', 1, 'norm', inf);
coco(prob, 'hopf_eq1', [], 1, parbd{:});

%% continuation along branch P2

prob = coco_prob;
prob = osc_equivHB2po(prob, '', 'run_eq', 2, {[2,1],2; [3,4],1});
prob = coco_add_func(prob, 'ustab', @po_USTAB, [], 'discrete', ...
  'po.USTAB', 'passChart', 'requires', {'po.orb.coll.test' 'po.test'});
prob = coco_add_event(prob, 'UST', 'po.USTAB');

prob = coco_set(prob, 'cont', 'PtMX', [100 0], 'NPR', inf, ...
  'NAdapt', 1, 'norm', inf);
coco(prob, 'hopf_eq2', [], 1, parbd{:});

%% continuation along branch P3

prob = coco_prob;
prob = osc_equivHB2po(prob, '', 'run_eq', 2, {1,1; [2,3,4],1});
prob = coco_add_func(prob, 'ustab', @po_USTAB, [], 'discrete', ...
  'po.USTAB', 'passChart', 'requires', {'po.orb.coll.test' 'po.test'});
prob = coco_add_event(prob, 'UST', 'po.USTAB');

prob = coco_set(prob, 'cont', 'PtMX', [90 0], 'NPR', inf, ...
  'NAdapt', 1, 'norm', inf);
coco(prob, 'hopf_eq3', [], 1, parbd{:});

%% graphical representation

if figs
  srun = @(i)sprintf('hopf_eq%d', i);
  hopf_eq = arrayfun(@(i){coco_bd_table(srun(i))}, 1:3);
  
  nruns = length(hopf_eq);
  names = { ...
    '$P_1$: $\Pi\sim(1\,2\,4\,3)_4$', ...
    '$P_2$: $\Pi\sim(2\,1)_2(3\,4)_1$', ...
    '$P_3$: $\Pi\sim(3\,4\,2)_1$'};
  
  proflab = {{'UST',2}, {'FP',1}, {'UST',1}};
  clri = [1,2,4];
  
  figure(1); clf;
  tl = tiledlayout(3,2);
  tl.TileSpacing = 'tight';
  nexttile([3,1]);
  ax1 = gca; hold(ax1,'on');
  
  lw  = {'linewidth', 2};
  txt = {'FontSize', 20, 'FontName', 'courier', 'FontWeight', 'bold'};
  ltx = {'Interpreter', 'latex'};
  clr = lines();
  ty  = @(r,s)strcmp(r{:,'TYPE'}, s);
  
  for i=1:nruns
    pl(i) = plot(ax1, hopf_eq{i}.lambda, hopf_eq{i}.amplitude,lw{:}, ...
      'DisplayName', names{i}, 'color', clr(clri(i),:)); %#ok<*AGROW>
    specsel{i} = ty(hopf_eq{i}, 'UST') | ty(hopf_eq{i}, 'BP') | ...
      ty(hopf_eq{i}, 'FP');
    uz = plot(ax1, hopf_eq{i}.lambda(specsel{i}), ...
      hopf_eq{i}.amplitude(specsel{i}), 'ks', lw{:});
  end
  
  hbi = find(ty(hopf_eq{2}, 'EP') & hopf_eq{2}.amplitude<1e-2);
  hb  = plot(ax1, hopf_eq{2}.lambda(hbi), hopf_eq{2}.amplitude(hbi), ...
    'kd', lw{:}, 'MarkerSize', 12, 'MarkerFaceColor', 'g', ...
    'DisplayName', 'equivariant Hopf bifurcation');
  uz.DisplayName = sprintf('UST (number of unstable$^{\\phantom{\\int}}$\nFloquet multipliers changes)');
  legend(ax1, [pl,uz,hb], ltx{:}, txt{:}, 'Location', 'North');
  xlim(ax1, [18,22]);
  ylim(ax1, [0,1]);
  xlabel(ax1, '$\epsilon\times10^{3}$', ltx{:}, txt{:})
  ylabel(ax1, 'amplitude', ltx{:}, txt{:})
  set(ax1, txt{:}, lw{:}, 'box', 'on');
  
  for i=1:nruns
    ustlabs{i} = coco_bd_labs(srun(i), proflab{i}{1});
    hopf_eqprof{i} = po_read_solution(srun(i), ustlabs{i}(proflab{i}{2}));
    nexttile(2*i); ax(i)=gca; hold(ax(i), 'on');
    t = hopf_eqprof{i}.tbp/hopf_eqprof{i}.T;
    tprofs{i} = plot(ax(i), t, hopf_eqprof{i}.xbp(:,1:2:end), lw{:});
    for k=1:4
      set(tprofs{i}(k), 'DisplayName', sprintf('$x_%d$', k));
    end
    pllab{i} = plot(ax(i), NaN, NaN, 'o', 'MarkerEdgeColor', 'k', ...
      lw{:}, 'MarkerFaceColor', clr(clri(i),:), ...
      'DisplayName', sprintf('$P_%d$', i));
    legend(ax(i), [pllab{i};tprofs{i}], txt{:}, ltx{:}, ...
      'Location','EastOutside');
    set(ax(i), txt{:}, lw{:}, 'box', 'on');
    if i<3
      ax(i).XTickLabel = {};
    end
    ax(i).YTickLabel = {};
    ax(i).YLimitMethod = 'padded';
    ax(i).XLimitMethod = 'padded';
  end
  set(tprofs{2}(4), 'LineStyle', '--');
  set(tprofs{3}(4), 'LineStyle', '--');
  set(tprofs{3}(2), 'LineStyle', ':', 'LineWidth', 4);
  xlabel(ax(3), 'time $t/T$', txt{:}, ltx{:});

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
% equivariant Hopf bifurcation and impose the symmetry to eliminate branch
% points.

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

for i=1:size(cycle,1)
  [c,irot] = deal(cycle{i,:});
  name = ['cycle',sprintf('%d',c)];
  pmat = cycle2perm([n_osc,basedim],{c,irot});
  inds = sort(sub2ind([basedim,n_osc], repmat(1:basedim,1,length(c)), ...
    reshape(repmat(c,basedim,1),1,[])));
  prob = po_add_func(prob, '', name, @symmetry, @d_symmetry, ...
    struct('i', inds, 'tau', 1/irot, 'gamma',inv(pmat)), ...
    arrayfun(@(j){sprintf('sym%01d_%01d',i,j)},1:length(inds)), ...
    'inactive');
end

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
