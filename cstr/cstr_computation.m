%% Queries
%
% Demos corresponding to Section 2.1 in the Getting Fancy tutorial
% for COCO. For more information about the 'ep' and 'po' toolboxes, see
% EP-Tutorial.pdf and PO-Tutorial.pdf in the coco/help folder. 

% Copyright (C) Harry Dankowicz, Jan Sieber

%% Bifurcation analysis of equilibria 

if exist('sys_cstr', 'file')~=2
  syms x y beta gamma delta sigma
  F = sco_sym2funcs([(1-x)*exp(y/(1+beta*y))-x/delta; ...
    ((1-x)*exp(y/(1+beta*y))-y/sigma)/gamma], ...
    {[x; y], [beta; gamma; delta; sigma]}, ...
    {'x', 'p'}, 'filename', 'sys_cstr', 'maxorder', 3);
else
  F = sco_gen(@sys_cstr);
end

funcs = struct('f', F(''), 'dfdx', F('x'), 'dfdp', F('p'), ...
  'dfdx_dx', F('x*v'), 'dfdxdx_dx', F({'x','x*v'}), ...
  'dfdpdx_dx', F({'p','x*v'}), 'Dfdxdx', F({'x*v','x*v'}), ...
  'Dfdxdxdx', F({'x*v','x*v','x*v'}));
%% Two-parameter run but y fixed, equilibria
prob = coco_prob;
prob = ode_isol2ep(prob, '', funcs, [0.9; 2], ...
  {'beta', 'gamma', 'delta', 'sigma'}, ...
  [0; 1/10; 9/exp(2); 20/exp(2)]);
[data, uidx] = coco_get_func_data(prob, 'ep', 'data', 'uidx');
prob = coco_add_pars(prob, 'pars', uidx(data.ep_eqn.x_idx), {'x', 'y'});
prob = coco_set(prob, 'cont', 'NPR', 0);
coco(prob, 'ep_run', [], 1, {'x' 'sigma' 'delta'}, [0.2 0.9]);
%% Continuation of saddle node
SN = coco_bd_labs('ep_run', 'SN');
prob = coco_prob;
prob = ode_SN2SN(prob, '', 'ep_run', SN);
[data, uidx] = coco_get_func_data(prob, 'ep', 'data', 'uidx');
prob = coco_add_pars(prob, 'pars', uidx(data.ep_eqn.x_idx), {'x', 'y'});
prob = coco_set(prob, 'cont', 'NPR', 0, 'PtMX', 500, 'norm', inf,'h_max',0.001);
coco(prob, 'sn_run', [], {'sigma' 'delta' 'x' 'y' 'gamma' 'beta'}, {[0 20], [], [0.01 0.9]});
%% two-parameter bifurcation of Hopf bifurcation
HB = coco_bd_labs('ep_run', 'HB');
prob = coco_prob;
prob = ode_HB2HB(prob, '', 'ep_run', HB);
[data, uidx] = coco_get_func_data(prob, 'ep', 'data', 'uidx');
prob = coco_add_pars(prob, 'pars', uidx(data.ep_eqn.x_idx), {'x', 'y'});
data = struct('dfdxhan', funcs.dfdx, 'Dfdxdxhan', funcs.Dfdxdx, ...
  'Dfdxdxdxhan', funcs.Dfdxdxdx, 'nanflag', 1);
prob = ep_HB_add_func(prob, '', 'lyap', @lyapunov, data, 'regular', 'L1');
prob = coco_add_event(prob, 'UZ', 'delta', [0.15,0.2,0.24,0.3]);
prob = coco_add_event(prob, 'DH', 'L1', 0);
prob = coco_add_event(prob, 'SP', 'delta', 1.1426e-01);
prob = coco_add_event(prob, 'HOM', 'sigma', 0.44);
prob = coco_set(prob, 'cont', 'NPR', 0, 'PtMX', 1000, 'h_max', 1000, 'norm', inf, 'MaxRes', 1);
coco(prob, 'hb_run1', [], 1, {'sigma' 'delta' 'x' 'y', 'L1'}, {[0 20], [], [0.01 0.9]});
%% surface of Hopf bifurcations in 3 parameters beta, delta, sigma
HB = coco_bd_labs('ep_run', 'HB');
prob = coco_prob;
prob = coco_set(prob, 'ep', 'BTP', 'off');
prob = ode_HB2HB(prob, '', 'ep_run', HB);
[data, uidx] = coco_get_func_data(prob, 'ep', 'data', 'uidx');
prob = coco_add_pars(prob, 'pars', uidx(data.ep_eqn.x_idx), {'x', 'y'});
prob = coco_set(prob, 'cont', 'NPR', 0, 'PtMX', 6000, 'R', 0.1, 'R_max', 10, 'R_min', 0.001, 'MaxRes', 1, 'almax', 15);
coco(prob, 'hb_run2', [], 2, {'sigma' 'delta' 'x' 'y' 'beta' 'gamma' 'L1'}, {[0 20],[],[0.2 0.9],[],[0,1/4]});
%% provisional plot of Hopf surface in 3d
figure(1)
clf
atlas = coco_bd_read('hb_run2', 'atlas');
plot_atlas_kd(atlas.charts, 1,3,5)
hold on
thm =  struct('lspec', {{{'r', 'LineWidth', 3},{'r-.', 'LineWidth', 3}}},'ustab','');
coco_plot_bd(thm, 'hb_run1', 'sigma', 'x', 'beta')
grid
axis([0 inf 0.2 0.9 0 inf])
hold off
%%
figure(1)
clf
hold on
grid on
thm = struct('lspec', {{{'b--', 'LineWidth', 3},{'b-.', 'LineWidth', 3}}});
coco_plot_bd(thm, 'sn_run', 'x', {'x', 'y'}, @(x,y) x.*exp(-y)./(1-x), ...
{'x', 'y'}, @(x,y) y.*exp(-y)./(1-x))
coco_plot_bd(thm, 'hb_run1', 'x', {'x', 'y'}, @(x,y) x.*exp(-y)./(1-x), ...
{'x', 'y'}, @(x,y) y.*exp(-y)./(1-x))
thm = struct('lspec', {{{'r', 'LineWidth', 1},{'r', 'LineWidth', 1}}});
coco_plot_bd(thm, 'sn_run', 'x', 'delta', 'sigma')
coco_plot_bd(thm, 'hb_run1', 'x', 'delta', 'sigma')
%% track degenerate Hopf bifurcation in 3 parameters sigma-delta-beta (incl in 2d surface)
DH = coco_bd_labs('hb_run1', 'DH');
prob = coco_prob;
prob = ode_HB2HB(prob, '', 'hb_run1', DH(1));
[data, uidx] = coco_get_func_data(prob, 'ep', 'data', 'uidx');
prob = coco_add_pars(prob, 'pars', uidx(data.ep_eqn.x_idx), {'x', 'y'});
data = struct('dfdxhan', funcs.dfdx, 'Dfdxdxhan', funcs.Dfdxdx, ...
  'Dfdxdxdxhan', funcs.Dfdxdxdx, 'nanflag', 0);
prob = ep_HB_add_func(prob, '', 'lyap', @lyapunov, data, 'inactive', 'L1');
prob = coco_set_parival(prob, 'L1', 0);
prob = coco_set(prob, 'cont', 'NPR', 10, 'PtMX', 200, 'h_max', 0.4, 'norm', inf, 'al_max', 30,'MaxRes',Inf);
coco(prob, 'dh_run1', [], {'sigma' 'delta' 'beta' 'x' 'y' 'gamma' 'L1'}, {[0 1],[0,0.3],[0,0.06]});
%%
prob = coco_prob;
prob = ode_HB2HB(prob, '', 'hb_run1', DH(2));
[data, uidx] = coco_get_func_data(prob, 'ep', 'data', 'uidx');
prob = coco_add_pars(prob, 'pars', uidx(data.ep_eqn.x_idx), {'x', 'y'});
data = struct('dfdxhan', funcs.dfdx, 'Dfdxdxhan', funcs.Dfdxdx, ...
  'Dfdxdxdxhan', funcs.Dfdxdxdx, 'nanflag', 0);
prob = ep_HB_add_func(prob, '', 'lyap', @lyapunov, data, 'inactive', 'L1');
prob = coco_set_parival(prob, 'L1', 0);
prob = coco_set(prob, 'cont', 'NPR', 10, 'PtMX', 200, 'h_max', 100, 'norm', inf, 'al_max', 30,'MaxRes',Inf);
coco(prob, 'dh_run2', [],  {'sigma' 'delta' 'beta' 'x' 'y' 'gamma' 'L1'}, {[0 1],[0,0.3],[-1e-4,0.06]})
%% track degenerate Hopf bifurcation in 3 parameters sigma-delta-gamma
DH = coco_bd_labs('hb_run1', 'DH');
prob = coco_prob;
prob = ode_HB2HB(prob, '', 'hb_run1', DH(1));
[data, uidx] = coco_get_func_data(prob, 'ep', 'data', 'uidx');
prob = coco_add_pars(prob, 'pars', uidx(data.ep_eqn.x_idx), {'x', 'y'});
data = struct('dfdxhan', funcs.dfdx, 'Dfdxdxhan', funcs.Dfdxdx, ...
  'Dfdxdxdxhan', funcs.Dfdxdxdx, 'nanflag', 0);
prob = ep_HB_add_func(prob, '', 'lyap', @lyapunov, data, 'inactive', 'L1');
prob = coco_set_parival(prob, 'L1', 0);
prob = coco_set(prob, 'cont', 'NPR', 1000, 'PtMX', 200, 'h_max', 0.4, 'norm', inf, 'al_max', 30,'MaxRes',Inf);
coco(prob, 'dh_run3', [], {'sigma' 'delta' 'gamma' 'x' 'y' 'beta' 'L1'}, {[0 1],[0,0.3],[0,1/8]});
%% track degenerate Hopf bifurcation in 3 parameters sigma-delta-gamma
prob=coco_prob();
prob = ode_HB2HB(prob, '', 'hb_run1', DH(2));
[data, uidx] = coco_get_func_data(prob, 'ep', 'data', 'uidx');
prob = coco_add_pars(prob, 'pars', uidx(data.ep_eqn.x_idx), {'x', 'y'});
data = struct('dfdxhan', funcs.dfdx, 'Dfdxdxhan', funcs.Dfdxdx, ...
  'Dfdxdxdxhan', funcs.Dfdxdxdx, 'nanflag', 0);
prob = ep_HB_add_func(prob, '', 'lyap', @lyapunov, data, 'inactive', 'L1');
prob = coco_set_parival(prob, 'L1', 0);
prob = coco_set(prob, 'cont', 'NPR', 1000, 'PtMX', 200, 'h_max', 10, 'norm', inf, 'al_max', 30,'MaxRes',Inf);
coco(prob, 'dh_run4', [], {'sigma' 'delta' 'gamma' 'x' 'y' 'beta' 'L1'}, {[0 1],[0,0.3],[0,(7-3*sqrt(5))/2]});
%% Bifurcation analysis of periodic orbits
% run to SNIC
SP = coco_bd_labs('hb_run1', 'SP');
prob = coco_prob;
prob = coco_set(prob, 'coll', 'NTST', 100);
prob = coco_set(prob, 'po', 'bifus', 'off');
prob = ode_HB2po(prob, '', 'hb_run1', SP(2));
prob = po_add_func(prob, '', 'amplitude', ...
    @(~,xbp,~,~,~) max(xbp(1:2:end))-min(xbp(1:2:end)), 'amplitude', 'regular');
prob = coco_set(prob, 'cont', 'NPR', 0, 'NAdapt', 1, 'PtMX', [1000 0], 'norm', inf, 'FP', 'on');
prob=coco_add_event(prob, 'UZ', 'po.period', (0:5:20)+1e-6);
coco(prob, 'po_run', [], 1, {'po.period' 'sigma' 'amplitude'}, {[0 20], [0 10]});
SP = coco_bd_labs('hb_run1', 'SP');
%% PO run to homoclinic
HOM=coco_bd_labs('hb_run1','HOM');
prob = coco_prob;
prob = coco_set(prob, 'coll', 'NTST', 100);
prob = coco_set(prob, 'po', 'bifus', 'off');
prob = ode_HB2po(prob, '', 'hb_run1', HOM(2));
prob = po_add_func(prob, '', 'amplitude', ...
    @(~,xbp,~,~,~) max(xbp(1:2:end))-min(xbp(1:2:end)), 'amplitude', 'regular');
prob = coco_set(prob, 'cont', 'NPR', 0, 'NAdapt', 1, 'PtMX', [1000 0], 'norm', inf, 'FP', 'on');
prob=coco_add_event(prob, 'UZ', 'po.period', (0:5:20)+1e-6);
coco(prob, 'po_run2', [], 1, {'po.period' 'sigma' 'amplitude'}, {[0 20], [0 10]});

%% Hopf bubbles
UZH = coco_bd_labs('hb_run1', 'UZ');
for i=1:4
    prob=coco_prob();
    prob = coco_set(prob, 'coll', 'NTST', 100);
    prob = ode_HB2po(prob, '', 'hb_run1', UZH(i));
    prob = po_add_func(prob, '', 'amplitude', ...
        @(~,xbp,~,~,~) max(xbp(1:2:end))-min(xbp(1:2:end)), 'amplitude', 'regular');
    prob = coco_set(prob, 'cont', 'NPR', 0, 'NAdapt', 1, 'PtMX', [100 0], 'norm', inf, 'FP', 'on');
    coco(prob, sprintf('bubble%d',i), [], 1, {'sigma','po.period', 'amplitude'});
end
%% provisional plot
figure(3);clf;ax3=gca;hold(ax3,'on');
for i=1:4
    coco_plot_bd(sprintf('bubble%d',i),'sigma','amplitude');
end
%% track SN of periodic orbits
SN = coco_bd_labs('bubble1', 'SN');
prob = coco_prob;
prob = coco_set(prob, 'coll', 'NTST', 100);
prob = coco_set(prob, 'po', 'bifus', 'off');
prob = ode_po2SN(prob, '', 'bubble1', SN);
prob = po_add_func(prob, '', 'amplitude', ...
    @(~,xbp,~,~,~) max(xbp(1:2:end))-min(xbp(1:2:end)), 'amplitude', 'regular');
%prob = coco_set(prob, 'corr', 'TOL', 1e-4, 'ResTOL', 1e-6);
prob = coco_set(prob, 'cont', 'NPR', 10, 'NAdapt', 1, 'PtMX', [1000,50],'h_max',100, 'norm', inf, 'FP', 'on');
coco(prob, 'po_sn_run', [], 1, {'sigma' 'delta' 'po.period' 'amplitude'}, {[0.4,0.8],[00.1,0.3],[0,20],[0 10]});

%% homoclinic

EP = coco_bd_labs('po_run', 'EP');
[sol, data] = coll_read_solution('po.orb', 'po_run', max(EP)); % Periodic orbit with T=1.8652e+01
f = feval(funcs.f, sol.xbp', repmat(sol.p, [1 size(sol.xbp, 1)])); % Evaluate vector field at basepoints
[~, idx] = min(sqrt(sum(f.*f, 1))); % Find basepoint closest to equilibrium

scale = 25;
T  = sol.T;
t0 = [sol.tbp(1:idx,1) ; T*(scale-1)+sol.tbp(idx+1:end,1)]; % Crank up period by factor scale
x0 = sol.xbp;
p0 = sol.p;
pnames = {'beta', 'gamma', 'delta', 'sigma'};

% Initialize continuation problem structure with the same number of
% intervals as in previous run.

prob = coco_prob;
prob = coco_set(prob, 'coll', 'NTST', data.coll.NTST);
prob = coco_set(prob, 'po', 'bifus', 'off');
prob = ode_isol2po(prob, '', funcs, t0, x0, pnames, p0);
prob = coco_set(prob, 'cont', 'NAdapt', 10);
prob = coco_xchg_pars(prob, 'delta', 'po.period');

coco(prob, 'po_long_find', [], 0, {'sigma' 'po.period'});

prob = coco_prob;
prob = coco_set(prob, 'po', 'bifus', 'off');
prob = coco_set(prob, 'coll', 'NTST', 100);
prob = ode_po2po(prob, '', 'po_long_find', 2);
prob = coco_xchg_pars(prob, 'delta', 'po.period');
prob = po_add_func(prob, '', 'amplitude', ...
  @(~,xbp,~,~,~) max(xbp(1:2:end))-min(xbp(1:2:end)), 'amplitude', 'regular');
prob = po_add_func(prob, '', 'saddle_q', ...
  @(data,xbp,T0,T,p)saddle_q(data,xbp,T0,T,p,funcs), {'xeq','yeq','fmin','det','tr'}, 'regular');
prob=coco_add_event(prob, 'NSA', 'tr', 0);
prob=coco_add_event(prob, 'NCS', 'det', -0.05);
prob = coco_set(prob, 'cont', 'NPR', 0, 'NAdapt', 1, 'PtMX', [192 102], 'norm', inf, 'h_max', 10, 'h_min', 1e-4);

coco(prob, 'po_long_run', [], 1, {'sigma' 'delta' 'amplitude' 'det' 'tr' 'xeq' 'yeq' 'po.period'}, [0 20]);
%%

%%
figure(5)
hold on
thm.ustab = '';
thm.zlab = '\delta';
thm.special = {};
coco_plot_bd(thm, 'po_long_run', 'sigma', 'MAX(x)', 'delta')
coco_plot_bd(thm, 'po_long_run', 'sigma', 'MIN(x)', 'delta')
hold off
%%
figure(1)
clf
hold on
 thm =  struct('lspec', {{{'k.-', 'LineWidth', 1},{'b-.', 'LineWidth', 3}}},'ustab','');
coco_plot_bd(thm,'po_long_run', 'sigma', 'delta', 'MAX(x)')
coco_plot_bd('sn_run', 'sigma', 'delta', 'MAX(x)')
coco_plot_bd('hb_run1', 'sigma', 'delta', 'MAX(x)')
hold off
%%
figure(1)
clf
hold on
 thm =  struct('lspec', {{{'b.', 'MarkerSize', 12},{'b-.', 'LineWidth', 3}}},'ustab','');
coco_plot_bd(thm,'po_long_run', 'sigma', 'delta')
coco_plot_bd('sn_run', 'sigma', 'delta')
coco_plot_bd('hb_run1', 'sigma', 'delta')
hold off
