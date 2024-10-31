function prob=hopf_equiv_po(prob,oid,omega,n_osc,cycles,epsilon)
period = 2*pi/omega;
t0 = 0:period/400:period;
[re, im,sol,data] = equiv_hopf_var(prob,oid,omega,n_osc,cycles);
alldim=length(sol.x);
x0 = sol.x + epsilon*...
  (re.*cos(t0*2*pi/period)-im.*sin(t0*2*pi/period));
prob = coco_prob;
prob = coco_set(prob, 'coll', 'MXCL', 'off');
prob = coco_set(prob, 'po', 'SN', 'off', 'PD', 'off', 'TR', 'off');
funcs={data.pr.fhan,data.pr.dfdxhan,data.pr.dfdphan};
prob = ode_isol2po(prob, '', funcs{:}, t0', x0', data.pr.pnames, sol.p);
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
prob = coco_add_func(prob, 'amplitude', ...
  @(p,d,u) deal(d, max(u(1:2:end))-min(u(1:2:end))), ...
  struct('cid', 'po.orb.coll'), 'singular', ...
  'amplitude', 'uidx', uidx(data.coll_seg.maps.xbp_idx), ...
  'remesh', @ampremesh);
basedim=alldim/n_osc;
if ~iscell(cycles)
    cycles={cycles,length(cycles)};
end
for i=1:size(cycles,1)
    [c,irot]=deal(cycles{i,:});
    name=['cycle',sprintf('%d',c)];
    pmat=cycle2perm([n_osc,basedim],{c,irot});
    inds=sort(sub2ind([basedim,n_osc],repmat(1:basedim,1,length(c)),...
        reshape(repmat(c,basedim,1),1,[])));
    prob = po_add_func(prob, '', name, @symmetry, @d_symmetry, ...
        struct('i', inds, 'tau', 1/irot, 'gamma',inv(pmat)), ...
        arrayfun(@(j){sprintf('sym%01d_%01d',i,j)},1:length(inds)), 'inactive');
end
%prob = po_construct_sbtest(prob, 'po', 5);
%prob = coco_add_event(prob, 'UST', 'po.test.USTAB', (0:5-1)'+0.5);
prob = coco_add_event(prob, 'UST', 'po.test.USTAB', (0:5-1)'+0.5);
end