function [prob, data] = po_construct_sbtest(prob, poid, maxdim)
%PO_CONSTRUCT_SBTST  check for symmetry breaking, works only for single
% branches branching off

data = init_data(poid, maxdim);
data = coco_func_data(data);
data.sh.po_M = [];
prob = coco_add_chart_data(prob, data.po_tst.fid, [], []);

void = coco_get_id(coco_get_id(poid, 'orb.coll'), 'test');

prob = coco_add_func(prob, data.po_tst.fid, @test_SB, data, ...
  'regular', data.po_tst.pids, 'requires', void, 'passChart');
for i=1:maxdim
  prob = coco_add_event(prob, @evhan_SB, data, 'SP', data.po_tst.pids{i}, 0);
end

end

%%
function data = init_data(poid, maxdim)

tst.ustnames = strcat('UST', arrayfun(@(i){num2str(i)},1:maxdim));
tst.maxdim   = maxdim;
tst.poid     = poid;
tst.fid      = coco_get_id(poid, 'sb');
tst.pids     = coco_get_id(tst.fid, tst.ustnames);

data.po_tst  = tst;

end

%%
function [data, chart, y] = test_SB(prob, data, chart, ~)
%TEST_SB  Monitor function for degree of stability change
%
% y(i) : number of unstable eigenvalues>=i-1

tst = data.po_tst;

cdata = coco_get_chart_data(chart, tst.fid);
if ~isempty(cdata) && isfield(cdata, 'la') % check this!
  la = cdata.la;
else
  fdata  = coco_get_func_data(prob, coco_get_id(tst.poid, 'orb.coll'), 'data');
  ctst   = fdata.coll_tst;
  M0     = ctst.M(ctst.M0_idx,:);
  M1     = ctst.M(ctst.M1_idx,:);
  M      = M1/M0;
  [D, la] = eig(full(M)); %#ok<ASGLU>
  la     = diag(la); % check not getting D!
%  if data.ode.autonomous
    [~, idx] = sort(abs(la-1));
    la = la(idx(2:end));
%  end
  cdata.la = la;
  chart  = coco_set_chart_data(chart, tst.fid, cdata);
  data.po_M = M;
end
% Stability degree
nunst = sum(abs(la)>=1);
y = nunst-(0:tst.maxdim-1)'-0.5;

end

%%
function [data, cseg, msg] = evhan_SB(prob, data, cseg, cmd, msg)
%EVHAN_SB  symmetry breaking event handler: store nullvectors

tst = data.po_tst;
fid = tst.fid;
switch cmd
  case 'init'
    if isfield(msg, 'finish') || strcmp(msg.action, 'warn')
      msg.action = 'finish';
    elseif strcmp(msg.action, 'locate')
      msg.action = 'warn';
    else
      cdata = coco_get_chart_data(cseg.ptlist{1}, fid);
      la0   = cdata.la;
      cdata = coco_get_chart_data(cseg.ptlist{end}, fid);
      la1   = cdata.la;
      tst.degree_change = abs(sum(abs(la1)>1)-sum(abs(la0)>1));
      msg.point_type = tst.ustnames{tst.degree_change};
      msg.action = 'locate';
      msg.idx    = 1;
    end
  case 'check'
    cdata = coco_get_chart_data(cseg.curr_chart, fid);

      [fdata, uidx] = coco_get_func_data(prob, coco_get_id(tst.poid, 'orb.coll'), 'data', 'uidx');
      u     = cseg.curr_chart.x(uidx);
      maps  = fdata.coll_seg.maps;
      f0    = fdata.ode_F(fdata, 0, u(maps.x0_idx), u(maps.p_idx));
      M     = [ data.po_M-eye(numel(f0)) f0 ; f0' 0];
      [V, D]   = eig(M);
      [~, idx] = min(abs(diag(D)));
      V     = V(:,idx); % generalized eigenvector for double eigenvalue at 1
      V     = V/norm(V(1:end-1));
      sn.v = V(1:end-1);
      sn.b = V(end);
      sn.M = data.po_M;

    cdata.sn   = sn;
    cseg.curr_chart = coco_set_chart_data(cseg.curr_chart, fid, cdata);

    msg.action = 'add';
    msg.finish = true;
end

data.po_tst = tst;
end
