function [prob, status, xtr] = ampremesh(prob, data, chart, old_u, old_V) %#ok<INUSD>

[fdata, uidx] = coco_get_func_data(prob, data.cid, 'data', 'uidx');
seg  = fdata.coll_seg;
maps = seg.maps;
data.coll_seg = seg;

xtr       = [];
prob      = coco_change_func(prob, data, 'uidx', uidx(maps.xbp_idx));
status    = 'success';

end
