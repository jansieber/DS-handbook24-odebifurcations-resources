function tb=prob_fcn_info(p,varargin)
%% extract information about functions already defined in coco problem
ids=p.efunc.identifyers;
if isempty(varargin)
    vnames={'type','fidx','uidx','midx','data'};
else
    vnames=varargin;
end
for i=length(vnames):-1:1
    tbc(i,:)=cellfun(@(s)coco_get_func_data(p,s,vnames{i}),ids,'uniformoutput',false);
end
tb=cell2table(tbc', 'RowNames',ids,'VariableNames',vnames);
end
