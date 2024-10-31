function [pmat,irot,cycles]=cycle2perm(dims,cycles)
%% create permutation matrix and rotation numbers from cycle representation
% dims is vector of dimensions:[number_of_nodes,dim_of_node]
% cycles is ncyc x 2 cell array, where cycles{i,1} is cycle for
% permutation, cycles{i,2} is time shift
alldim=prod(dims);
n_osc=dims(1);
pmat=zeros(n_osc);
if ~iscell(cycles)
    cycles={cycles,length(cycles)};
end
unchanged=setdiff(1:n_osc,cat(2,cycles{:,1}));
pmat(sub2ind([n_osc,n_osc],unchanged,unchanged))=1;
for i=1:size(cycles,1)
    c=cycles{i,1};
    pmat(sub2ind([n_osc,n_osc],[c(2:end),c(1)],c))=1;
end
pmat=kron(pmat,eye(dims(2)));
if nargout==1
    return
end
irot=ones(1,n_osc);
for i=1:size(cycles,1)
    irot(cycles{i,1})=cycles{i,2};
end
irot=reshape(repmat(irot,dims(2),1),1,alldim);
end