function dy = d_symmetry(data, ~, ~, ~, ~)
%% derivative of monitoring function for symmetry
% ys=x(0)-gamma*x(tau)
% y=ys(i);
%

m    = data.coll_seg.int.NCOL;
xdim = data.coll_seg.int.dim;
N    = data.coll_seg.maps.NTST;
pdim = data.coll_seg.maps.pdim;
tmi  = data.coll_seg.mesh.tmi;
tbp  = data.coll_seg.mesh.tbp;

dy = zeros(xdim, xdim*(m+1)*N+2+pdim);

%%first point

dy(:,1:xdim) = eye(xdim);

%%second point

% find interval
trs    = data.tau;
J = min(N, find(N*trs>=tmi, 1, 'last'));
tint   = tbp((J-1)*(m+1)+(1:m+1));
ts     = linspace(-1, 1, m+1)';

% interpolated point
s  = repmat(2*(trs-tint(1))/(tint(end)-tint(1))-1, [1 m+1 m+1]);
sj = repmat(reshape(ts, [1 m+1 1]), [1 1 m+1]);
sk = repmat(reshape(ts, [1 1 m+1]), [1 m+1 1]);

t1 = s-sk;
t2 = sj-sk;
idx = find(abs(t2)<=eps);
t1(idx) = 1;
t2(idx) = 1;

dy(:, (J-1)*xdim*(m+1)+(1:xdim*(m+1))) = -data.gamma*kron(prod(t1./t2, 3), eye(xdim));
dy = dy(data.i,:);

end
