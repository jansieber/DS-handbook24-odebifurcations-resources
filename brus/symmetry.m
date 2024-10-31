function y = symmetry(data, xbp, ~, ~, ~)
%% add monitoring function
% ys=x(0)-gamma*x(tau)
% y=ys(i);
%
m    = data.coll_seg.int.NCOL;
xdim = data.coll_seg.int.dim;
N    = data.coll_seg.maps.NTST;
tmi  = data.coll_seg.mesh.tmi;
tbp  = data.coll_seg.mesh.tbp;

%%first point

x1 = xbp(1:xdim);

%%second point

% find interval
trs    = data.tau; %data.tau = 1/2
J = min(N, find(N*trs>=tmi, 1, 'last'));
xbpint = xbp((J-1)*xdim*(m+1)+(1:xdim*(m+1)));
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

x2 = kron(prod(t1./t2, 3), eye(xdim))*xbpint;

y = x1-data.gamma*x2; % data.gamma = [-1 0 0;0 -1 0;0 0 0]
y = y(data.i(:));

end
