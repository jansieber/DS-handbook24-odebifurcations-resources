function  y = saddle_q(data, xbp, T0, T, par,funcs) %#ok<INUSD>
m    = data.coll_seg.int.NCOL;
xdim = data.coll_seg.int.dim;
xbp=reshape(xbp,xdim,[]);
f = feval(funcs.f, xbp, repmat(par, [1 size(xbp, 2)])); % Evaluate vector field at basepoints
[fmin, idx] = min(sqrt(sum(f.*f, 1))); % Find basepoint closest to equilibrium
xeq=xbp(:,idx);
df=feval(funcs.dfdx, xeq, par); % Evaluate Jacobian at equilibrium;
y=[xeq;fmin;det(df);trace(df)];
end