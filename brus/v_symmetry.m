function y=v_symmetry(data,x,p,v,k)
om = sqrt(k);
if om<1e-9 % if neutral saddle
  y = NaN;
  return
end
A  = data.dfdxhan(x,p);
vc = om*v-1i*A*v;
S=data.gamma-diag(exp(1i*2*pi./data.tau));
yc=S*vc;
y=[real(yc);imag(yc)];
end
