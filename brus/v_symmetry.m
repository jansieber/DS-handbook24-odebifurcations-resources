function y=v_symmetry(data,x,p,v,k)
om = sqrt(k);
if om<1e-6 % if neutral saddle
  y = NaN;
  return
end
A  = data.dfdxhan(x,p);
va = v-1i*A*v/om;
id=eye(length(v));
S=data.gamma-id*exp(1i*2*pi/data.tau);
yc=S*va;
y=[real(yc);imag(yc)];
end
