function [re, im,sol,data,pmat,irot] = equiv_hopf_var(prob,oid,omega,n_osc,cycles)
id=coco_get_id(oid,'ep');
[data,u0,t0]=coco_get_func_data(prob,id,'data','u0','t0');
xhopf=u0(data.pr.ep_eqn.x_idx);
phopf=u0(data.pr.ep_eqn.p_idx);
alldim=length(xhopf);
basedim=alldim/n_osc;
id=eye(alldim);
dfdx=data.pr.ode_DFDX;
Jhopf_om=dfdx(data,prob,xhopf,phopf)-1i*omega*id;
[pmat,irot]=cycle2perm([n_osc,basedim],cycles);
nvec=null([Jhopf_om;pmat-diag(exp(1i*2*pi./irot))]);
assert(size(nvec,2)==1,'dimension of nullspace=%d\neq1',size(nvec,2));
[re,im]=deal(real(nvec),imag(nvec));
if nargout<3
    return
end
[v,w]=deal(re,im);
k=omega^2;
al = [-w'*v v'*v];
al = al/norm(al);
nv = al(1)*v + al(2)*w;
nv = nv/norm(nv);
sol.x  = xhopf;
sol.p  = phopf;
sol.u0 = u0;
sol.t0 = t0;

sol.var.v  = [ v/norm(v) -sqrt(k)*w/norm(v) ];
sol.var.w  = [ -sqrt(k)*w/norm(v) -k*v/norm(v) ];
sol.var.u0 = [ sol.var.v(:) ; sol.var.w(:) ];
sol.var.t0 = [];

sol.hb.k  = k;
sol.hb.nv = nv;
sol.hb.u0 = k;
sol.hb.t0 = [];
end