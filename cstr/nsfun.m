function [y,J]=nsfun(f,dfx,u,iv,uref,sel,output)
y=locfun(f,dfx,u,iv,uref,sel);
if nargout<2 && nargin<7
    return
end
J=ScJacobian(@(ua)locfun(f,dfx,ua,iv,uref,sel),u);
if nargin>=7 && strcmp(output,'J')
    [y,J]=deal(J,y);
end
end
function y=locfun(f,dfx,u,iv,uref,sel)
x=u(iv.x);
p=u(iv.p);
k=u(iv.k);
v=u(iv.v);
y1=f(x,p);
J=dfx(x,p);
y2=J^2*v+k*v;
y3=v'*v-1;
vref=uref(iv.v);
y4=vref'*[0,-1;1,0]*v;
y5=u(sel)-uref(sel);
y=[y1;y2;y3;y4;y5];
end