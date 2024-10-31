%% test branching off from Hopf bifurcation
clear
%% define r.h.s.
if exist('sys_brus', 'file')~=2
  syms y1 y2 y3 y4 y5 y6 y7 y8 A B lambda
  F = sco_sym2funcs([
    A-(B+1)*y1+y1^2*y2+lambda/1000*(-4*y1+y1+y3+y5+y7); ...
    B*y1-y1^2*y2+lambda/100*(-4*y2+y2+y4+y6+y8); ...
    A-(B+1)*y3+y3^2*y4+lambda/1000*(-4*y3+y1+y3+y5+y7); ...
    B*y3-y3^2*y4+lambda/100*(-4*y4+y2+y4+y6+y8); ...
    A-(B+1)*y5+y5^2*y6+lambda/1000*(-4*y5+y1+y3+y5+y7); ...
    B*y5-y5^2*y6+lambda/100*(-4*y6+y2+y4+y6+y8); ...
    A-(B+1)*y7+y7^2*y8+lambda/1000*(-4*y7+y1+y3+y5+y7); ...
    B*y7-y7^2*y8+lambda/100*(-4*y8+y2+y4+y6+y8)], ...
    {[y1; y2; y3; y4; y5; y6; y7; y8], [A; B; lambda]}, ...
    {'x', 'p'}, 'filename', 'sys_brus');
  funcs = {F(''), F('x'), F('p')};
else
  F = sco_gen(@sys_brus);
  funcs = {F(''), F('x'), F('p')};
end
%% Set common parameters
basedim=2;
n_osc=4;
alldim=n_osc*basedim;
A = 2;
B = 5.9;
%% create initial periodic solution guess
delta = B-1-A^2;
omega = sqrt(11*A^2*(11-9*delta)-100*delta^2)/11;

period = 2*pi/omega;
epsilon = 1e-3;

mu = A^2;
nu = 1i*omega-A^2-10*delta/11;

vec = [mu; nu; (1i)^3*mu; (1i)^3*nu; 1i*mu; 1i*nu; (1i)^2*mu; (1i)^2*nu];
re = real(vec);
im = imag(vec);

t0 = 0:period/400:period;
xhopf=[A; B/A; A; B/A; A; B/A; A; B/A];
[pnames,phopf]=deal(...
    {'A', 'B',   'lambda'},...
     [A;   B; 1000*delta/44]);
%%
%cycles={[1,2,4,3],4};
%cycles={[2,1],2; [3,4],1};
cycles={[2,3,4],1};
prob=hopf_equiv_hb('',funcs,xhopf,pnames,phopf,omega,n_osc,cycles);
[data,uidx]=coco_get_func_data(prob,'ep.HB','data','uidx');
prob = coco_add_pars(prob, 'kmon', uidx(data.pr.ep_hb.k_idx), 'k');
prob = coco_set(prob, 'cont', 'PtMX', [20 30], 'NPR', 100,'norm', inf);
coco(prob, 'runHB', [], 1, {'lambda','B','k'}, [15 30]);

%%
