%% Use of symbolic derivatives  for coco demo bistable
% (see coco demo bistable). The right-hand side for the ODE (the Duffing
% oscillator) was generated in <gen_sym_bistable.html>
%
% $$
%   \frac{\mathrm{d}}{\mathrm{d} t}\vec{x}=\vec{f}(t,\vec{x},\vec{p})
% $$ 
%
% where 
%
% $$
%  \vec{x}=
%   \left[\begin{array}{c}
%     x\\ v
%   \end{array}\right],\quad 
%   \vec{p}=\left[
%   \begin{array}{c}
%     T\\ a\\ \gamma
%   \end{array}\right],\quad
%   \vec{f}=\left[
%   \begin{array}{l}
%     v\\ -\gamma v-x-x^3+a \cos(2 \pi t/T)
%   \end{array}\right]
% $$
%% Load path
clear
format compact
%% Extract functions with sco_gen
% After the function file |sym_bistable.m| has been generated, the symbolic
% toolbox is no longer needed. The function file is usually not directly
% called by the user. Instead, one uses the wrapper |sco_gen| to define the
% right-hand side and its partial derivatives.  The input format of the
% right-hand side was set by the 2nd and 3rd arguments of the call to
% |sco_sym2funcs|.
%
% The first line calls |sco_gen| to extract the right-hand side $f$ from
% |sym_bistable|. The empty-string second argument indicates that no
% partial derivatives are requested such that the returned function handle
% |f| is will return the values of $f$. Partial derivatives of $f$ with
% respect to its arguments are generated by giving the name of the argument
% of $f$ with respect to which differentiation should occur as the second
% input. Line 2 demonstrates this, defining $\partial
% f_x(t,x,p)$.
f=sco_gen(@sym_bistable,'');    % define f: y=f(t,x,p)
dfx=sco_gen(@sym_bistable,'x'); % define df/dx: J=df(t,x,p)
%%  Initial run with IVP solver
% This section runs an initial-value problem ODE solver for the Duffing
% oscillator problem using |ode45|.  Line 2 and 3 call the standard Matlab
% function |ode45| to solve the ODE using the constructed |f| as its
% argument (with parameters |p0| as fixed in line 1). 
p0=[2*pi; 0.015; 0.04];         % set initial parameter values
[~,x0]=ode45(@(t,x)f(t,x,p0), [0 20*pi], [0;1]); % Transients, starting at x=[0;1]
[t0,x0]=ode45(@(t,x)f(t,x,p0), [0 2*pi], x0(end,:)'); % Approximate periodic orbit
%% Demonstrate vectorized calls to function and derivatives
% The function |f| may be called with arguments |t| of shape |1xN|, |x| of
% shape |2xN|, |p| of shape |3xN|, returning, for example, in line 1, |f_t|
% of shape |2xN|, where |N| is an arbitrary number, permitting vectorized
% calls to |f|. Line 1 demonstrates this call by evaluating |f| along the
% solution obtained by the IVP solver (note the transposing of the inputs
% |t0| and |x0|).
%
% Inputs format of |dfx| is the same as for |f|, but the output has format
% |2x2xN|, returning the |2x2| partial derivative $\partial_xf$ in the
% arguments. Line 2 below demonstrates a vectorized call to |dfx|
% along the solution, where output |df_t| has shape
% |[2,2,size(x0',2)]|.
%
% The plot shows the position, velocity, force and instantaneous spring
% stiffness, thus extracted.
f_t=f(t0',x0',p0);              % values of f(t,x(t),p) along orbit
dfx_t=dfx(t0',x0',p0);           % values of df/dx(t,x(t),p) along orbit
fprintf('size(f_t)=[%d,%d], ',size(f_t));
fprintf('size(dfx_t)=[%d,%d,%d]\n',size(dfx_t));
clf;plot(t0,x0(:,1),'DisplayName','position');hold on;legend;
plot(t0,x0(:,2),'DisplayName','velocity');
plot(t0,f_t(2,:),'DisplayName','force');
plot(t0,squeeze(dfx_t(2,1,:)),'DisplayName','spring stiffness')
hold off
xlabel('time');
%% Second-order derivatives
% Second order partial derivatives are generated by giving a cell array
% with two names as the second argument of |sco_gen|. The input format of
% the resulting function is the same as for |f|, but, for example for
% |dfxp|$=\partial^2_{x,p}f(t,x,p)$, the output has shape |2x2x3xN|. In
% this example, only one of the mixed second derivatives is non-zero (equal
% |-1|): $\partial_{v,\gamma}f_2$, corresponding to |dfxp_t(2,2,3,:)|. 
dfxp=sco_gen(@sym_bistable,{'x','p'});
dfxp_t=dfxp(t0',x0',p0);
fprintf('size(dfxp_t)=[%d,%d,%d,%d]\n',size(dfxp_t));
df2_vgam=squeeze(dfxp_t(2,2,3,:)).' %#ok<*NOPTS>
%% Scalar arguments
% Outputs of derivatives with respect to inputs that were declared as
% scalar in the call to |sco_sym2funcs| have an output format of reduced
% shape. In this example, the 1st argument of |f|, |'t'|, is of this type.
% This exception accommodates conventions of the COLL toolbox.
dft=sco_gen(@sym_bistable,'t');
dftt=sco_gen(@sym_bistable,{'t','t'});
dftp=sco_gen(@sym_bistable,{'t','p'});
dfxt=sco_gen(@sym_bistable,{'x','t'});
fprintf('size(dft_t)=[%d,%d]\n',size(dft(t0',x0',p0)));
fprintf('size(dftt_t)=[%d,%d]\n',size(dftt(t0',x0',p0)));
fprintf('size(dfxt_t)=[%d,%d,%d]\n',size(dfxt(t0',x0',p0)));
%% Shortcut
% The call to |sco_gen| with only its first input is a shortcut for the
% creating the function in line 2 below. This shortcut makes the notation
% more succinct when a sequence of partial derivatives is passed on as an
% argument to coco constructors. For example, the cell array |funcs|
% collects $f$ and its first partial derivatives wit hrespet to |'x'| and
% |'p'|. This cell array is the argument required for the coco demo
% |bistable|.
F=sco_gen(@sym_bistable);                        % F and F2 are the same
F2=@(varargin)sco_gen(@sym_bistable,varargin{:});% F and F2 are the same
funcs = {F(''),F('x'),F('p')}; %#ok<NASGU>
%% Test mixed and high derivatives
% These are not currently used in the demo, but are useful for normal form
% coefficient monitoring. The examples below creates
% 
% $$\mathtt{dfxpdir}:(t,x,p,\delta_t,\delta_x,\delta_p)\to\partial_{xpt}f(t,x,p)\delta_t+
%    \partial_{xpx}f(t,x,p)\delta_x+\partial_{xpp}f(t,x,p)\delta_p
%  $$
%
% with result in $n_f\times n_x\times n_p\times n_\mathrm{vec}$ and.
% 
% $$\mathtt{dfxxvp}:(t,x,p,\delta_x)\to\partial_{xxp}f(t,x,p)\delta_x
% $$
%
% also with result in $n_f\times n_x\times n_p\times n_\mathrm{vec}$.
% 
df3=F(3); % third-order directional derivative
dfxpdir=@(t,x,p,dt,dx,dp)feval(F({3}),t,x,p,{0,0,dt},{'I',0,dx},{0,'I',dp})
argbase={t0',x0',p0};
dev={1,0*x0'+1,0*p0+2};
dfxpdirval=dfxpdir(argbase{:},dev{:});
size(dfxpdirval)
dfxxvp=F({'x','x*v','p'});
xdev=0*x0'+1;
dfxxpval=dfxxvp(argbase{:},xdev)
size(dfxxpval)
% compare to total derivative
J3=df3(argbase{:});
dev_exp=sco_argexpand(dev);
df3uval=sum(J3.*reshape(cat(1,dev_exp{:}),1,1,1,[],length(t0)),4);
dfxpdirval2=squeeze(df3uval(:,2:3,4:6,1,:));
max(abs(dfxpdirval-dfxpdirval2),[],'all')
xdev_exp=sco_argexpand({0,xdev,0*p0});
dfxxpuval=sum(J3.*reshape(cat(1,xdev_exp{:}),1,1,1,[],length(t0)),4);
dfxxpval2=squeeze(dfxxpuval(:,2:3,4:6,1,:));
max(abs(dfxxpval-dfxxpval2),[],'all')
%% Run of coco demo
% Most of the above parts were only an illustration of using symboilc
% derivatives. The coco demo bistable only needs the outputs of the IVP
% solution and  the definition of the array |funcs|.
% Below is the remainder of the original coco demo bistable. At the
% beginning we repeat the commands needed from above.
%% Define r.h.s. and derivatives, and use IVP solution as initial guess
clear
F=sco_gen(@sym_bistable);       %shortcut for function generator
funcs = {F(''),F('x'),F('p')};  % r.h.s. and partial derivatives
p0=[2*pi; 0.015; 0.04];         % set initial parameter values
[~,x0]=ode45(@(t,x)funcs{1}(t,x,p0),[0 20*pi],[0;1]); % Transients, starting at x=[0;1]
[t0,x0]=ode45(@(t,x)funcs{1}(t,x,p0),[0 2*pi],x0(end,:)'); % Approximate periodic orbit
%% Define coco problem for initial response curve
prob = coco_prob();
prob = coco_set(prob, 'ode', 'autonomous', false);
coll_args = [funcs, {t0, x0, {'T' 'A' 'd'}, p0}];
prob = ode_isol2po(prob, '', coll_args{:});
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_glue(prob, 'glue', uidx(maps.T_idx), uidx(maps.p_idx(1)));
prob = po_mult_add(prob, 'po.orb');  % Store Floquet multipliers with bifurcation data
cont_args = {1, {'po.period' 'T'}, [2*pi/1.3 2*pi/0.7]};

fprintf('\n Run=''%s'': Continue primary family of periodic orbits.\n', ...
  'freq_resp');

bd1  = coco(prob, 'freq_resp', [], cont_args{:});

%% Continue curve of saddle-node bifurcations
%
% After imposition of gluing conditions, the continuation problem encoded
% below has a dimensional deficit of -2. A one-dimensional family of
% saddle-node bifurcations of periodic orbits is obtained by releasing
% 'po.period', 'T', and 'A', and allowing these to vary during
% continuation.

labs = coco_bd_labs(bd1, 'SN');

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 25);
prob = ode_SN2SN(prob, '', 'freq_resp', labs(1));
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_glue(prob, 'glue', uidx(maps.T_idx), uidx(maps.p_idx(1)));
prob = coco_set(prob, 'cont', 'NAdapt', 5);
prob = coco_add_event(prob, 'UZ', 'A', 0.01:0.005:0.1);
cont_args = {1, {'po.period' 'T' 'A'}, [2*pi/1.3 2*pi/0.7]};

fprintf(...
  '\n Run=''%s'': Continue saddle-node bifurcations from point %d in run ''%s''.\n', ...
  'saddle-node', labs(1), 'freq_resp');

bd2  = coco(prob, 'saddle-node', [], cont_args{:});

%% Sweep frequency-response curves
%
% For each of a sample of points on the curve of saddle-node bifurcations,
% the continuation problem encoded below has dimensional deficit -1.  A
% one-dimensional family of periodic orbits is obtained by releasing
% 'po.period' and 'T', and allowing these to vary during continuation.

labs  = coco_bd_labs(bd2, 'UZ');

for lab=labs
  prob = coco_prob();
  prob = coco_set(prob, 'ode', 'autonomous', false);
  prob = coco_set(prob, 'po', 'bifus', 'off');
  prob = ode_po2po(prob, '', 'saddle-node', lab, '-no-var');
  [data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
  maps = data.coll_seg.maps;
  prob = coco_add_glue(prob, 'glue', uidx(maps.T_idx), uidx(maps.p_idx(1)));
  prob = coco_set(prob, 'cont', 'NAdapt', 10);
  cont_args = {1, {'po.period' 'T'}, [2*pi/1.3 2*pi/0.7]};
  
  fprintf(...
  '\n Run=''%s'': Continue family of periodic orbits from point %d in run ''%s''.\n', ...
  sprintf('lab=%d', lab), lab, 'saddle-node');

  coco(prob, sprintf('lab=%d', lab), [], cont_args{:});
end


%% Free vibrations - backbone
%
% In this continuation problem, we seek periodic responses in the absence
% of excitation. These only exist for zero damping, and constitute a
% one-dimensional family of orbits parameterized by the orbit period. The
% encoding includes four monitor functions that evaluate to the problem
% parameters and the first component of the initial point on the periodic
% orbit, and the corresponding inactive continuation parameters 'T', 'A',
% 'd', and 'y0'. Holding 'y0' fixed corresponds to the imposition of a
% Poincare condition, which results in an overall dimensional deficit of
% -1. A one-dimensional family of periodic orbits is obtained by releasing
% 'po.period' and 'd', and allowing these to vary during continuation.

t0 = (0:0.01:2*pi)';
x0 = 2e-2*[sin(t0) cos(t0)];
p0 = [2*pi; 0; 0];
prob = coco_prob();
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'po', 'bifus', 'off');
funcs = {F(''),F('x'),F('p')};
coll_args = [funcs, {t0, x0, {'T' 'A' 'd'}, p0}];
prob = ode_isol2po(prob, '', coll_args{:});
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_pars(prob, 'section', uidx(maps.x0_idx(1)), 'y0');
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'PtMX', [0 100]);
cont_args = {1, {'po.period'  'd'}, [2*pi/1.3 2*pi/0.7]};

fprintf('\n Run=''%s'': Continue backbone curve of periodic orbits.\n', ...
  'backbone');

coco(prob, 'backbone', [], cont_args{:});

%% Graphical representation of stored solutions

figure(1); clf; hold on; grid on; box on
coco_plot_bd('saddle-node', 'T', @(T) 2*pi./T, ...
    '||po.orb.x||_{L_2[0,T]}')
thm = struct('ustab', '', 'lspec', {{'r', 'LineWidth', 2}});
for lab=labs
  coco_plot_bd(thm, sprintf('lab=%d', lab), 'T', @(T) 2*pi./T, ...
    '||po.orb.x||_{L_2[0,T]}')
end
thm = struct('ustab', '', 'lspec', {{'b--', 'LineWidth', 2}});
coco_plot_bd(thm, 'backbone', 'po.period', @(T) 2*pi./T, ...
    '||po.orb.x||_{L_2[0,T]}')
axis([0.7 1.3 0 inf]); hold off
