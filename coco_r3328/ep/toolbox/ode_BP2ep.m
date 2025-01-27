function prob = ode_BP2ep(prob, oid, varargin)
%ODE_BP2EP   Switch to secondary branch of equilibrium points at branch point.
%
% PROB = ODE_BP2EP(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-ep-end' | '-end-ep' | '-var' VECS | '-no-var' | '-no-pars'}
%
% Start a continuation of equilibrium points along a secondary branch
% intersecting a previously computed branch with name RUN in a branch
% point. To start from a saved branch point, at least the name RUN of the
% continuation run and the solution label LAB must be given. The label LAB
% must be the label of a branch point.
%
% The arguments and their meaning are identical to ODE_EP2EP.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string).
% RUN  : Run identifier (string or cell-array of strings).
% SOID : Source object instance identifier (string, optional).
% LAB  : Solution label (integer).
%
% OPTS : '-ep-end', '-end-ep',  and '-var' VECS (optional, multiple options
%        may be given). Either '-ep-end' or '-end-ep' mark the end of input
%        to ODE_BP2EP. The option '-var' indicates the inclusion of the
%        variational problem J*v=w where the initial solution guess for v
%        is given by the content of VECS.
% 
% See also: ODE_EP2EP, EP_READ_SOLUTION, EP_ADD, EP_ADD_VAR
%
% Use of this function is deprecated.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_BP2ep.m 3317 2025-01-07 21:10:49Z hdankowicz $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}  
  };
opts_spec = {
   '-ep-end',       '',    '',    'end', {}
   '-end-ep',       '',    '',    'end', {}
      '-var',   'vecs',    [],   'read', {}
   '-no-var',  'novar', false, 'toggle', {}
  '-no-pars', 'nopars', false, 'toggle', {}
  };
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:});

prob = coco_set(prob, 'cont', 'branch', 'switch');
warning('Use of ODE_BP2EP is deprecated. Reset ''branch'' property');

[sol, data] = ep_read_solution(args.soid, args.run, args.lab);

if opts.nopars
  data.pnames = {};
end

data = ode_init_data(prob, data, oid, 'ep');

if ~isempty(opts.vecs)
  assert(isnumeric(opts.vecs) && data.xdim == size(opts.vecs,1), ...
    '%s: incompatible specification of vectors of perturbations', ...
    mfilename);
  [prob, data] = ep_add(prob, data, sol, '-cache-jac');
  sol.var = struct('v', opts.vecs, 'w', data.ep_Jx*opts.vecs);
  sol.var.u0 = [sol.var.v; sol.var.w];
  sol.var.t0 = [];
  prob = ep_add_var(prob, data, sol);
  prob = ode_add_tb_info(prob, oid, 'ep', 'ep', 'ep', ep_sol_info('VAR'));
elseif isfield(sol, 'var') && ~opts.novar
  [prob, data] = ep_add(prob, data, sol, '-cache-jac');
  prob = ep_add_var(prob, data, sol.var.v);
  prob = ode_add_tb_info(prob, oid, 'ep', 'ep', 'ep', ep_sol_info('VAR'));
else
  prob = ep_add(prob, data, sol);
  prob = ode_add_tb_info(prob, oid, 'ep', 'ep', 'ep', ep_sol_info());
end

end
