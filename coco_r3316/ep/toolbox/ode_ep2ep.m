function prob = ode_ep2ep(prob, oid, varargin)
%ODE_EP2EP   Reconstruct equilibrium point continuation problem from saved solution.
%
% PROB = ODE_EP2EP(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-switch' | '-ep-end' | '-end-ep' | '-var' VECS | '-no-var' | '-no-pars'}
%
% Reconstruct continuation problem for dynamic or static equilibria from an
% equilibrium point associated with the same vector field that was computed
% and saved to disk in a previous continuation run. At least the name RUN
% of the previous continuation run and the solution label LAB denoting the
% previously computed equilibrium point must be given.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string). Pass the empty string
%        '' for a simple continuation of equilibria.
%
% See ODE_ISOL2EP for more details on PROB and OID.
%
% RUN  : Run identifier (string or cell-array of strings). Name of the run
%        from which to reconstruct the continuation problem.
% SOID : Source object instance identifier (string, optional). If the
%        argument SOID is omitted, OID is used. Pass the empty string ''
%        for OID and omit SOID for a simple continuation of equilibria.
%        Pass non-trivial object identifiers if an instance of the EP
%        toolbox is part of a composite continuation problem.
% LAB  : Solution label (integer). The integer label assigned by COCO to an
%        equilibrium point during the continuation run RUN.
%
% OPTS : '-switch', '-ep-end', '-end-ep', '-var' VECS, and '-no-var'
%        (optional, multiple options may be given). '-switch' instructs
%        ODE_EP2EP to switch branches at the solution point, which must
%        then be a branch point. Either '-ep-end' or '-end-ep' mark the end
%        of input to ODE_EP2EP. The option '-var' indicates the inclusion
%        of the variational problem J*v=w, where the initial solution guess
%        for v is given by the content of VECS. Alternatively, in the
%        absence of the option '-var', the option '-no-var' indicates the
%        exclusion of the variational problem even if this were included in
%        the continuation run RUN.
%
% See also: ODE_ISOL2EP, EP_READ_SOLUTION, EP_ADD, EP_ADD_VAR

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_ep2ep.m 3314 2024-09-02 15:31:13Z hdankowicz $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {
  '-ep-end',       '',    '',    'end', {}
  '-end-ep',       '',    '',    'end', {}
  '-switch', 'switch', false, 'toggle', {}
     '-var',   'vecs',    [],   'read', {}
  '-no-var',  'novar', false, 'toggle', {}
 '-no-pars', 'nopars', false, 'toggle', {}
  };
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:});

if opts.switch
  prob = ode_BP2ep(prob, oid, args.run, args.soid, args.lab, ...
    '-var', opts.vecs);
  return
end

[sol, data] = ep_read_solution(args.soid, args.run, args.lab);
sol.t0 = [];

if opts.nopars
  data.pnames = {};
end

if ~isempty(sol.t0)
  % We restart at branch point but stay on original solution manifold =>
  % reset start direction to tangent vector of primary branch.
  sol.t0 = sol.t;
  if ~coco_exist('NullItMX', 'class_prop', prob, 'cont')
    % Compute improved tangent to new branch on same solution manifold to
    % allow change of parameters on restart at a branch-point.
    prob = coco_set(prob, 'cont', 'NullItMX', 1);
  end
end

data = ode_init_data(prob, data, oid, 'ep');

if ~isempty(opts.vecs)
  assert(isnumeric(opts.vecs) && data.xdim == size(opts.vecs,1), ...
    '%s: incompatible specification of vectors of perturbations', ...
    mfilename);
  [prob, data] = ep_add(prob, data, sol, '-cache-jac');
  prob = ep_add_var(prob, data, opts.vecs);
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
