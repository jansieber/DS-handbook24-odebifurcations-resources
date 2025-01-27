function prob = ode_HB2po(prob, oid, varargin)
%ODE_HB2PO   Switch to branch of periodic orbits at Hopf bifurcation.
%
% PROB = ODE_HB2PO(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-po-end' | '-end-po' | '-var' VECS }
%
% Start a continuation of periodic orbits emanating from a Hopf bifurcation
% that was obtained and saved to dis in a previous continuation along a
% branch of equilibrium points with name RUN. To start from a saved Hopf
% bifurcation, at least the name RUN of the continuation run and the
% solution label LAB must be given. The label LAB must be the label of a
% Hopf bifurcation.
%
% The arguments and their meaning are identical to ODE_PO2PO.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string).
% RUN  : Run identifier (string or cell-array of strings).
% SOID : Source object instance identifier (string, optional).
% LAB  : Solution label (integer).
%
% OPTS : '-po-end', '-end-po', and '-var' VECS (optional, multiple options
%        may be given). Either '-po-end' or '-end-po' mark the end of input
%        to ODE_HB2PO. The option '-var' indicates the inclusion of the
%        variational problem for the corresponding trajectory segment,
%        where the initial solution guess for the perturbations to the
%        initial condition of the orbit is given by the content of VECS.
% 
% See also: ODE_PO2PO, PO_READ_SOLUTION, PO_ADD

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_HB2po.m 3317 2025-01-07 21:10:49Z hdankowicz $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {
  '-po-end',     '', '',  'end', {}
  '-end-po',     '', '',  'end', {}
     '-var', 'vecs', [], 'read', {}
  };
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:});

tbid = coco_get_id(oid, 'po');
po  = po_get_settings(prob, tbid, []);

[sol, data] = ep_read_solution(args.soid, args.run, args.lab);
om = sqrt(sol.hb.k);
vR = sol.var.v(:,1);
vI = sol.var.w(:,1)/(-om);
t0 = linspace(0,2*pi/om, 100)';
x0 = repmat(sol.x', size(t0)) + po.h0*(cos(om*t0)*vR'-sin(om*t0)*vI');
p0 = sol.p;

func = struct('f', data.fhan);
fnames = fieldnames(data);
for i=1:length(fnames)
  if strncmp(fnames(i), 'dfd', 3)
    func.(fnames{i}(1:end-3)) = data.(fnames{i});
  end
end

prob = coco_set(prob, 'cont', 'Valpha', 0);
prob = ode_isol2po(prob, oid, func, t0, x0, data.pnames, p0, ...
  '-var', opts.vecs);

end
