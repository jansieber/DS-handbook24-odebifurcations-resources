function varargout=sym_sphere_constraint(action,varargin)
%% Automatically generated with matlabFunction
%#ok<*DEFNU,*INUSD,*INUSL>

switch action
  case 'nargs'
   varargout{1}=2;
   return
  case 'nout'
   varargout{1}=1;
   return
  case 'argrange'
   varargout{1}=struct('x',1:1,'p',2:4);
   return
  case 'argsize'
   varargout{1}=struct('x',1,'p',3);
   return
  case 'vector'
   varargout{1}=struct('x',1,'p',1);
   return
  case 'extension'
   varargout{1}='rhs';
   return
  case 'maxorder'
   varargout{1}=2;
   return
end
nout=1;
order=varargin{1};
f=str2func(sprintf('sym_sphere_constraint_%s_%d',action,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{2:end});
end



function out1 = sym_sphere_constraint_rhs_0(in1,in2,in3,in4,in5,in6,in7,in8)
%SYM_SPHERE_CONSTRAINT_RHS_0
%    OUT1 = SYM_SPHERE_CONSTRAINT_RHS_0(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    11-Aug-2023 01:11:11

out1 = in1.^2+in2.^2+in3.^2+in4.^2-1.0;
end


function out1 = sym_sphere_constraint_rhs_1(in1,in2,in3,in4,in5,in6,in7,in8)
%SYM_SPHERE_CONSTRAINT_RHS_1
%    OUT1 = SYM_SPHERE_CONSTRAINT_RHS_1(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    11-Aug-2023 01:11:11

out1 = in1.*in5.*2.0+in2.*in6.*2.0+in3.*in7.*2.0+in4.*in8.*2.0;
end


function out1 = sym_sphere_constraint_rhs_2(in1,in2,in3,in4,in5,in6,in7,in8)
%SYM_SPHERE_CONSTRAINT_RHS_2
%    OUT1 = SYM_SPHERE_CONSTRAINT_RHS_2(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    11-Aug-2023 01:11:11

out1 = in5.^2.*2.0+in6.^2.*2.0+in7.^2.*2.0+in8.^2.*2.0;
end

