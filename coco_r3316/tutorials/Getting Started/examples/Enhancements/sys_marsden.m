function varargout=sys_marsden(action,varargin)
%% Automatically generated with matlabFunction
%#ok<*DEFNU,*INUSD,*INUSL>

switch action
  case 'nargs'
   varargout{1}=2;
   return
  case 'nout'
   varargout{1}=3;
   return
  case 'argrange'
   varargout{1}=struct('x',1:3,'p',4:5);
   return
  case 'argsize'
   varargout{1}=struct('x',3,'p',2);
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
nout=3;
order=varargin{1};
f=str2func(sprintf('sys_marsden_%s_%d',action,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{2:end});




function [out1,out2,out3] = sys_marsden_rhs_0(x1,x2,x3,p1,p2,x1_dev,x2_dev,x3_dev,p1_dev,p2_dev)
%SYS_MARSDEN_RHS_0
%    [OUT1,OUT2,OUT3] = SYS_MARSDEN_RHS_0(X1,X2,X3,P1,P2,X1_DEV,X2_DEV,X3_DEV,P1_DEV,P2_DEV)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    15-Jun-2023 02:51:10

t2 = x1.^2;
out1 = x2+p2.*t2+p1.*x1;
if nargout > 1
    t3 = -x1;
    out2 = t3+p1.*x2+x2.*x3;
end
if nargout > 2
    out3 = t2+t3-x3+x2.*(p1.^2-1.0);
end


function [out1,out2,out3] = sys_marsden_rhs_1(x1,x2,x3,p1,p2,x1_dev,x2_dev,x3_dev,p1_dev,p2_dev)
%SYS_MARSDEN_RHS_1
%    [OUT1,OUT2,OUT3] = SYS_MARSDEN_RHS_1(X1,X2,X3,P1,P2,X1_DEV,X2_DEV,X3_DEV,P1_DEV,P2_DEV)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    15-Jun-2023 02:51:10

out1 = x2_dev+p1.*x1_dev+p1_dev.*x1+p2_dev.*x1.^2+p2.*x1.*x1_dev.*2.0;
if nargout > 1
    t2 = -x1_dev;
    out2 = t2+p1.*x2_dev+p1_dev.*x2+x2.*x3_dev+x3.*x2_dev;
end
if nargout > 2
    out3 = t2-x3_dev+x2_dev.*(p1.^2-1.0)+x1.*x1_dev.*2.0+p1.*p1_dev.*x2.*2.0;
end


function [out1,out2,out3] = sys_marsden_rhs_2(x1,x2,x3,p1,p2,x1_dev,x2_dev,x3_dev,p1_dev,p2_dev)
%SYS_MARSDEN_RHS_2
%    [OUT1,OUT2,OUT3] = SYS_MARSDEN_RHS_2(X1,X2,X3,P1,P2,X1_DEV,X2_DEV,X3_DEV,P1_DEV,P2_DEV)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    15-Jun-2023 02:51:10

t2 = x1_dev.^2;
out1 = p2.*t2.*2.0+p1_dev.*x1_dev.*2.0+p2_dev.*x1.*x1_dev.*4.0;
if nargout > 1
    out2 = p1_dev.*x2_dev.*2.0+x2_dev.*x3_dev.*2.0;
end
if nargout > 2
    out3 = t2.*2.0+p1_dev.^2.*x2.*2.0+p1.*p1_dev.*x2_dev.*4.0;
end
