function varargout=sys_brus(action,varargin)
%% Automatically generated with matlabFunction
%#ok<*DEFNU,*INUSD,*INUSL>

switch action
  case 'nargs'
   varargout{1}=2;
   return
  case 'nout'
   varargout{1}=8;
   return
  case 'argrange'
   varargout{1}=struct('x',1:8,'p',9:11);
   return
  case 'argsize'
   varargout{1}=struct('x',8,'p',3);
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
nout=8;
order=varargin{1};
f=str2func(sprintf('sys_brus_%s_%d',action,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{2:end});




function [out1,out2,out3,out4,out5,out6,out7,out8] = sys_brus_rhs_0(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20,in21,in22)
%SYS_BRUS_RHS_0
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6,OUT7,OUT8] = SYS_BRUS_RHS_0(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20,IN21,IN22)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    08-Aug-2024 16:01:54

t2 = in10+1.0;
t3 = in1.^2;
t4 = in3.^2;
t5 = in5.^2;
t6 = in7.^2;
t7 = in2.*t3;
out1 = in9+t7+(in11.*(in1.*-3.0+in3+in5+in7))./1.0e+3-in1.*t2;
if nargout > 1
    out2 = -t7+(in11.*(in2.*-3.0+in4+in6+in8))./1.0e+2+in1.*in10;
end
if nargout > 2
    t8 = in4.*t4;
    out3 = in9+t8+(in11.*(in1-in3.*3.0+in5+in7))./1.0e+3-in3.*t2;
end
if nargout > 3
    out4 = -t8+(in11.*(in2-in4.*3.0+in6+in8))./1.0e+2+in3.*in10;
end
if nargout > 4
    t9 = in6.*t5;
    out5 = in9+t9+(in11.*(in1+in3-in5.*3.0+in7))./1.0e+3-in5.*t2;
end
if nargout > 5
    out6 = -t9+(in11.*(in2+in4-in6.*3.0+in8))./1.0e+2+in5.*in10;
end
if nargout > 6
    t10 = in8.*t6;
    out7 = in9+t10+(in11.*(in1+in3+in5-in7.*3.0))./1.0e+3-in7.*t2;
end
if nargout > 7
    out8 = -t10+(in11.*(in2+in4+in6-in8.*3.0))./1.0e+2+in7.*in10;
end


function [out1,out2,out3,out4,out5,out6,out7,out8] = sys_brus_rhs_1(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20,in21,in22)
%SYS_BRUS_RHS_1
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6,OUT7,OUT8] = SYS_BRUS_RHS_1(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20,IN21,IN22)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    08-Aug-2024 16:01:54

t2 = in1.*in21;
t3 = in3.*in21;
t4 = in5.*in21;
t5 = in7.*in21;
t6 = in10+1.0;
t7 = in1.^2;
t8 = in3.^2;
t9 = in5.^2;
t10 = in7.^2;
t15 = in1.*in2.*in12.*2.0;
t16 = in3.*in4.*in14.*2.0;
t17 = in5.*in6.*in16.*2.0;
t18 = in7.*in8.*in18.*2.0;
t11 = in13.*t7;
out1 = in20-t2+t11+t15+(in22.*(in1.*-3.0+in3+in5+in7))./1.0e+3+(in11.*(in12.*-3.0+in14+in16+in18))./1.0e+3-in12.*t6;
if nargout > 1
    out2 = t2-t11-t15+(in22.*(in2.*-3.0+in4+in6+in8))./1.0e+2+(in11.*(in13.*-3.0+in15+in17+in19))./1.0e+2+in10.*in12;
end
if nargout > 2
    t12 = in15.*t8;
    out3 = in20-t3+t12+t16+(in22.*(in1-in3.*3.0+in5+in7))./1.0e+3+(in11.*(in12-in14.*3.0+in16+in18))./1.0e+3-in14.*t6;
end
if nargout > 3
    out4 = t3-t12-t16+(in22.*(in2-in4.*3.0+in6+in8))./1.0e+2+(in11.*(in13-in15.*3.0+in17+in19))./1.0e+2+in10.*in14;
end
if nargout > 4
    t13 = in17.*t9;
    out5 = in20-t4+t13+t17+(in22.*(in1+in3-in5.*3.0+in7))./1.0e+3+(in11.*(in12+in14-in16.*3.0+in18))./1.0e+3-in16.*t6;
end
if nargout > 5
    out6 = t4-t13-t17+(in22.*(in2+in4-in6.*3.0+in8))./1.0e+2+(in11.*(in13+in15-in17.*3.0+in19))./1.0e+2+in10.*in16;
end
if nargout > 6
    t14 = in19.*t10;
    out7 = in20-t5+t14+t18+(in22.*(in1+in3+in5-in7.*3.0))./1.0e+3+(in11.*(in12+in14+in16-in18.*3.0))./1.0e+3-in18.*t6;
end
if nargout > 7
    out8 = t5-t14-t18+(in22.*(in2+in4+in6-in8.*3.0))./1.0e+2+(in11.*(in13+in15+in17-in19.*3.0))./1.0e+2+in10.*in18;
end


function [out1,out2,out3,out4,out5,out6,out7,out8] = sys_brus_rhs_2(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20,in21,in22)
%SYS_BRUS_RHS_2
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6,OUT7,OUT8] = SYS_BRUS_RHS_2(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20,IN21,IN22)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    08-Aug-2024 16:01:54

t2 = in12.^2;
t3 = in14.^2;
t4 = in16.^2;
t5 = in18.^2;
t6 = in12.*in21.*2.0;
t7 = in14.*in21.*2.0;
t8 = in16.*in21.*2.0;
t9 = in18.*in21.*2.0;
t10 = in1.*in12.*in13.*4.0;
t11 = in3.*in14.*in15.*4.0;
t12 = in5.*in16.*in17.*4.0;
t13 = in7.*in18.*in19.*4.0;
t14 = in2.*t2.*2.0;
out1 = -t6+t10+t14+(in22.*(in12.*-3.0+in14+in16+in18))./5.0e+2;
if nargout > 1
    out2 = t6-t10-t14+(in22.*(in13.*-3.0+in15+in17+in19))./5.0e+1;
end
if nargout > 2
    t15 = in4.*t3.*2.0;
    out3 = -t7+t11+t15+(in22.*(in12-in14.*3.0+in16+in18))./5.0e+2;
end
if nargout > 3
    out4 = t7-t11-t15+(in22.*(in13-in15.*3.0+in17+in19))./5.0e+1;
end
if nargout > 4
    t16 = in6.*t4.*2.0;
    out5 = -t8+t12+t16+(in22.*(in12+in14-in16.*3.0+in18))./5.0e+2;
end
if nargout > 5
    out6 = t8-t12-t16+(in22.*(in13+in15-in17.*3.0+in19))./5.0e+1;
end
if nargout > 6
    t17 = in8.*t5.*2.0;
    out7 = -t9+t13+t17+(in22.*(in12+in14+in16-in18.*3.0))./5.0e+2;
end
if nargout > 7
    out8 = t9-t13-t17+(in22.*(in13+in15+in17-in19.*3.0))./5.0e+1;
end

