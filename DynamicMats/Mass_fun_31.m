function out1 = Mass_fun_31(in1,in2)
%MASS_FUN_31
%    OUT1 = MASS_FUN_31(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    20-Jul-2020 14:25:58

qd1_2 = in2(3);
qd1_4 = in2(7);
qd1_5 = in2(9);
qd1_6 = in2(11);
qd2_2 = in2(4);
qd2_3 = in2(6);
qd2_4 = in2(8);
qd2_5 = in2(10);
qd2_6 = in2(12);
qt1 = in1(1,:);
qt2 = in1(2,:);
t2 = cos(qd1_4);
t3 = cos(qd1_5);
t4 = sin(qd1_4);
t5 = sin(qt1);
t6 = sin(qd1_5);
t7 = cos(qd1_6);
t8 = cos(qt1);
t9 = sin(qd1_6);
t10 = sin(qd2_4);
t11 = sin(qd2_5);
t12 = cos(qd2_4);
t13 = cos(qd2_5);
t14 = sin(qd2_6);
t15 = cos(qt2);
t16 = t4.*t5.*t6;
t17 = t2.*t7.*t8;
t37 = t2.*t3.*t5.*t9;
t18 = t16+t17-t37;
t19 = sin(qt2);
t20 = t8.*t9;
t21 = t3.*t5.*t7;
t22 = t20+t21;
t23 = t10.*t11;
t79 = t12.*t13.*t14;
t24 = t23-t79;
t25 = t10.*t13;
t26 = t11.*t12.*t14;
t27 = t25+t26;
t28 = cos(qd2_6);
t29 = t2.*t5.*t7;
t30 = t2.*t3.*t8.*t9;
t52 = t4.*t6.*t8;
t31 = t29+t30-t52;
t32 = t5.*t9;
t54 = t3.*t7.*t8;
t33 = t32-t54;
t34 = t2.*t5.*t6;
t35 = t3.*t4.*t5.*t9;
t36 = t34+t35-t4.*t7.*t8;
t38 = t18.*t19;
t39 = t15.*t22;
t40 = t38+t39;
t41 = t15.*t18;
t42 = t41-t19.*t22;
t43 = t8.*(3.0./2.0);
t44 = qd2_3.*t2.*t6.*t8;
t45 = qd2_3.*t4.*t5.*t7;
t46 = qd2_2.*t5.*t9.*t19;
t47 = t3.*t7.*t8.*t15.*(3.0./2.0);
t48 = t4.*t6.*t8.*t19.*(3.0./2.0);
t49 = qd2_3.*t3.*t4.*t8.*t9;
t50 = qd2_2.*t4.*t6.*t8.*t15;
t51 = t43+t44+t45+t46+t47+t48+t49+t50-qd1_2.*t5-t5.*t9.*t15.*(3.0./2.0)-t2.*t5.*t7.*t19.*(3.0./2.0)-qd2_2.*t2.*t5.*t7.*t15-qd2_2.*t3.*t7.*t8.*t19-t2.*t3.*t8.*t9.*t19.*(3.0./2.0)-qd2_2.*t2.*t3.*t8.*t9.*t15;
t53 = t19.*t31;
t55 = t15.*t33;
t56 = t53+t55;
t57 = t11.*t12;
t58 = t10.*t13.*t14;
t59 = t57+t58;
t60 = t12.*t13;
t80 = t10.*t11.*t14;
t61 = t60-t80;
t62 = t2.*t6.*t8;
t63 = t4.*t5.*t7;
t64 = t3.*t4.*t8.*t9;
t65 = t62+t63+t64;
t66 = t15.*t31;
t67 = t66-t19.*t33;
t68 = t5.*(3.0./2.0);
t69 = qd1_2.*t8;
t70 = t8.*t9.*t15.*(3.0./2.0);
t71 = qd2_3.*t2.*t5.*t6;
t72 = t2.*t7.*t8.*t19.*(3.0./2.0);
t73 = t3.*t5.*t7.*t15.*(3.0./2.0);
t74 = t4.*t5.*t6.*t19.*(3.0./2.0);
t75 = qd2_2.*t2.*t7.*t8.*t15;
t76 = qd2_3.*t3.*t4.*t5.*t9;
t77 = qd2_2.*t4.*t5.*t6.*t15;
t78 = t68+t69+t70+t71+t72+t73+t74+t75+t76+t77-qd2_3.*t4.*t7.*t8-qd2_2.*t8.*t9.*t19-qd2_2.*t3.*t5.*t7.*t19-t2.*t3.*t5.*t9.*t19.*(3.0./2.0)-qd2_2.*t2.*t3.*t5.*t9.*t15;
t81 = t3.*t4.*t19;
t82 = t2.*t6.*t9.*t19;
t90 = t6.*t7.*t15;
t83 = t81+t82-t90;
t84 = t2.*t3;
t91 = t4.*t6.*t9;
t85 = t84-t91;
t86 = t3.*t4.*t15;
t87 = t6.*t7.*t19;
t88 = t2.*t6.*t9.*t15;
t89 = t86+t87+t88;
out1 = reshape([(t14.*(3.0./2.0)-qd2_2.*t13.*t28).*(t2.*t3.*t14.*3.0-t9.*t14.*t19.*3.0-t4.*t6.*t9.*t14.*3.0+t2.*t7.*t14.*t15.*3.0+t4.*t7.*t11.*t28.*3.0+t9.*t13.*t15.*t28.*3.0-qd1_2.*t4.*t6.*t14.*t15.*2.0+qd1_2.*t3.*t7.*t14.*t19.*2.0+qd1_2.*t2.*t6.*t11.*t28.*2.0-qd2_3.*t3.*t4.*t14.*t19.*2.0+qd2_3.*t6.*t7.*t14.*t15.*2.0-qd2_2.*t2.*t3.*t13.*t28.*2.0+t3.*t4.*t11.*t15.*t28.*3.0+t2.*t7.*t13.*t19.*t28.*3.0+t6.*t7.*t11.*t19.*t28.*3.0+qd1_2.*t2.*t3.*t9.*t14.*t15.*2.0+qd1_2.*t3.*t4.*t9.*t11.*t28.*2.0-qd2_3.*t2.*t6.*t9.*t14.*t19.*2.0-qd1_2.*t3.*t7.*t13.*t15.*t28.*2.0-qd1_2.*t4.*t6.*t13.*t19.*t28.*2.0+qd2_2.*t4.*t6.*t9.*t13.*t28.*2.0+qd2_3.*t3.*t4.*t13.*t15.*t28.*2.0-qd2_2.*t3.*t4.*t11.*t19.*t28.*2.0+qd2_2.*t6.*t7.*t11.*t15.*t28.*2.0+qd2_3.*t6.*t7.*t13.*t19.*t28.*2.0+t2.*t6.*t9.*t11.*t15.*t28.*3.0+qd1_2.*t2.*t3.*t9.*t13.*t19.*t28.*2.0+qd2_3.*t2.*t6.*t9.*t13.*t15.*t28.*2.0-qd2_2.*t2.*t6.*t9.*t11.*t19.*t28.*2.0)+(t51.*(t27.*t36+t24.*t40+t12.*t28.*t42).*2.0+t78.*(t24.*t56-t27.*t65+t12.*t28.*t67).*2.0).*(t12.*t28.*(3.0./2.0)-qd2_2.*t10.*t11+qd2_2.*t12.*t13.*t14)-(t51.*(t36.*t61+t40.*t59-t10.*t28.*t42).*2.0-t78.*(-t56.*t59+t61.*t65+t10.*t28.*t67).*2.0).*(t10.*t28.*(3.0./2.0)+qd2_2.*t11.*t12+qd2_2.*t10.*t13.*t14)+t2.*t3.*(2.5e1./1.6e1)-t9.*t19.*(9.0./8.0)+t27.*(t24.*t83.*(1.0./1.0e1)+t27.*t85.*(1.0./1.0e1)+t12.*t28.*t89.*(1.0./1.0e1))+t61.*(t59.*t83.*(1.0./1.0e1)+t61.*t85.*(1.0./1.0e1)-t10.*t28.*t89.*(1.0./1.0e1))-t11.*t28.*(t14.*t89.*(1.0./1.0e1)-t11.*t28.*t85.*(1.0./1.0e1)+t13.*t28.*t83.*(1.0./1.0e1))-t4.*t6.*t9.*(2.5e1./1.6e1)+t2.*t7.*t15.*(9.0./8.0)-qd1_2.*t4.*t6.*t15.*(3.0./4.0)+qd1_2.*t3.*t7.*t19.*(3.0./4.0)+qd1_2.*t2.*t3.*t9.*t15.*(3.0./4.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[6,6]);