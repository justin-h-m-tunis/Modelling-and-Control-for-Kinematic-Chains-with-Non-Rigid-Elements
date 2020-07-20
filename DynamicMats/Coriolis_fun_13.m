function out1 = Coriolis_fun_13(in1,in2,in3,in4)
%CORIOLIS_FUN_13
%    OUT1 = CORIOLIS_FUN_13(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    14-May-2020 18:16:36

dqt1 = in3(1,:);
dqt2 = in3(2,:);
qd1_2 = in2(3);
qd1_4 = in2(7);
qd1_5 = in2(9);
qd1_6 = in2(11);
qd2_2 = in2(4);
qd2_3 = in2(6);
qd2_4 = in2(8);
qd2_5 = in2(10);
qd2_6 = in2(12);
qt2 = in1(2,:);
t2 = cos(qt2);
t3 = cos(qd1_5);
t4 = cos(qd2_6);
t5 = sin(qd1_4);
t6 = sin(qd2_4);
t7 = cos(qd1_6);
t8 = cos(qd2_4);
t9 = sin(qd1_5);
t10 = sin(qd2_5);
t11 = sin(qt2);
t12 = cos(qd1_4);
t13 = sin(qd1_6);
t14 = cos(qd2_5);
t15 = sin(qd2_6);
t16 = t2.*t13.*(3.0./2.0);
t17 = t7.*t11.*t12.*(3.0./2.0);
t18 = qd1_2.*t3.*t11.*t12.*t13;
t19 = t16+t17+t18-qd1_2.*t2.*t3.*t7-qd1_2.*t5.*t9.*t11;
t20 = t2.*t7.*t12.*(3.0./2.0);
t21 = qd1_2.*t3.*t7.*t11;
t22 = qd1_2.*t2.*t3.*t12.*t13;
t23 = t2.*t4.*t6.*t13.*(3.0./2.0);
t24 = t2.*t7.*t8.*t10.*t12.*(3.0./2.0);
t25 = t3.*t5.*t8.*t11.*t14.*(3.0./2.0);
t26 = t4.*t6.*t7.*t11.*t12.*(3.0./2.0);
t27 = qd2_2.*t2.*t3.*t5.*t8.*t14;
t28 = qd1_2.*t3.*t7.*t8.*t10.*t11;
t29 = qd2_2.*t7.*t8.*t9.*t11.*t14;
t30 = qd2_3.*t2.*t3.*t4.*t5.*t6;
t31 = qd2_3.*t2.*t7.*t8.*t9.*t10;
t32 = qd2_3.*t4.*t6.*t7.*t9.*t11;
t33 = t2.*t6.*t7.*t12.*t14.*t15.*(3.0./2.0);
t34 = t8.*t9.*t11.*t12.*t13.*t14.*(3.0./2.0);
t35 = t2.*t6.*t7.*t9.*t10.*t15.*(3.0./2.0);
t36 = qd1_2.*t2.*t3.*t8.*t10.*t12.*t13;
t37 = qd2_2.*t2.*t8.*t9.*t12.*t13.*t14;
t38 = qd1_2.*t3.*t4.*t6.*t11.*t12.*t13;
t39 = qd1_2.*t3.*t6.*t7.*t11.*t14.*t15;
t40 = qd2_3.*t2.*t4.*t6.*t9.*t12.*t13;
t41 = qd2_3.*t2.*t6.*t7.*t9.*t14.*t15;
t42 = qd1_2.*t2.*t3.*t6.*t12.*t13.*t14.*t15;
t43 = t23+t24+t25+t26+t27+t28+t29+t30+t31+t32+t33+t34+t35+t36+t37+t38+t39+t40+t41+t42-t8.*t10.*t11.*t13.*(3.0./2.0)-t2.*t7.*t8.*t9.*t14.*(3.0./2.0)-t6.*t11.*t13.*t14.*t15.*(3.0./2.0)-qd1_2.*t2.*t3.*t4.*t6.*t7-qd1_2.*t2.*t5.*t8.*t9.*t10-qd1_2.*t4.*t5.*t6.*t9.*t11-qd2_3.*t3.*t5.*t8.*t10.*t11-t3.*t5.*t6.*t10.*t11.*t15.*(3.0./2.0)-qd1_2.*t2.*t5.*t6.*t9.*t14.*t15-qd2_2.*t2.*t3.*t5.*t6.*t10.*t15-qd2_3.*t3.*t5.*t6.*t11.*t14.*t15-qd2_2.*t6.*t7.*t9.*t10.*t11.*t15-qd2_3.*t8.*t9.*t10.*t11.*t12.*t13-t6.*t9.*t10.*t11.*t12.*t13.*t15.*(3.0./2.0)-qd2_2.*t2.*t6.*t9.*t10.*t12.*t13.*t15-qd2_3.*t6.*t9.*t11.*t12.*t13.*t14.*t15;
t44 = t3.*t4.*t5.*t8.*t11;
t45 = t2.*t3.*t5.*t8.*t14.*t15;
t46 = t4.*t8.*t9.*t11.*t12.*t13;
t47 = t7.*t8.*t9.*t11.*t14.*t15;
t48 = t2.*t8.*t9.*t12.*t13.*t14.*t15;
t49 = t44+t45+t46+t47+t48-t2.*t3.*t5.*t6.*t10-t2.*t4.*t7.*t8.*t9-t6.*t7.*t9.*t10.*t11-t2.*t6.*t9.*t10.*t12.*t13;
t50 = t2.*t6.*t7.*t10.*t12.*(3.0./2.0);
t51 = t3.*t5.*t6.*t11.*t14.*(3.0./2.0);
t52 = t8.*t11.*t13.*t14.*t15.*(3.0./2.0);
t53 = qd2_2.*t2.*t3.*t5.*t6.*t14;
t54 = qd1_2.*t4.*t5.*t8.*t9.*t11;
t55 = qd1_2.*t3.*t6.*t7.*t10.*t11;
t56 = qd2_2.*t6.*t7.*t9.*t11.*t14;
t57 = qd2_3.*t2.*t6.*t7.*t9.*t10;
t58 = t6.*t9.*t11.*t12.*t13.*t14.*(3.0./2.0);
t59 = t3.*t5.*t8.*t10.*t11.*t15.*(3.0./2.0);
t60 = qd1_2.*t2.*t3.*t4.*t7.*t8;
t61 = qd2_2.*t7.*t8.*t9.*t10.*t11.*t15;
t62 = t8.*t9.*t10.*t11.*t12.*t13.*t15.*(3.0./2.0);
t63 = qd1_2.*t2.*t3.*t6.*t10.*t12.*t13;
t64 = qd1_2.*t2.*t5.*t8.*t9.*t14.*t15;
t65 = qd2_2.*t2.*t6.*t9.*t12.*t13.*t14;
t66 = qd2_2.*t2.*t3.*t5.*t8.*t10.*t15;
t67 = qd2_3.*t3.*t5.*t8.*t11.*t14.*t15;
t68 = qd2_2.*t2.*t8.*t9.*t10.*t12.*t13.*t15;
t69 = qd2_3.*t8.*t9.*t11.*t12.*t13.*t14.*t15;
t70 = t50+t51+t52+t53+t54+t55+t56+t57+t58+t59+t60+t61+t62+t63+t64+t65+t66+t67+t68+t69-t2.*t4.*t8.*t13.*(3.0./2.0)-t6.*t10.*t11.*t13.*(3.0./2.0)-t2.*t6.*t7.*t9.*t14.*(3.0./2.0)-t4.*t7.*t8.*t11.*t12.*(3.0./2.0)-qd1_2.*t2.*t5.*t6.*t9.*t10-qd2_3.*t2.*t3.*t4.*t5.*t8-qd2_3.*t3.*t5.*t6.*t10.*t11-qd2_3.*t4.*t7.*t8.*t9.*t11-t2.*t7.*t8.*t9.*t10.*t15.*(3.0./2.0)-t2.*t7.*t8.*t12.*t14.*t15.*(3.0./2.0)-qd1_2.*t3.*t4.*t8.*t11.*t12.*t13-qd1_2.*t3.*t7.*t8.*t11.*t14.*t15-qd2_3.*t2.*t4.*t8.*t9.*t12.*t13-qd2_3.*t2.*t7.*t8.*t9.*t14.*t15-qd2_3.*t6.*t9.*t10.*t11.*t12.*t13-qd1_2.*t2.*t3.*t8.*t12.*t13.*t14.*t15;
t71 = t2.*t3.*t5.*t8.*t10;
t72 = t3.*t4.*t5.*t6.*t11;
t73 = t7.*t8.*t9.*t10.*t11;
t74 = t2.*t8.*t9.*t10.*t12.*t13;
t75 = t2.*t3.*t5.*t6.*t14.*t15;
t76 = t4.*t6.*t9.*t11.*t12.*t13;
t77 = t6.*t7.*t9.*t11.*t14.*t15;
t78 = t2.*t6.*t9.*t12.*t13.*t14.*t15;
t79 = t71+t72+t73+t74+t75+t76+t77+t78-t2.*t4.*t6.*t7.*t9;
t80 = t2.*t13.*t15.*(3.0./2.0);
t81 = t7.*t11.*t12.*t15.*(3.0./2.0);
t82 = t4.*t11.*t13.*t14.*(3.0./2.0);
t83 = qd2_3.*t2.*t3.*t5.*t15;
t84 = qd2_3.*t7.*t9.*t11.*t15;
t85 = t3.*t4.*t5.*t10.*t11.*(3.0./2.0);
t86 = qd1_2.*t2.*t4.*t5.*t9.*t14;
t87 = qd2_2.*t2.*t3.*t4.*t5.*t10;
t88 = qd2_3.*t3.*t4.*t5.*t11.*t14;
t89 = qd1_2.*t3.*t11.*t12.*t13.*t15;
t90 = qd2_3.*t2.*t9.*t12.*t13.*t15;
t91 = qd2_2.*t4.*t7.*t9.*t10.*t11;
t92 = t4.*t9.*t10.*t11.*t12.*t13.*(3.0./2.0);
t93 = qd2_2.*t2.*t4.*t9.*t10.*t12.*t13;
t94 = qd2_3.*t4.*t9.*t11.*t12.*t13.*t14;
t95 = t80+t81+t82+t83+t84+t85+t86+t87+t88+t89+t90+t91+t92+t93+t94-qd1_2.*t2.*t3.*t7.*t15-qd1_2.*t5.*t9.*t11.*t15-t2.*t4.*t7.*t9.*t10.*(3.0./2.0)-t2.*t4.*t7.*t12.*t14.*(3.0./2.0)-qd1_2.*t3.*t4.*t7.*t11.*t14-qd2_3.*t2.*t4.*t7.*t9.*t14-qd1_2.*t2.*t3.*t4.*t12.*t13.*t14;
t96 = t2.*t7.*t9.*t15;
t97 = t2.*t3.*t4.*t5.*t14;
t98 = t4.*t7.*t9.*t11.*t14;
t99 = t2.*t4.*t9.*t12.*t13.*t14;
t100 = t96+t97+t98+t99-t3.*t5.*t11.*t15-t9.*t11.*t12.*t13.*t15;
out1 = reshape([-dqt1.*(-t19.*(t20+t21+t22-t11.*t13.*(3.0./2.0)-qd1_2.*t2.*t5.*t9)+t19.*(t20+t21+t22+t3.*t12.*(3.0./4.0)-t11.*t13.*(3.0./2.0)-t5.*t9.*t13.*(3.0./4.0)-qd1_2.*t2.*t5.*t9)+t49.*(t3.*t6.*t12.*t14.*(1.0./1.0e1)+t2.*t3.*t4.*t5.*t8.*(1.0./1.0e1)-t2.*t6.*t7.*t9.*t10.*(1.0./1.0e1)+t3.*t5.*t6.*t10.*t11.*(1.0./1.0e1)+t4.*t7.*t8.*t9.*t11.*(1.0./1.0e1)-t5.*t6.*t9.*t13.*t14.*(1.0./1.0e1)+t3.*t8.*t10.*t12.*t15.*(1.0./1.0e1)+t2.*t4.*t8.*t9.*t12.*t13.*(1.0./1.0e1)+t2.*t7.*t8.*t9.*t14.*t15.*(1.0./1.0e1)-t3.*t5.*t8.*t11.*t14.*t15.*(1.0./1.0e1)-t5.*t8.*t9.*t10.*t13.*t15.*(1.0./1.0e1)+t6.*t9.*t10.*t11.*t12.*t13.*(1.0./1.0e1)-t8.*t9.*t11.*t12.*t13.*t14.*t15.*(1.0./1.0e1)).*(1.0./2.0)+t79.*(t3.*t8.*t12.*t14.*(-1.0./1.0e1)+t2.*t3.*t4.*t5.*t6.*(1.0./1.0e1)+t2.*t7.*t8.*t9.*t10.*(1.0./1.0e1)-t3.*t5.*t8.*t10.*t11.*(1.0./1.0e1)+t4.*t6.*t7.*t9.*t11.*(1.0./1.0e1)+t3.*t6.*t10.*t12.*t15.*(1.0./1.0e1)+t5.*t8.*t9.*t13.*t14.*(1.0./1.0e1)+t2.*t4.*t6.*t9.*t12.*t13.*(1.0./1.0e1)+t2.*t6.*t7.*t9.*t14.*t15.*(1.0./1.0e1)-t3.*t5.*t6.*t11.*t14.*t15.*(1.0./1.0e1)-t5.*t6.*t9.*t10.*t13.*t15.*(1.0./1.0e1)-t8.*t9.*t10.*t11.*t12.*t13.*(1.0./1.0e1)-t6.*t9.*t11.*t12.*t13.*t14.*t15.*(1.0./1.0e1)).*(1.0./2.0)+(t2.*t13.*t15.*3.0+t4.*t11.*t13.*t14.*3.0+t7.*t11.*t12.*t15.*3.0-qd1_2.*t2.*t3.*t7.*t15.*2.0+qd2_3.*t2.*t3.*t5.*t15.*2.0-qd1_2.*t5.*t9.*t11.*t15.*2.0+qd2_3.*t7.*t9.*t11.*t15.*2.0-t2.*t4.*t7.*t9.*t10.*3.0+t3.*t4.*t5.*t10.*t11.*3.0-t2.*t4.*t7.*t12.*t14.*3.0+qd1_2.*t2.*t4.*t5.*t9.*t14.*2.0+qd2_2.*t2.*t3.*t4.*t5.*t10.*2.0-qd1_2.*t3.*t4.*t7.*t11.*t14.*2.0-qd2_3.*t2.*t4.*t7.*t9.*t14.*2.0+qd2_3.*t3.*t4.*t5.*t11.*t14.*2.0+qd2_2.*t4.*t7.*t9.*t10.*t11.*2.0+qd1_2.*t3.*t11.*t12.*t13.*t15.*2.0+qd2_3.*t2.*t9.*t12.*t13.*t15.*2.0+t4.*t9.*t10.*t11.*t12.*t13.*3.0-qd1_2.*t2.*t3.*t4.*t12.*t13.*t14.*2.0+qd2_2.*t2.*t4.*t9.*t10.*t12.*t13.*2.0+qd2_3.*t4.*t9.*t11.*t12.*t13.*t14.*2.0).*(t3.*t12.*t15.*(3.0./2.0)-t11.*t13.*t15.*(3.0./2.0)+t4.*t5.*t7.*t10.*(3.0./2.0)+t2.*t4.*t13.*t14.*(3.0./2.0)+t2.*t7.*t12.*t15.*(3.0./2.0)-t5.*t9.*t13.*t15.*(3.0./2.0)-qd1_2.*t2.*t5.*t9.*t15+qd1_2.*t4.*t9.*t10.*t12+qd1_2.*t3.*t7.*t11.*t15-qd2_2.*t3.*t4.*t12.*t14+qd2_3.*t2.*t7.*t9.*t15-qd2_3.*t3.*t5.*t11.*t15+t2.*t3.*t4.*t5.*t10.*(3.0./2.0)+t4.*t7.*t9.*t10.*t11.*(3.0./2.0)+t4.*t7.*t11.*t12.*t14.*(3.0./2.0)-qd1_2.*t2.*t3.*t4.*t7.*t14+qd1_2.*t3.*t4.*t5.*t10.*t13+qd2_3.*t2.*t3.*t4.*t5.*t14+qd2_2.*t2.*t4.*t7.*t9.*t10-qd1_2.*t4.*t5.*t9.*t11.*t14-qd2_2.*t3.*t4.*t5.*t10.*t11+qd1_2.*t2.*t3.*t12.*t13.*t15+qd2_2.*t4.*t5.*t9.*t13.*t14+qd2_3.*t4.*t7.*t9.*t11.*t14-qd2_3.*t9.*t11.*t12.*t13.*t15+t2.*t4.*t9.*t10.*t12.*t13.*(3.0./2.0)+qd1_2.*t3.*t4.*t11.*t12.*t13.*t14+qd2_3.*t2.*t4.*t9.*t12.*t13.*t14-qd2_2.*t4.*t9.*t10.*t11.*t12.*t13).*(1.0./2.0)-t100.*(t2.*t3.*t5.*t15.*(1.0./1.0e1)-t3.*t4.*t10.*t12.*(1.0./1.0e1)+t7.*t9.*t11.*t15.*(1.0./1.0e1)-t2.*t4.*t7.*t9.*t14.*(1.0./1.0e1)+t3.*t4.*t5.*t11.*t14.*(1.0./1.0e1)+t4.*t5.*t9.*t10.*t13.*(1.0./1.0e1)+t2.*t9.*t12.*t13.*t15.*(1.0./1.0e1)+t4.*t9.*t11.*t12.*t13.*t14.*(1.0./1.0e1)).*(1.0./2.0)+(t2.*t7.*t9.*(-3.0./4.0)+t3.*t5.*t11.*(3.0./4.0)+t9.*t11.*t12.*t13.*(3.0./4.0)).*(t5.*t7.*(3.0./2.0)+qd1_2.*t9.*t12+t2.*t3.*t5.*(3.0./4.0)+t7.*t9.*t11.*(3.0./4.0)+qd1_2.*t3.*t5.*t13+t2.*t9.*t12.*t13.*(3.0./4.0))+t95.*(t3.*t12.*t15.*3.0-t11.*t13.*t15.*3.0+t4.*t5.*t7.*t10.*3.0+t2.*t4.*t13.*t14.*3.0+t2.*t7.*t12.*t15.*3.0-t5.*t9.*t13.*t15.*3.0-qd1_2.*t2.*t5.*t9.*t15.*2.0+qd1_2.*t4.*t9.*t10.*t12.*2.0+qd1_2.*t3.*t7.*t11.*t15.*2.0-qd2_2.*t3.*t4.*t12.*t14.*2.0+qd2_3.*t2.*t7.*t9.*t15.*2.0-qd2_3.*t3.*t5.*t11.*t15.*2.0+t2.*t3.*t4.*t5.*t10.*3.0+t4.*t7.*t9.*t10.*t11.*3.0+t4.*t7.*t11.*t12.*t14.*3.0-qd1_2.*t2.*t3.*t4.*t7.*t14.*2.0+qd1_2.*t3.*t4.*t5.*t10.*t13.*2.0+qd2_3.*t2.*t3.*t4.*t5.*t14.*2.0+qd2_2.*t2.*t4.*t7.*t9.*t10.*2.0-qd1_2.*t4.*t5.*t9.*t11.*t14.*2.0-qd2_2.*t3.*t4.*t5.*t10.*t11.*2.0+qd1_2.*t2.*t3.*t12.*t13.*t15.*2.0+qd2_2.*t4.*t5.*t9.*t13.*t14.*2.0+qd2_3.*t4.*t7.*t9.*t11.*t14.*2.0-qd2_3.*t9.*t11.*t12.*t13.*t15.*2.0+t2.*t4.*t9.*t10.*t12.*t13.*3.0+qd1_2.*t3.*t4.*t11.*t12.*t13.*t14.*2.0+qd2_3.*t2.*t4.*t9.*t12.*t13.*t14.*2.0-qd2_2.*t4.*t9.*t10.*t11.*t12.*t13.*2.0).*(1.0./2.0)+(t2.*t3.*t5.*t8.*t10.*(1.0./1.0e1)-t2.*t4.*t6.*t7.*t9.*(1.0./1.0e1)+t3.*t4.*t5.*t6.*t11.*(1.0./1.0e1)+t7.*t8.*t9.*t10.*t11.*(1.0./1.0e1)+t2.*t3.*t5.*t6.*t14.*t15.*(1.0./1.0e1)+t2.*t8.*t9.*t10.*t12.*t13.*(1.0./1.0e1)+t4.*t6.*t9.*t11.*t12.*t13.*(1.0./1.0e1)+t6.*t7.*t9.*t11.*t14.*t15.*(1.0./1.0e1)+t2.*t6.*t9.*t12.*t13.*t14.*t15.*(1.0./1.0e1)).*(-t3.*t8.*t12.*t14+t2.*t3.*t4.*t5.*t6+t2.*t7.*t8.*t9.*t10-t3.*t5.*t8.*t10.*t11+t4.*t6.*t7.*t9.*t11+t3.*t6.*t10.*t12.*t15+t5.*t8.*t9.*t13.*t14+t2.*t4.*t6.*t9.*t12.*t13+t2.*t6.*t7.*t9.*t14.*t15-t3.*t5.*t6.*t11.*t14.*t15-t5.*t6.*t9.*t10.*t13.*t15-t8.*t9.*t10.*t11.*t12.*t13-t6.*t9.*t11.*t12.*t13.*t14.*t15).*(1.0./2.0)+(t2.*t3.*t5.*t6.*t10.*(-1.0./1.0e1)-t2.*t4.*t7.*t8.*t9.*(1.0./1.0e1)+t3.*t4.*t5.*t8.*t11.*(1.0./1.0e1)-t6.*t7.*t9.*t10.*t11.*(1.0./1.0e1)+t2.*t3.*t5.*t8.*t14.*t15.*(1.0./1.0e1)-t2.*t6.*t9.*t10.*t12.*t13.*(1.0./1.0e1)+t4.*t8.*t9.*t11.*t12.*t13.*(1.0./1.0e1)+t7.*t8.*t9.*t11.*t14.*t15.*(1.0./1.0e1)+t2.*t8.*t9.*t12.*t13.*t14.*t15.*(1.0./1.0e1)).*(t3.*t6.*t12.*t14+t2.*t3.*t4.*t5.*t8-t2.*t6.*t7.*t9.*t10+t3.*t5.*t6.*t10.*t11+t4.*t7.*t8.*t9.*t11-t5.*t6.*t9.*t13.*t14+t3.*t8.*t10.*t12.*t15+t2.*t4.*t8.*t9.*t12.*t13+t2.*t7.*t8.*t9.*t14.*t15-t3.*t5.*t8.*t11.*t14.*t15-t5.*t8.*t9.*t10.*t13.*t15+t6.*t9.*t10.*t11.*t12.*t13-t8.*t9.*t11.*t12.*t13.*t14.*t15).*(1.0./2.0)-t43.*(t3.*t4.*t6.*t12.*-3.0+t2.*t8.*t10.*t13.*3.0+t4.*t6.*t11.*t13.*3.0-t5.*t7.*t8.*t14.*3.0-qd1_2.*t8.*t9.*t12.*t14.*2.0-qd2_2.*t3.*t8.*t10.*t12.*2.0-t2.*t4.*t6.*t7.*t12.*3.0-t2.*t3.*t5.*t8.*t14.*3.0+t4.*t5.*t6.*t9.*t13.*3.0+t5.*t6.*t7.*t10.*t15.*3.0+t7.*t8.*t10.*t11.*t12.*3.0-t7.*t8.*t9.*t11.*t14.*3.0+t2.*t6.*t13.*t14.*t15.*3.0+qd1_2.*t2.*t4.*t5.*t6.*t9.*2.0-qd1_2.*t2.*t3.*t7.*t8.*t10.*2.0-qd1_2.*t3.*t4.*t6.*t7.*t11.*2.0+qd2_3.*t2.*t3.*t5.*t8.*t10.*2.0-qd2_3.*t2.*t4.*t6.*t7.*t9.*2.0+qd2_3.*t3.*t4.*t5.*t6.*t11.*2.0-qd1_2.*t3.*t5.*t8.*t13.*t14.*2.0-qd1_2.*t5.*t8.*t9.*t10.*t11.*2.0-qd2_2.*t2.*t7.*t8.*t9.*t14.*2.0+qd2_2.*t3.*t5.*t8.*t11.*t14.*2.0+qd1_2.*t6.*t9.*t10.*t12.*t15.*2.0+qd2_2.*t5.*t8.*t9.*t10.*t13.*2.0+qd2_3.*t7.*t8.*t9.*t10.*t11.*2.0-qd2_2.*t3.*t6.*t12.*t14.*t15.*2.0+t2.*t3.*t5.*t6.*t10.*t15.*3.0-t2.*t8.*t9.*t12.*t13.*t14.*3.0+t6.*t7.*t9.*t10.*t11.*t15.*3.0+t6.*t7.*t11.*t12.*t14.*t15.*3.0-qd1_2.*t2.*t3.*t4.*t6.*t12.*t13.*2.0-qd1_2.*t2.*t3.*t6.*t7.*t14.*t15.*2.0+qd1_2.*t3.*t5.*t6.*t10.*t13.*t15.*2.0+qd2_3.*t2.*t3.*t5.*t6.*t14.*t15.*2.0+qd1_2.*t3.*t8.*t10.*t11.*t12.*t13.*2.0+qd2_2.*t2.*t6.*t7.*t9.*t10.*t15.*2.0-qd1_2.*t5.*t6.*t9.*t11.*t14.*t15.*2.0-qd2_2.*t3.*t5.*t6.*t10.*t11.*t15.*2.0+qd2_3.*t2.*t8.*t9.*t10.*t12.*t13.*2.0+qd2_3.*t4.*t6.*t9.*t11.*t12.*t13.*2.0+qd2_2.*t5.*t6.*t9.*t13.*t14.*t15.*2.0+qd2_3.*t6.*t7.*t9.*t11.*t14.*t15.*2.0+qd2_2.*t8.*t9.*t11.*t12.*t13.*t14.*2.0+t2.*t6.*t9.*t10.*t12.*t13.*t15.*3.0+qd1_2.*t3.*t6.*t11.*t12.*t13.*t14.*t15.*2.0+qd2_3.*t2.*t6.*t9.*t12.*t13.*t14.*t15.*2.0-qd2_2.*t6.*t9.*t10.*t11.*t12.*t13.*t15.*2.0).*(1.0./2.0)+t70.*(t3.*t4.*t8.*t12.*-3.0-t2.*t6.*t10.*t13.*3.0+t5.*t6.*t7.*t14.*3.0+t4.*t8.*t11.*t13.*3.0+qd1_2.*t6.*t9.*t12.*t14.*2.0+qd2_2.*t3.*t6.*t10.*t12.*2.0+t2.*t3.*t5.*t6.*t14.*3.0-t2.*t4.*t7.*t8.*t12.*3.0+t4.*t5.*t8.*t9.*t13.*3.0+t5.*t7.*t8.*t10.*t15.*3.0-t6.*t7.*t10.*t11.*t12.*3.0+t6.*t7.*t9.*t11.*t14.*3.0+t2.*t8.*t13.*t14.*t15.*3.0+qd1_2.*t2.*t3.*t6.*t7.*t10.*2.0+qd1_2.*t2.*t4.*t5.*t8.*t9.*2.0-qd1_2.*t3.*t4.*t7.*t8.*t11.*2.0-qd2_3.*t2.*t3.*t5.*t6.*t10.*2.0+qd1_2.*t3.*t5.*t6.*t13.*t14.*2.0+qd1_2.*t5.*t6.*t9.*t10.*t11.*2.0-qd2_3.*t2.*t4.*t7.*t8.*t9.*2.0+qd2_3.*t3.*t4.*t5.*t8.*t11.*2.0+qd2_2.*t2.*t6.*t7.*t9.*t14.*2.0-qd2_2.*t3.*t5.*t6.*t11.*t14.*2.0-qd2_2.*t5.*t6.*t9.*t10.*t13.*2.0+qd1_2.*t8.*t9.*t10.*t12.*t15.*2.0-qd2_3.*t6.*t7.*t9.*t10.*t11.*2.0-qd2_2.*t3.*t8.*t12.*t14.*t15.*2.0+t2.*t3.*t5.*t8.*t10.*t15.*3.0+t2.*t6.*t9.*t12.*t13.*t14.*3.0+t7.*t8.*t9.*t10.*t11.*t15.*3.0+t7.*t8.*t11.*t12.*t14.*t15.*3.0-qd1_2.*t2.*t3.*t4.*t8.*t12.*t13.*2.0-qd1_2.*t2.*t3.*t7.*t8.*t14.*t15.*2.0+qd1_2.*t3.*t5.*t8.*t10.*t13.*t15.*2.0-qd1_2.*t3.*t6.*t10.*t11.*t12.*t13.*2.0+qd2_3.*t2.*t3.*t5.*t8.*t14.*t15.*2.0+qd2_2.*t2.*t7.*t8.*t9.*t10.*t15.*2.0-qd1_2.*t5.*t8.*t9.*t11.*t14.*t15.*2.0-qd2_2.*t3.*t5.*t8.*t10.*t11.*t15.*2.0-qd2_3.*t2.*t6.*t9.*t10.*t12.*t13.*2.0+qd2_3.*t4.*t8.*t9.*t11.*t12.*t13.*2.0+qd2_2.*t5.*t8.*t9.*t13.*t14.*t15.*2.0-qd2_2.*t6.*t9.*t11.*t12.*t13.*t14.*2.0+qd2_3.*t7.*t8.*t9.*t11.*t14.*t15.*2.0+t2.*t8.*t9.*t10.*t12.*t13.*t15.*3.0+qd1_2.*t3.*t8.*t11.*t12.*t13.*t14.*t15.*2.0+qd2_3.*t2.*t8.*t9.*t12.*t13.*t14.*t15.*2.0-qd2_2.*t8.*t9.*t10.*t11.*t12.*t13.*t15.*2.0).*(1.0./2.0)-(t2.*t4.*t6.*t13.*3.0-t8.*t10.*t11.*t13.*3.0+t2.*t7.*t8.*t10.*t12.*3.0-t2.*t7.*t8.*t9.*t14.*3.0+t4.*t6.*t7.*t11.*t12.*3.0+t3.*t5.*t8.*t11.*t14.*3.0-t6.*t11.*t13.*t14.*t15.*3.0-qd1_2.*t2.*t3.*t4.*t6.*t7.*2.0+qd2_3.*t2.*t3.*t4.*t5.*t6.*2.0-qd1_2.*t2.*t5.*t8.*t9.*t10.*2.0-qd1_2.*t4.*t5.*t6.*t9.*t11.*2.0+qd1_2.*t3.*t7.*t8.*t10.*t11.*2.0+qd2_2.*t2.*t3.*t5.*t8.*t14.*2.0+qd2_3.*t2.*t7.*t8.*t9.*t10.*2.0-qd2_3.*t3.*t5.*t8.*t10.*t11.*2.0+qd2_3.*t4.*t6.*t7.*t9.*t11.*2.0+qd2_2.*t7.*t8.*t9.*t11.*t14.*2.0+t2.*t6.*t7.*t9.*t10.*t15.*3.0-t3.*t5.*t6.*t10.*t11.*t15.*3.0+t2.*t6.*t7.*t12.*t14.*t15.*3.0+t8.*t9.*t11.*t12.*t13.*t14.*3.0+qd1_2.*t2.*t3.*t8.*t10.*t12.*t13.*2.0+qd1_2.*t3.*t4.*t6.*t11.*t12.*t13.*2.0-qd1_2.*t2.*t5.*t6.*t9.*t14.*t15.*2.0-qd2_2.*t2.*t3.*t5.*t6.*t10.*t15.*2.0+qd1_2.*t3.*t6.*t7.*t11.*t14.*t15.*2.0+qd2_3.*t2.*t4.*t6.*t9.*t12.*t13.*2.0+qd2_3.*t2.*t6.*t7.*t9.*t14.*t15.*2.0-qd2_3.*t3.*t5.*t6.*t11.*t14.*t15.*2.0+qd2_2.*t2.*t8.*t9.*t12.*t13.*t14.*2.0-qd2_2.*t6.*t7.*t9.*t10.*t11.*t15.*2.0-qd2_3.*t8.*t9.*t10.*t11.*t12.*t13.*2.0-t6.*t9.*t10.*t11.*t12.*t13.*t15.*3.0+qd1_2.*t2.*t3.*t6.*t12.*t13.*t14.*t15.*2.0-qd2_2.*t2.*t6.*t9.*t10.*t12.*t13.*t15.*2.0-qd2_3.*t6.*t9.*t11.*t12.*t13.*t14.*t15.*2.0).*(t3.*t4.*t6.*t12.*(-3.0./2.0)+t2.*t8.*t10.*t13.*(3.0./2.0)+t4.*t6.*t11.*t13.*(3.0./2.0)-t5.*t7.*t8.*t14.*(3.0./2.0)-qd1_2.*t8.*t9.*t12.*t14-qd2_2.*t3.*t8.*t10.*t12-t2.*t4.*t6.*t7.*t12.*(3.0./2.0)-t2.*t3.*t5.*t8.*t14.*(3.0./2.0)+t4.*t5.*t6.*t9.*t13.*(3.0./2.0)+t5.*t6.*t7.*t10.*t15.*(3.0./2.0)+t7.*t8.*t10.*t11.*t12.*(3.0./2.0)-t7.*t8.*t9.*t11.*t14.*(3.0./2.0)+t2.*t6.*t13.*t14.*t15.*(3.0./2.0)+qd1_2.*t2.*t4.*t5.*t6.*t9-qd1_2.*t2.*t3.*t7.*t8.*t10-qd1_2.*t3.*t4.*t6.*t7.*t11+qd2_3.*t2.*t3.*t5.*t8.*t10-qd2_3.*t2.*t4.*t6.*t7.*t9+qd2_3.*t3.*t4.*t5.*t6.*t11-qd1_2.*t3.*t5.*t8.*t13.*t14-qd1_2.*t5.*t8.*t9.*t10.*t11-qd2_2.*t2.*t7.*t8.*t9.*t14+qd2_2.*t3.*t5.*t8.*t11.*t14+qd1_2.*t6.*t9.*t10.*t12.*t15+qd2_2.*t5.*t8.*t9.*t10.*t13+qd2_3.*t7.*t8.*t9.*t10.*t11-qd2_2.*t3.*t6.*t12.*t14.*t15+t2.*t3.*t5.*t6.*t10.*t15.*(3.0./2.0)-t2.*t8.*t9.*t12.*t13.*t14.*(3.0./2.0)+t6.*t7.*t9.*t10.*t11.*t15.*(3.0./2.0)+t6.*t7.*t11.*t12.*t14.*t15.*(3.0./2.0)-qd1_2.*t2.*t3.*t4.*t6.*t12.*t13-qd1_2.*t2.*t3.*t6.*t7.*t14.*t15+qd1_2.*t3.*t5.*t6.*t10.*t13.*t15+qd2_3.*t2.*t3.*t5.*t6.*t14.*t15+qd1_2.*t3.*t8.*t10.*t11.*t12.*t13+qd2_2.*t2.*t6.*t7.*t9.*t10.*t15-qd1_2.*t5.*t6.*t9.*t11.*t14.*t15-qd2_2.*t3.*t5.*t6.*t10.*t11.*t15+qd2_3.*t2.*t8.*t9.*t10.*t12.*t13+qd2_3.*t4.*t6.*t9.*t11.*t12.*t13+qd2_2.*t5.*t6.*t9.*t13.*t14.*t15+qd2_3.*t6.*t7.*t9.*t11.*t14.*t15+qd2_2.*t8.*t9.*t11.*t12.*t13.*t14+t2.*t6.*t9.*t10.*t12.*t13.*t15.*(3.0./2.0)+qd1_2.*t3.*t6.*t11.*t12.*t13.*t14.*t15+qd2_3.*t2.*t6.*t9.*t12.*t13.*t14.*t15-qd2_2.*t6.*t9.*t10.*t11.*t12.*t13.*t15).*(1.0./2.0)+(t2.*t4.*t8.*t13.*-3.0-t6.*t10.*t11.*t13.*3.0+t2.*t6.*t7.*t10.*t12.*3.0-t2.*t6.*t7.*t9.*t14.*3.0+t3.*t5.*t6.*t11.*t14.*3.0-t4.*t7.*t8.*t11.*t12.*3.0+t8.*t11.*t13.*t14.*t15.*3.0+qd1_2.*t2.*t3.*t4.*t7.*t8.*2.0-qd1_2.*t2.*t5.*t6.*t9.*t10.*2.0-qd2_3.*t2.*t3.*t4.*t5.*t8.*2.0+qd1_2.*t3.*t6.*t7.*t10.*t11.*2.0+qd1_2.*t4.*t5.*t8.*t9.*t11.*2.0+qd2_2.*t2.*t3.*t5.*t6.*t14.*2.0+qd2_3.*t2.*t6.*t7.*t9.*t10.*2.0-qd2_3.*t3.*t5.*t6.*t10.*t11.*2.0-qd2_3.*t4.*t7.*t8.*t9.*t11.*2.0+qd2_2.*t6.*t7.*t9.*t11.*t14.*2.0-t2.*t7.*t8.*t9.*t10.*t15.*3.0+t3.*t5.*t8.*t10.*t11.*t15.*3.0-t2.*t7.*t8.*t12.*t14.*t15.*3.0+t6.*t9.*t11.*t12.*t13.*t14.*3.0+qd1_2.*t2.*t3.*t6.*t10.*t12.*t13.*2.0-qd1_2.*t3.*t4.*t8.*t11.*t12.*t13.*2.0+qd1_2.*t2.*t5.*t8.*t9.*t14.*t15.*2.0+qd2_2.*t2.*t3.*t5.*t8.*t10.*t15.*2.0-qd1_2.*t3.*t7.*t8.*t11.*t14.*t15.*2.0-qd2_3.*t2.*t4.*t8.*t9.*t12.*t13.*2.0+qd2_2.*t2.*t6.*t9.*t12.*t13.*t14.*2.0-qd2_3.*t2.*t7.*t8.*t9.*t14.*t15.*2.0+qd2_3.*t3.*t5.*t8.*t11.*t14.*t15.*2.0+qd2_2.*t7.*t8.*t9.*t10.*t11.*t15.*2.0-qd2_3.*t6.*t9.*t10.*t11.*t12.*t13.*2.0+t8.*t9.*t10.*t11.*t12.*t13.*t15.*3.0-qd1_2.*t2.*t3.*t8.*t12.*t13.*t14.*t15.*2.0+qd2_2.*t2.*t8.*t9.*t10.*t12.*t13.*t15.*2.0+qd2_3.*t8.*t9.*t11.*t12.*t13.*t14.*t15.*2.0).*(t3.*t4.*t8.*t12.*(-3.0./2.0)-t2.*t6.*t10.*t13.*(3.0./2.0)+t5.*t6.*t7.*t14.*(3.0./2.0)+t4.*t8.*t11.*t13.*(3.0./2.0)+qd1_2.*t6.*t9.*t12.*t14+qd2_2.*t3.*t6.*t10.*t12+t2.*t3.*t5.*t6.*t14.*(3.0./2.0)-t2.*t4.*t7.*t8.*t12.*(3.0./2.0)+t4.*t5.*t8.*t9.*t13.*(3.0./2.0)+t5.*t7.*t8.*t10.*t15.*(3.0./2.0)-t6.*t7.*t10.*t11.*t12.*(3.0./2.0)+t6.*t7.*t9.*t11.*t14.*(3.0./2.0)+t2.*t8.*t13.*t14.*t15.*(3.0./2.0)+qd1_2.*t2.*t3.*t6.*t7.*t10+qd1_2.*t2.*t4.*t5.*t8.*t9-qd1_2.*t3.*t4.*t7.*t8.*t11-qd2_3.*t2.*t3.*t5.*t6.*t10+qd1_2.*t3.*t5.*t6.*t13.*t14+qd1_2.*t5.*t6.*t9.*t10.*t11-qd2_3.*t2.*t4.*t7.*t8.*t9+qd2_3.*t3.*t4.*t5.*t8.*t11+qd2_2.*t2.*t6.*t7.*t9.*t14-qd2_2.*t3.*t5.*t6.*t11.*t14-qd2_2.*t5.*t6.*t9.*t10.*t13+qd1_2.*t8.*t9.*t10.*t12.*t15-qd2_3.*t6.*t7.*t9.*t10.*t11-qd2_2.*t3.*t8.*t12.*t14.*t15+t2.*t3.*t5.*t8.*t10.*t15.*(3.0./2.0)+t2.*t6.*t9.*t12.*t13.*t14.*(3.0./2.0)+t7.*t8.*t9.*t10.*t11.*t15.*(3.0./2.0)+t7.*t8.*t11.*t12.*t14.*t15.*(3.0./2.0)-qd1_2.*t2.*t3.*t4.*t8.*t12.*t13-qd1_2.*t2.*t3.*t7.*t8.*t14.*t15+qd1_2.*t3.*t5.*t8.*t10.*t13.*t15-qd1_2.*t3.*t6.*t10.*t11.*t12.*t13+qd2_3.*t2.*t3.*t5.*t8.*t14.*t15+qd2_2.*t2.*t7.*t8.*t9.*t10.*t15-qd1_2.*t5.*t8.*t9.*t11.*t14.*t15-qd2_2.*t3.*t5.*t8.*t10.*t11.*t15-qd2_3.*t2.*t6.*t9.*t10.*t12.*t13+qd2_3.*t4.*t8.*t9.*t11.*t12.*t13+qd2_2.*t5.*t8.*t9.*t13.*t14.*t15-qd2_2.*t6.*t9.*t11.*t12.*t13.*t14+qd2_3.*t7.*t8.*t9.*t11.*t14.*t15+t2.*t8.*t9.*t10.*t12.*t13.*t15.*(3.0./2.0)+qd1_2.*t3.*t8.*t11.*t12.*t13.*t14.*t15+qd2_3.*t2.*t8.*t9.*t12.*t13.*t14.*t15-qd2_2.*t8.*t9.*t10.*t11.*t12.*t13.*t15).*(1.0./2.0)-(t2.*t7.*t9.*t15.*(1.0./1.0e1)-t3.*t5.*t11.*t15.*(1.0./1.0e1)+t2.*t3.*t4.*t5.*t14.*(1.0./1.0e1)+t4.*t7.*t9.*t11.*t14.*(1.0./1.0e1)-t9.*t11.*t12.*t13.*t15.*(1.0./1.0e1)+t2.*t4.*t9.*t12.*t13.*t14.*(1.0./1.0e1)).*(t2.*t3.*t5.*t15-t3.*t4.*t10.*t12+t7.*t9.*t11.*t15-t2.*t4.*t7.*t9.*t14+t3.*t4.*t5.*t11.*t14+t4.*t5.*t9.*t10.*t13+t2.*t9.*t12.*t13.*t15+t4.*t9.*t11.*t12.*t13.*t14).*(1.0./2.0))-dqt2.*(t49.*(t6.*t14.*(1.0./1.0e1)+t8.*t10.*t15.*(1.0./1.0e1))-t79.*(t8.*t14.*(1.0./1.0e1)-t6.*t10.*t15.*(1.0./1.0e1))+t2.*t13.*(9.0./8.0)+t43.*(t4.*t6.*3.0+qd2_2.*t8.*t10.*2.0+qd2_2.*t6.*t14.*t15.*2.0)-t70.*(t4.*t8.*3.0-qd2_2.*t6.*t10.*2.0+qd2_2.*t8.*t14.*t15.*2.0)+t95.*(t15.*3.0-qd2_2.*t4.*t14.*2.0)+t7.*t11.*t12.*(9.0./8.0)+t4.*t10.*t100.*(1.0./1.0e1)-qd1_2.*t2.*t3.*t7.*(3.0./4.0)-qd1_2.*t5.*t9.*t11.*(3.0./4.0)+qd1_2.*t3.*t11.*t12.*t13.*(3.0./4.0)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[6,6]);
