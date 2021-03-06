function out1 = Mass_fun_11(in1,in2)
%MASS_FUN_11
%    OUT1 = MASS_FUN_11(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    20-Jul-2020 14:22:18

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
t2 = sin(qd1_4);
t3 = sin(qt1);
t4 = sin(qd1_5);
t5 = cos(qd1_4);
t6 = cos(qd1_6);
t7 = cos(qt1);
t8 = cos(qd1_5);
t9 = sin(qd1_6);
t10 = sin(qd2_4);
t11 = sin(qd2_5);
t12 = cos(qd2_4);
t13 = cos(qd2_5);
t14 = sin(qd2_6);
t15 = cos(qt2);
t16 = t2.*t3.*t4;
t17 = t5.*t6.*t7;
t38 = t3.*t5.*t8.*t9;
t18 = t16+t17-t38;
t19 = sin(qt2);
t20 = t7.*t9;
t21 = t3.*t6.*t8;
t22 = t20+t21;
t23 = t10.*t11;
t42 = t12.*t13.*t14;
t24 = t23-t42;
t25 = t10.*t13;
t26 = t11.*t12.*t14;
t27 = t25+t26;
t28 = cos(qd2_6);
t29 = t3.*t5.*t6;
t30 = t5.*t7.*t8.*t9;
t57 = t2.*t4.*t7;
t31 = t29+t30-t57;
t32 = t3.*t9;
t59 = t6.*t7.*t8;
t33 = t32-t59;
t34 = t3.*t4.*t5;
t35 = t2.*t3.*t8.*t9;
t99 = t2.*t6.*t7;
t36 = t34+t35-t99;
t37 = t27.*t36;
t39 = t18.*t19;
t40 = t15.*t22;
t41 = t39+t40;
t43 = t24.*t41;
t44 = t15.*t18;
t100 = t19.*t22;
t45 = t44-t100;
t46 = t12.*t28.*t45;
t47 = t37+t43+t46;
t48 = t7.*(3.0./2.0);
t49 = qd2_3.*t4.*t5.*t7;
t50 = qd2_3.*t2.*t3.*t6;
t51 = qd2_2.*t3.*t9.*t19;
t52 = t6.*t7.*t8.*t15.*(3.0./2.0);
t53 = t2.*t4.*t7.*t19.*(3.0./2.0);
t54 = qd2_3.*t2.*t7.*t8.*t9;
t55 = qd2_2.*t2.*t4.*t7.*t15;
t101 = qd1_2.*t3;
t102 = t3.*t9.*t15.*(3.0./2.0);
t103 = t3.*t5.*t6.*t19.*(3.0./2.0);
t104 = qd2_2.*t3.*t5.*t6.*t15;
t105 = qd2_2.*t6.*t7.*t8.*t19;
t106 = t5.*t7.*t8.*t9.*t19.*(3.0./2.0);
t107 = qd2_2.*t5.*t7.*t8.*t9.*t15;
t56 = t48+t49+t50+t51+t52+t53+t54+t55-t101-t102-t103-t104-t105-t106-t107;
t58 = t19.*t31;
t60 = t15.*t33;
t61 = t58+t60;
t62 = t24.*t61;
t63 = t4.*t5.*t7;
t64 = t2.*t3.*t6;
t65 = t2.*t7.*t8.*t9;
t66 = t63+t64+t65;
t67 = t15.*t31;
t108 = t19.*t33;
t68 = t67-t108;
t69 = t12.*t28.*t68;
t70 = t62+t69-t27.*t66;
t71 = t3.*(3.0./2.0);
t72 = qd1_2.*t7;
t73 = t7.*t9.*t15.*(3.0./2.0);
t74 = qd2_3.*t3.*t4.*t5;
t75 = t5.*t6.*t7.*t19.*(3.0./2.0);
t76 = t3.*t6.*t8.*t15.*(3.0./2.0);
t77 = t2.*t3.*t4.*t19.*(3.0./2.0);
t78 = qd2_2.*t5.*t6.*t7.*t15;
t79 = qd2_3.*t2.*t3.*t8.*t9;
t80 = qd2_2.*t2.*t3.*t4.*t15;
t109 = qd2_3.*t2.*t6.*t7;
t110 = qd2_2.*t7.*t9.*t19;
t111 = qd2_2.*t3.*t6.*t8.*t19;
t112 = t3.*t5.*t8.*t9.*t19.*(3.0./2.0);
t113 = qd2_2.*t3.*t5.*t8.*t9.*t15;
t81 = t71+t72+t73+t74+t75+t76+t77+t78+t79+t80-t109-t110-t111-t112-t113;
t82 = t2.*t8.*t19;
t83 = t4.*t5.*t9.*t19;
t91 = t4.*t6.*t15;
t84 = t82+t83-t91;
t85 = t5.*t8;
t92 = t2.*t4.*t9;
t86 = t85-t92;
t87 = t2.*t8.*t15;
t88 = t4.*t6.*t19;
t89 = t4.*t5.*t9.*t15;
t90 = t87+t88+t89;
t93 = t11.*t12;
t94 = t10.*t13.*t14;
t95 = t93+t94;
t96 = t12.*t13;
t98 = t10.*t11.*t14;
t97 = t96-t98;
t114 = t36.*t97;
t115 = t41.*t95;
t116 = t114+t115-t10.*t28.*t45;
t117 = t66.*t97;
t118 = t10.*t28.*t68;
t119 = t117+t118-t61.*t95;
t120 = t2.*t6.*(3.0./2.0);
t121 = qd1_2.*t4.*t5;
t122 = qd1_2.*t2.*t8.*t9;
t123 = t120+t121+t122+t2.*t8.*t15.*(3.0./4.0)+t4.*t6.*t19.*(3.0./4.0)+t4.*t5.*t9.*t15.*(3.0./4.0);
t124 = t9.*t15.*(3.0./2.0)+t5.*t6.*t19.*(3.0./2.0)-qd1_2.*t2.*t4.*t19-qd1_2.*t6.*t8.*t15+qd1_2.*t5.*t8.*t9.*t19;
t125 = t5.*t8.*(3.0./4.0)-t9.*t19.*(3.0./2.0)-t2.*t4.*t9.*(3.0./4.0)+t5.*t6.*t15.*(3.0./2.0)-qd1_2.*t2.*t4.*t15+qd1_2.*t6.*t8.*t19+qd1_2.*t5.*t8.*t9.*t15;
out1 = reshape([t4.^2.*t6.^2.*(1.0./1.0e1)+(t2.*t8+t4.*t5.*t9).*(t2.*t8.*(1.0./1.0e1)+t4.*t5.*t9.*(1.0./1.0e1))+t86.*(t5.*t8.*(1.0./1.0e1)-t2.*t4.*t9.*(1.0./1.0e1))+(t9.*(3.0./2.0)-qd1_2.*t6.*t8).*(t9.*(3.0./2.0e1)-qd1_2.*t6.*t8.*(1.0./1.0e1))+(t47.*t56+t70.*t81).*(t47.*t56.*2.0+t70.*t81.*2.0)+(t56.*t116-t81.*t119).*(t56.*t116.*2.0-t81.*t119.*2.0)+(t5.*t6.*(3.0./2.0)-qd1_2.*t2.*t4+qd1_2.*t5.*t8.*t9).*(t5.*t6.*(3.0./2.0e1)-qd1_2.*t2.*t4.*(1.0./1.0e1)+qd1_2.*t5.*t8.*t9.*(1.0./1.0e1))+(t14.*t90-t11.*t28.*t86+t13.*t28.*t84).*(t14.*t90.*(1.0./1.0e1)-t11.*t28.*t86.*(1.0./1.0e1)+t13.*t28.*t84.*(1.0./1.0e1))+(t2.*t6.*(3.0./2.0e1)+qd1_2.*t4.*t5.*(1.0./1.0e1)+qd1_2.*t2.*t8.*t9.*(1.0./1.0e1)).*(t120+t121+t122)+(t24.*t84+t27.*t86+t12.*t28.*t90).*(t24.*t84.*(1.0./1.0e1)+t27.*t86.*(1.0./1.0e1)+t12.*t28.*t90.*(1.0./1.0e1))+(t84.*t95+t86.*t97-t10.*t28.*t90).*(t84.*t95.*(1.0./1.0e1)+t86.*t97.*(1.0./1.0e1)-t10.*t28.*t90.*(1.0./1.0e1))+t84.^2+t86.^2+t90.^2+t123.^2+t124.^2+t125.^2+(t5.*t8.*t14.*(3.0./2.0)-t9.*t14.*t19.*(3.0./2.0)-t2.*t4.*t9.*t14.*(3.0./2.0)+t5.*t6.*t14.*t15.*(3.0./2.0)+t2.*t6.*t11.*t28.*(3.0./2.0)+t9.*t13.*t15.*t28.*(3.0./2.0)-qd1_2.*t2.*t4.*t14.*t15+qd1_2.*t6.*t8.*t14.*t19+qd1_2.*t4.*t5.*t11.*t28+qd2_3.*t4.*t6.*t14.*t15-qd2_3.*t2.*t8.*t14.*t19-qd2_2.*t5.*t8.*t13.*t28+t2.*t8.*t11.*t15.*t28.*(3.0./2.0)+t4.*t6.*t11.*t19.*t28.*(3.0./2.0)+t5.*t6.*t13.*t19.*t28.*(3.0./2.0)+qd1_2.*t5.*t8.*t9.*t14.*t15+qd1_2.*t2.*t8.*t9.*t11.*t28-qd2_3.*t4.*t5.*t9.*t14.*t19-qd1_2.*t2.*t4.*t13.*t19.*t28+qd2_2.*t2.*t4.*t9.*t13.*t28-qd1_2.*t6.*t8.*t13.*t15.*t28+qd2_2.*t4.*t6.*t11.*t15.*t28+qd2_3.*t2.*t8.*t13.*t15.*t28-qd2_2.*t2.*t8.*t11.*t19.*t28+qd2_3.*t4.*t6.*t13.*t19.*t28+t4.*t5.*t9.*t11.*t15.*t28.*(3.0./2.0)+qd1_2.*t5.*t8.*t9.*t13.*t19.*t28+qd2_3.*t4.*t5.*t9.*t13.*t15.*t28-qd2_2.*t4.*t5.*t9.*t11.*t19.*t28).*(t5.*t8.*t14.*3.0-t9.*t14.*t19.*3.0-t2.*t4.*t9.*t14.*3.0+t5.*t6.*t14.*t15.*3.0+t2.*t6.*t11.*t28.*3.0+t9.*t13.*t15.*t28.*3.0-qd1_2.*t2.*t4.*t14.*t15.*2.0+qd1_2.*t6.*t8.*t14.*t19.*2.0+qd1_2.*t4.*t5.*t11.*t28.*2.0+qd2_3.*t4.*t6.*t14.*t15.*2.0-qd2_3.*t2.*t8.*t14.*t19.*2.0-qd2_2.*t5.*t8.*t13.*t28.*2.0+t2.*t8.*t11.*t15.*t28.*3.0+t4.*t6.*t11.*t19.*t28.*3.0+t5.*t6.*t13.*t19.*t28.*3.0+qd1_2.*t5.*t8.*t9.*t14.*t15.*2.0+qd1_2.*t2.*t8.*t9.*t11.*t28.*2.0-qd2_3.*t4.*t5.*t9.*t14.*t19.*2.0-qd1_2.*t2.*t4.*t13.*t19.*t28.*2.0+qd2_2.*t2.*t4.*t9.*t13.*t28.*2.0-qd1_2.*t6.*t8.*t13.*t15.*t28.*2.0+qd2_2.*t4.*t6.*t11.*t15.*t28.*2.0+qd2_3.*t2.*t8.*t13.*t15.*t28.*2.0-qd2_2.*t2.*t8.*t11.*t19.*t28.*2.0+qd2_3.*t4.*t6.*t13.*t19.*t28.*2.0+t4.*t5.*t9.*t11.*t15.*t28.*3.0+qd1_2.*t5.*t8.*t9.*t13.*t19.*t28.*2.0+qd2_3.*t4.*t5.*t9.*t13.*t15.*t28.*2.0-qd2_2.*t4.*t5.*t9.*t11.*t19.*t28.*2.0)+2.5e1./1.6e1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[6,6]);
