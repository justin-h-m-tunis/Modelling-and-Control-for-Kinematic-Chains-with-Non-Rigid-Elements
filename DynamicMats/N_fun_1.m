function out1 = N_fun_1(in1,in2)
%N_FUN_1
%    OUT1 = N_FUN_1(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    20-Jul-2020 14:27:19

qd1_2 = in2(3);
qd1_4 = in2(7);
qd1_5 = in2(9);
qd1_6 = in2(11);
qd2_2 = in2(4);
qd2_3 = in2(6);
qd2_5 = in2(10);
qd2_6 = in2(12);
qt1 = in1(1,:);
qt2 = in1(2,:);
t2 = cos(qt1);
t3 = sin(qd1_6);
t4 = cos(qd1_4);
t5 = sin(qt1);
t6 = sin(qd1_5);
t7 = cos(qd1_5);
t8 = sin(qd1_4);
t9 = cos(qd1_6);
t10 = cos(qt2);
t11 = t6.*t8;
t20 = t3.*t4.*t7;
t12 = t11-t20;
t13 = t2.*t12;
t21 = t4.*t5.*t9;
t14 = t13-t21;
t15 = sin(qt2);
t16 = t3.*t5;
t18 = t2.*t7.*t9;
t17 = t16-t18;
t19 = sin(qd2_6);
t22 = t10.*t14;
t23 = t15.*t17;
t24 = t22+t23;
t25 = cos(qd2_6);
t26 = sin(qd2_5);
t27 = t4.*t6;
t28 = t3.*t7.*t8;
t29 = t27+t28;
t30 = t2.*t29;
t31 = t5.*t8.*t9;
t32 = t30+t31;
t33 = cos(qd2_5);
t34 = t10.*t17;
t35 = t34-t14.*t15;
out1 = [t2.*7.3575+t17.*(t10.*(3.0./2.0)-3.0./2.0).*2.943e1-t3.*t5.*1.4715-t10.*t17.*2.20725e1-t14.*t15.*2.20725e1+t19.*t24.*5.886e1-t5.*(qd1_2-t3.*(3.0./2.0)).*3.0411e1+t24.*(qd2_2-t19.*3.0).*(9.81e2./5.0e1)-t2.*(t7.*t9.*(3.0./2.0)-3.0./2.0).*3.0411e1+t35.*(t25.*t33.*3.0-3.0).*(9.81e2./5.0e1)+t32.*(qd2_3+t25.*t26.*3.0).*(9.81e2./5.0e1)+t2.*t7.*t9.*1.4715-t25.*t26.*t32.*5.886e1-t25.*t33.*t35.*5.886e1,0.0,0.0,0.0,0.0,0.0];
