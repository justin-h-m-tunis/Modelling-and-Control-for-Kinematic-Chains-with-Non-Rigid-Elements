function out1 = Mass_fun_21(in1,in2)
%MASS_FUN_21
%    OUT1 = MASS_FUN_21(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    20-Jul-2020 14:23:25

qd1_2 = in2(3);
qd1_3 = in2(5);
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
t2 = cos(qt2);
t3 = sin(qd1_6);
t4 = cos(qd1_4);
t5 = cos(qd1_6);
t6 = sin(qt2);
t7 = cos(qd1_5);
t8 = sin(qd1_4);
t9 = sin(qd1_5);
t10 = sin(qt1);
t11 = cos(qt1);
t12 = sin(qd2_4);
t13 = sin(qd2_5);
t14 = cos(qd2_4);
t15 = cos(qd2_5);
t16 = sin(qd2_6);
t17 = t8.*t9.*t10;
t18 = t4.*t5.*t11;
t43 = t3.*t4.*t7.*t10;
t19 = t17+t18-t43;
t20 = t3.*t11;
t21 = t5.*t7.*t10;
t22 = t20+t21;
t23 = t12.*t13;
t37 = t14.*t15.*t16;
t24 = t23-t37;
t25 = t12.*t15;
t26 = t13.*t14.*t16;
t27 = t25+t26;
t28 = cos(qd2_6);
t29 = t4.*t5.*t10;
t30 = t3.*t4.*t7.*t11;
t58 = t8.*t9.*t11;
t31 = t29+t30-t58;
t32 = t3.*t10;
t60 = t5.*t7.*t11;
t33 = t32-t60;
t34 = t2.*t3;
t35 = t4.*t5.*t6;
t36 = t34+t35;
t38 = t3.*t6;
t85 = t2.*t4.*t5;
t39 = t38-t85;
t40 = t4.*t9.*t10;
t41 = t3.*t7.*t8.*t10;
t107 = t5.*t8.*t11;
t42 = t40+t41-t107;
t44 = t6.*t19;
t45 = t2.*t22;
t46 = t44+t45;
t47 = t2.*t19;
t110 = t6.*t22;
t48 = t47-t110;
t49 = t11.*(3.0./2.0);
t50 = qd2_3.*t4.*t9.*t11;
t51 = qd2_3.*t5.*t8.*t10;
t52 = qd2_2.*t3.*t6.*t10;
t53 = t2.*t5.*t7.*t11.*(3.0./2.0);
t54 = t6.*t8.*t9.*t11.*(3.0./2.0);
t55 = qd2_3.*t3.*t7.*t8.*t11;
t56 = qd2_2.*t2.*t8.*t9.*t11;
t113 = qd1_2.*t10;
t114 = t2.*t3.*t10.*(3.0./2.0);
t115 = t4.*t5.*t6.*t10.*(3.0./2.0);
t116 = qd2_2.*t2.*t4.*t5.*t10;
t117 = qd2_2.*t5.*t6.*t7.*t11;
t118 = t3.*t4.*t6.*t7.*t11.*(3.0./2.0);
t119 = qd2_2.*t2.*t3.*t4.*t7.*t11;
t57 = t49+t50+t51+t52+t53+t54+t55+t56-t113-t114-t115-t116-t117-t118-t119;
t59 = t6.*t31;
t61 = t2.*t33;
t62 = t59+t61;
t63 = t13.*t14;
t64 = t12.*t15.*t16;
t65 = t63+t64;
t66 = t14.*t15;
t86 = t12.*t13.*t16;
t67 = t66-t86;
t68 = t4.*t9.*t11;
t69 = t5.*t8.*t10;
t70 = t3.*t7.*t8.*t11;
t71 = t68+t69+t70;
t72 = t2.*t31;
t122 = t6.*t33;
t73 = t72-t122;
t74 = t10.*(3.0./2.0);
t75 = qd1_2.*t11;
t76 = t2.*t3.*t11.*(3.0./2.0);
t77 = qd2_3.*t4.*t9.*t10;
t78 = t4.*t5.*t6.*t11.*(3.0./2.0);
t79 = t2.*t5.*t7.*t10.*(3.0./2.0);
t80 = t6.*t8.*t9.*t10.*(3.0./2.0);
t81 = qd2_2.*t2.*t4.*t5.*t11;
t82 = qd2_3.*t3.*t7.*t8.*t10;
t83 = qd2_2.*t2.*t8.*t9.*t10;
t125 = qd2_3.*t5.*t8.*t11;
t126 = qd2_2.*t3.*t6.*t11;
t127 = qd2_2.*t5.*t6.*t7.*t10;
t128 = t3.*t4.*t6.*t7.*t10.*(3.0./2.0);
t129 = qd2_2.*t2.*t3.*t4.*t7.*t10;
t84 = t74+t75+t76+t77+t78+t79+t80+t81+t82+t83-t125-t126-t127-t128-t129;
t87 = t5.*t8.*(3.0./2.0);
t88 = qd1_2.*t4.*t9;
t89 = t2.*t7.*t8.*(3.0./4.0);
t90 = t5.*t6.*t9.*(3.0./4.0);
t91 = qd1_2.*t3.*t7.*t8;
t92 = t2.*t3.*t4.*t9.*(3.0./4.0);
t93 = t87+t88+t89+t90+t91+t92;
t94 = t2.*t3.*(3.0./2.0);
t95 = t4.*t5.*t6.*(3.0./2.0);
t96 = qd1_2.*t3.*t4.*t6.*t7;
t210 = qd1_2.*t2.*t5.*t7;
t211 = qd1_2.*t6.*t8.*t9;
t97 = t94+t95+t96-t210-t211;
t98 = t4.*t5.*(3.0./2.0e1);
t99 = qd1_2.*t3.*t4.*t7.*(1.0./1.0e1);
t238 = qd1_2.*t8.*t9.*(1.0./1.0e1);
t100 = t98+t99-t238;
t101 = t4.*t7;
t135 = t3.*t8.*t9;
t102 = t101-t135;
t103 = t5.*t8.*(3.0./2.0e1);
t104 = qd1_2.*t4.*t9.*(1.0./1.0e1);
t105 = qd1_2.*t3.*t7.*t8.*(1.0./1.0e1);
t106 = t103+t104+t105;
t108 = t27.*t42;
t109 = t24.*t46;
t111 = t14.*t28.*t48;
t112 = t108+t109+t111;
t120 = t57.*t112.*2.0;
t121 = t24.*t62;
t123 = t14.*t28.*t73;
t193 = t27.*t71;
t124 = t121+t123-t193;
t130 = t84.*t124.*2.0;
t131 = t120+t130;
t132 = t6.*t7.*t8;
t133 = t3.*t4.*t6.*t9;
t144 = t2.*t5.*t9;
t134 = t132+t133-t144;
t136 = t42.*t67;
t137 = t46.*t65;
t224 = t12.*t28.*t48;
t138 = t136+t137-t224;
t139 = t57.*t138.*2.0;
t140 = t67.*t71;
t141 = t12.*t28.*t73;
t225 = t62.*t65;
t142 = t140+t141-t225;
t226 = t84.*t142.*2.0;
t143 = t139-t226;
t145 = t2.*t7.*t8;
t146 = t5.*t6.*t9;
t147 = t2.*t3.*t4.*t9;
t148 = t145+t146+t147;
t149 = t4.*t7.*t16.*3.0;
t150 = t2.*t4.*t5.*t16.*3.0;
t151 = t2.*t3.*t15.*t28.*3.0;
t152 = t5.*t8.*t13.*t28.*3.0;
t153 = qd1_2.*t4.*t9.*t13.*t28.*2.0;
t154 = qd1_2.*t5.*t6.*t7.*t16.*2.0;
t155 = qd2_3.*t2.*t5.*t9.*t16.*2.0;
t156 = t4.*t5.*t6.*t15.*t28.*3.0;
t157 = t2.*t7.*t8.*t13.*t28.*3.0;
t158 = t5.*t6.*t9.*t13.*t28.*3.0;
t159 = qd2_3.*t2.*t7.*t8.*t15.*t28.*2.0;
t160 = qd1_2.*t2.*t3.*t4.*t7.*t16.*2.0;
t161 = qd2_2.*t2.*t5.*t9.*t13.*t28.*2.0;
t162 = qd2_3.*t5.*t6.*t9.*t15.*t28.*2.0;
t163 = qd1_2.*t3.*t7.*t8.*t13.*t28.*2.0;
t164 = qd2_2.*t3.*t8.*t9.*t15.*t28.*2.0;
t165 = t2.*t3.*t4.*t9.*t13.*t28.*3.0;
t166 = qd1_2.*t3.*t4.*t6.*t7.*t15.*t28.*2.0;
t167 = qd2_3.*t2.*t3.*t4.*t9.*t15.*t28.*2.0;
t183 = t3.*t6.*t16.*3.0;
t184 = t3.*t8.*t9.*t16.*3.0;
t185 = qd2_2.*t4.*t7.*t15.*t28.*2.0;
t186 = qd1_2.*t2.*t8.*t9.*t16.*2.0;
t187 = qd2_3.*t6.*t7.*t8.*t16.*2.0;
t188 = qd1_2.*t6.*t8.*t9.*t15.*t28.*2.0;
t189 = qd2_2.*t6.*t7.*t8.*t13.*t28.*2.0;
t190 = qd2_3.*t3.*t4.*t6.*t9.*t16.*2.0;
t191 = qd1_2.*t2.*t5.*t7.*t15.*t28.*2.0;
t192 = qd2_2.*t3.*t4.*t6.*t9.*t13.*t28.*2.0;
t168 = t149+t150+t151+t152+t153+t154+t155+t156+t157+t158+t159+t160+t161+t162+t163+t164+t165+t166+t167-t183-t184-t185-t186-t187-t188-t189-t190-t191-t192;
t169 = t4.*t7.*(3.0./4.0);
t170 = t2.*t4.*t5.*(3.0./2.0);
t171 = qd1_2.*t5.*t6.*t7;
t172 = qd1_2.*t2.*t3.*t4.*t7;
t221 = t3.*t6.*(3.0./2.0);
t222 = t3.*t8.*t9.*(3.0./4.0);
t223 = qd1_2.*t2.*t8.*t9;
t173 = t169+t170+t171+t172-t221-t222-t223;
t174 = t3.*(3.0./2.0e1);
t176 = qd1_2.*t5.*t7.*(1.0./1.0e1);
t175 = t174-t176;
t177 = t16.*t148;
t178 = t15.*t28.*t134;
t287 = t13.*t28.*t102;
t179 = t177+t178-t287;
t180 = t15.*t28.*t36;
t181 = t5.*t8.*t13.*t28;
t242 = t16.*t39;
t182 = t180+t181-t242;
t194 = t24.*t134;
t195 = t27.*t102;
t196 = t14.*t28.*t148;
t197 = t194+t195+t196;
t198 = qd2_2.*t2.*t4.*t5;
t230 = qd2_3.*t5.*t8;
t231 = qd2_2.*t3.*t6;
t199 = qd1_2+t94+t95+t198-t230-t231;
t200 = t14.*t28.*t39;
t201 = t5.*t8.*t27;
t249 = t24.*t36;
t202 = t200+t201-t249;
t203 = qd2_3.*t4.*t7;
t204 = t6.*t7.*t8.*(3.0./2.0);
t205 = qd2_2.*t2.*t7.*t8;
t206 = qd2_2.*t5.*t6.*t9;
t207 = t3.*t4.*t6.*t9.*(3.0./2.0);
t208 = qd2_2.*t2.*t3.*t4.*t9;
t235 = t2.*t5.*t9.*(3.0./2.0);
t236 = qd2_3.*t3.*t8.*t9;
t209 = qd1_3+t203+t204+t205+t206+t207+t208-t235-t236;
t212 = t4.*t9;
t213 = t3.*t7.*t8;
t214 = t212+t213;
t215 = t5.*t6.*t7;
t216 = t2.*t3.*t4.*t7;
t237 = t2.*t8.*t9;
t217 = t215+t216-t237;
t218 = t6.*t8.*t9;
t219 = t2.*t5.*t7;
t239 = t3.*t4.*t6.*t7;
t220 = t218+t219-t239;
t227 = t65.*t134;
t228 = t67.*t102;
t254 = t12.*t28.*t148;
t229 = t227+t228-t254;
t232 = t36.*t65;
t233 = t12.*t28.*t39;
t265 = t5.*t8.*t67;
t234 = t232+t233-t265;
t240 = qd1_2.*t4.*t7;
t241 = qd1_3.*t5.*t8;
t243 = t16.*t148.*(1.0./1.0e1);
t244 = t15.*t28.*t134.*(1.0./1.0e1);
t288 = t13.*t28.*t102.*(1.0./1.0e1);
t245 = t243+t244-t288;
t246 = t27.*t214;
t247 = t24.*t220;
t289 = t14.*t28.*t217;
t248 = t246+t247-t289;
t250 = t24.*t134.*(1.0./1.0e1);
t251 = t27.*t102.*(1.0./1.0e1);
t252 = t14.*t28.*t148.*(1.0./1.0e1);
t253 = t250+t251+t252;
t255 = t6.*t8.*t9.*(3.0./2.0);
t256 = qd2_3.*t4.*t9;
t257 = t2.*t5.*t7.*(3.0./2.0);
t258 = qd2_3.*t3.*t7.*t8;
t259 = qd2_2.*t2.*t8.*t9;
t276 = qd2_2.*t5.*t6.*t7;
t277 = t3.*t4.*t6.*t7.*(3.0./2.0);
t278 = qd2_2.*t2.*t3.*t4.*t7;
t260 = t255+t256+t257+t258+t259-t276-t277-t278+3.0./2.0;
t261 = t67.*t214;
t262 = t65.*t220;
t263 = t12.*t28.*t217;
t264 = t261+t262+t263;
t266 = t65.*t134.*(1.0./1.0e1);
t267 = t67.*t102.*(1.0./1.0e1);
t279 = t12.*t28.*t148.*(1.0./1.0e1);
t268 = t266+t267-t279;
t269 = t7.*t8.*(1.0./1.0e1);
t270 = t3.*t4.*t9.*(1.0./1.0e1);
t271 = t269+t270;
t272 = t4.*t7.*(1.0./1.0e1);
t283 = t3.*t8.*t9.*(1.0./1.0e1);
t273 = t272-t283;
t274 = t8.*t9;
t275 = t274-t3.*t4.*t7;
t280 = t7.*t8;
t281 = t3.*t4.*t9;
t282 = t280+t281;
t284 = t16.*t217;
t285 = t13.*t28.*t214;
t286 = t284+t285-t15.*t28.*t220;
t290 = t5.^2;
out1 = reshape([0.0,t36.*t97+t3.*t175-t39.*t173-t131.*t202+t168.*t182+t143.*t234+t5.*t8.*t93+t4.*t5.*t100+t5.*t8.*t106,-t93.*t102-t102.*t106+t97.*t134+t148.*t173+t131.*t197+t168.*t179+t143.*t229+t100.*t282-t5.*t9.*t175,-t97.*(qd1_3.*t2.*t3+qd1_2.*t2.*t5.*t9+qd1_3.*t4.*t5.*t6-qd1_2.*t6.*t7.*t8-qd1_2.*t3.*t4.*t6.*t9)+t173.*(t4.*t9.*(3.0./4.0)+qd1_3.*t3.*t6+t3.*t7.*t8.*(3.0./4.0)-qd1_3.*t2.*t4.*t5+qd1_2.*t2.*t7.*t8+qd1_2.*t5.*t6.*t9+qd1_2.*t2.*t3.*t4.*t9)-t175.*(qd1_3.*t3+qd1_2.*t5.*t9)+t102.*t214+t134.*t220-t148.*t217+t214.*t273+t248.*t253-t245.*t286+t264.*t268+t271.*t275-t106.*(t240+t241-qd1_2.*t3.*t8.*t9)+t131.*(t197.*t199+t202.*t209)+t168.*(t179.*t199-t182.*t209)+t143.*(t199.*t229-t209.*t234)+t100.*(-qd1_3.*t4.*t5+qd1_2.*t7.*t8+qd1_2.*t3.*t4.*t9)-t93.*(t240+t241+t5.*t6.*t7.*(3.0./4.0)-t2.*t8.*t9.*(3.0./4.0)-qd1_2.*t3.*t8.*t9+t2.*t3.*t4.*t7.*(3.0./4.0))-t7.*t9.*t290.*(1.0./1.0e1),t131.*(t27.*t102.*(3.0./2.0)+t24.*t134.*(3.0./2.0)-t197.*t260+t209.*t248+t14.*t28.*t148.*(3.0./2.0))+t143.*(t67.*t102.*(3.0./2.0)+t65.*t134.*(3.0./2.0)+t209.*t264-t229.*t260-t12.*t28.*t148.*(3.0./2.0))-t93.*(t3.*t6.*(3.0./4.0)+qd1_3.*t4.*t9-t2.*t4.*t5.*(3.0./4.0)+qd1_3.*t3.*t7.*t8)+t36.*t134-t39.*t148+t182.*t245-t202.*t253+t234.*t268-t168.*(t5.*t8.*t16.*(3.0./2.0)+qd2_3.*t2.*t3.*t16+t3.*t6.*t13.*t28.*(3.0./2.0)+qd1_3.*t5.*t6.*t7.*t16-qd1_3.*t2.*t8.*t9.*t16+qd2_3.*t4.*t5.*t6.*t16+qd1_3.*t4.*t9.*t13.*t28+qd2_2.*t2.*t3.*t13.*t28+qd2_3.*t3.*t6.*t15.*t28-qd2_2.*t5.*t8.*t15.*t28-t2.*t4.*t5.*t13.*t28.*(3.0./2.0)+qd1_3.*t2.*t3.*t4.*t7.*t16-qd1_3.*t2.*t5.*t7.*t15.*t28+qd1_3.*t3.*t7.*t8.*t13.*t28-qd2_3.*t2.*t4.*t5.*t15.*t28+qd2_2.*t4.*t5.*t6.*t13.*t28-qd1_3.*t6.*t8.*t9.*t15.*t28+qd1_3.*t3.*t4.*t6.*t7.*t15.*t28)-t173.*(t5.*t8.*(3.0./4.0)+qd1_3.*t5.*t6.*t7-qd1_3.*t2.*t8.*t9+qd1_3.*t2.*t3.*t4.*t7)+qd1_3.*t97.*t220-qd1_3.*t106.*t214+qd1_3.*t100.*t275-t3.*t5.*t9.*(1.0./1.0e1)-t5.*t8.*t102+t4.*t5.*t271-t5.*t8.*t273+qd1_3.*t5.*t7.*t175,-t100.*(qd1_2.*t8.*t9-qd1_2.*t3.*t4.*t7)+t173.*(t169+t171+t172-t222-t223)+t102.*t273+t179.*t245+t197.*t253+t229.*t268+t271.*t282+t93.*(t88+t89+t90+t91+t92)-t131.*(t24.*t36.*(3.0./2.0)+t199.*t248+t202.*t260-t5.*t8.*t27.*(3.0./2.0)-t14.*t28.*t39.*(3.0./2.0))-t143.*(t36.*t65.*(3.0./2.0)+t199.*t264-t234.*t260+t12.*t28.*t39.*(3.0./2.0)-t5.*t8.*t67.*(3.0./2.0))+t9.^2.*t290.*(1.0./1.0e1)+t102.^2+t134.^2+t148.^2+t168.*(t16.*t39.*(3.0./2.0)+t182.*t260+t199.*t286-t15.*t28.*t36.*(3.0./2.0)-t5.*t8.*t13.*t28.*(3.0./2.0))-qd1_2.*t97.*t220+qd1_2.*t106.*t214-qd1_2.*t5.*t7.*t175,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[6,6]);
