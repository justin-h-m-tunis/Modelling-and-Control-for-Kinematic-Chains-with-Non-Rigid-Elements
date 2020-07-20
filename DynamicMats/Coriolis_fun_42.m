function out1 = Coriolis_fun_42(in1,in2,in3,in4)
%CORIOLIS_FUN_42
%    OUT1 = CORIOLIS_FUN_42(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    14-May-2020 18:40:32

dqt2 = in3(2,:);
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
qt2 = in1(2,:);
t2 = sin(qd2_6);
t3 = sin(qd1_6);
t4 = sin(qt2);
t5 = cos(qd1_4);
t6 = cos(qd1_6);
t7 = cos(qd2_5);
t8 = cos(qd2_6);
t9 = cos(qt2);
t10 = cos(qd2_4);
t11 = sin(qd2_4);
t12 = sin(qd2_5);
t13 = cos(qd1_5);
t14 = sin(qd1_4);
t15 = sin(qd1_5);
t16 = t3.*t8.*t9.*t10.*2.0;
t17 = t3.*t4.*t11.*t12.*2.0;
t18 = t4.*t5.*t6.*t8.*t10.*2.0;
t19 = t2.*t5.*t6.*t7.*t9.*t10.*2.0;
t200 = t5.*t6.*t9.*t11.*t12.*2.0;
t201 = t2.*t3.*t4.*t7.*t10.*2.0;
t20 = t16+t17+t18+t19-t200-t201;
t21 = t3.*t8.*t9.*t11.*2.0;
t22 = t5.*t6.*t9.*t10.*t12.*2.0;
t23 = t4.*t5.*t6.*t8.*t11.*2.0;
t24 = t2.*t5.*t6.*t7.*t9.*t11.*2.0;
t202 = t3.*t4.*t10.*t12.*2.0;
t203 = t2.*t3.*t4.*t7.*t11.*2.0;
t25 = t21+t22+t23+t24-t202-t203;
t26 = t2.*t3.*t9.*2.0;
t27 = t2.*t4.*t5.*t6.*2.0;
t28 = t3.*t4.*t7.*t8.*2.0;
t199 = t5.*t6.*t7.*t8.*t9.*2.0;
t29 = t26+t27+t28-t199;
t30 = t7.*t11;
t31 = t2.*t10.*t12;
t32 = t30+t31;
t33 = t4.*t8.*t10.*t13.*t14.*2.0;
t34 = t2.*t7.*t9.*t10.*t13.*t14.*2.0;
t35 = t3.*t4.*t5.*t8.*t10.*t15.*2.0;
t36 = t2.*t4.*t6.*t7.*t10.*t15.*2.0;
t37 = t2.*t3.*t5.*t7.*t9.*t10.*t15.*2.0;
t216 = t6.*t8.*t9.*t10.*t15.*2.0;
t217 = t9.*t11.*t12.*t13.*t14.*2.0;
t218 = t4.*t6.*t11.*t12.*t15.*2.0;
t219 = t3.*t5.*t9.*t11.*t12.*t15.*2.0;
t38 = t33+t34+t35+t36+t37-t216-t217-t218-t219;
t39 = t7.*t10;
t75 = t2.*t11.*t12;
t40 = t39-t75;
t41 = t9.*t10.*t12.*t13.*t14.*2.0;
t42 = t4.*t8.*t11.*t13.*t14.*2.0;
t43 = t4.*t6.*t10.*t12.*t15.*2.0;
t44 = t3.*t5.*t9.*t10.*t12.*t15.*2.0;
t45 = t2.*t7.*t9.*t11.*t13.*t14.*2.0;
t46 = t3.*t4.*t5.*t8.*t11.*t15.*2.0;
t47 = t2.*t4.*t6.*t7.*t11.*t15.*2.0;
t48 = t2.*t3.*t5.*t7.*t9.*t11.*t15.*2.0;
t207 = t6.*t8.*t9.*t11.*t15.*2.0;
t49 = t41+t42+t43+t44+t45+t46+t47+t48-t207;
t50 = t2.*t6.*t9.*t15.*2.0;
t51 = t7.*t8.*t9.*t13.*t14.*2.0;
t52 = t4.*t6.*t7.*t8.*t15.*2.0;
t53 = t3.*t5.*t7.*t8.*t9.*t15.*2.0;
t211 = t2.*t4.*t13.*t14.*2.0;
t212 = t2.*t3.*t4.*t5.*t15.*2.0;
t54 = t50+t51+t52+t53-t211-t212;
t55 = qd1_3.*t3.*t8.*t9.*t10.*2.0;
t56 = qd1_3.*t3.*t4.*t11.*t12.*2.0;
t57 = t6.*t7.*t9.*t11.*t13.*3.0;
t58 = t4.*t7.*t11.*t14.*t15.*3.0;
t59 = qd1_2.*t6.*t8.*t9.*t10.*t15.*2.0;
t60 = qd1_3.*t4.*t5.*t6.*t8.*t10.*2.0;
t61 = qd2_3.*t4.*t6.*t8.*t10.*t13.*2.0;
t62 = qd1_2.*t9.*t11.*t12.*t13.*t14.*2.0;
t63 = qd2_2.*t7.*t9.*t11.*t14.*t15.*2.0;
t64 = qd1_2.*t4.*t6.*t11.*t12.*t15.*2.0;
t65 = t2.*t6.*t9.*t10.*t12.*t13.*3.0;
t66 = t2.*t4.*t10.*t12.*t14.*t15.*3.0;
t67 = qd1_2.*t3.*t5.*t9.*t11.*t12.*t15.*2.0;
t68 = qd2_3.*t3.*t4.*t5.*t11.*t12.*t13.*2.0;
t69 = qd2_2.*t2.*t9.*t10.*t12.*t14.*t15.*2.0;
t70 = qd2_3.*t2.*t4.*t7.*t10.*t14.*t15.*2.0;
t71 = qd1_3.*t2.*t5.*t6.*t7.*t9.*t10.*2.0;
t72 = qd2_3.*t3.*t5.*t8.*t9.*t10.*t13.*2.0;
t73 = qd2_3.*t2.*t6.*t7.*t9.*t10.*t13.*2.0;
t230 = qd1_2.*t4.*t8.*t10.*t13.*t14.*2.0;
t231 = qd1_3.*t5.*t6.*t9.*t11.*t12.*2.0;
t232 = qd2_2.*t4.*t6.*t7.*t11.*t13.*2.0;
t233 = qd2_3.*t8.*t9.*t10.*t14.*t15.*2.0;
t234 = qd2_3.*t6.*t9.*t11.*t12.*t13.*2.0;
t235 = qd1_3.*t2.*t3.*t4.*t7.*t10.*2.0;
t236 = qd2_3.*t4.*t11.*t12.*t14.*t15.*2.0;
t237 = t3.*t4.*t5.*t7.*t11.*t13.*3.0;
t238 = t2.*t3.*t4.*t5.*t10.*t12.*t13.*3.0;
t239 = qd1_2.*t2.*t7.*t9.*t10.*t13.*t14.*2.0;
t240 = qd2_2.*t3.*t5.*t7.*t9.*t11.*t13.*2.0;
t241 = qd1_2.*t3.*t4.*t5.*t8.*t10.*t15.*2.0;
t242 = qd1_2.*t2.*t4.*t6.*t7.*t10.*t15.*2.0;
t243 = qd2_2.*t2.*t4.*t6.*t10.*t12.*t13.*2.0;
t244 = qd1_2.*t2.*t3.*t5.*t7.*t9.*t10.*t15.*2.0;
t245 = qd2_2.*t2.*t3.*t5.*t9.*t10.*t12.*t13.*2.0;
t246 = qd2_3.*t2.*t3.*t4.*t5.*t7.*t10.*t13.*2.0;
t74 = t55+t56+t57+t58+t59+t60+t61+t62+t63+t64+t65+t66+t67+t68+t69+t70+t71+t72+t73-t230-t231-t232-t233-t234-t235-t236-t237-t238-t239-t240-t241-t242-t243-t244-t245-t246;
t76 = qd1_3.*t3.*t4.*t10.*t12.*2.0;
t77 = t6.*t7.*t9.*t10.*t13.*3.0;
t78 = t4.*t7.*t10.*t14.*t15.*3.0;
t79 = qd1_2.*t9.*t10.*t12.*t13.*t14.*2.0;
t80 = qd2_2.*t7.*t9.*t10.*t14.*t15.*2.0;
t81 = qd1_2.*t4.*t8.*t11.*t13.*t14.*2.0;
t82 = qd1_2.*t4.*t6.*t10.*t12.*t15.*2.0;
t83 = qd2_3.*t8.*t9.*t11.*t14.*t15.*2.0;
t84 = qd1_3.*t2.*t3.*t4.*t7.*t11.*2.0;
t85 = qd1_2.*t3.*t4.*t5.*t8.*t11.*t15.*2.0;
t86 = qd1_2.*t2.*t4.*t6.*t7.*t11.*t15.*2.0;
t87 = qd2_2.*t2.*t4.*t6.*t11.*t12.*t13.*2.0;
t88 = t2.*t3.*t4.*t5.*t11.*t12.*t13.*3.0;
t89 = qd1_2.*t3.*t5.*t9.*t10.*t12.*t15.*2.0;
t90 = qd1_2.*t2.*t7.*t9.*t11.*t13.*t14.*2.0;
t91 = qd2_3.*t3.*t4.*t5.*t10.*t12.*t13.*2.0;
t92 = qd1_2.*t2.*t3.*t5.*t7.*t9.*t11.*t15.*2.0;
t93 = qd2_2.*t2.*t3.*t5.*t9.*t11.*t12.*t13.*2.0;
t94 = qd2_3.*t2.*t3.*t4.*t5.*t7.*t11.*t13.*2.0;
t248 = qd1_3.*t3.*t8.*t9.*t11.*2.0;
t249 = qd1_3.*t5.*t6.*t9.*t10.*t12.*2.0;
t250 = qd2_2.*t4.*t6.*t7.*t10.*t13.*2.0;
t251 = qd2_3.*t6.*t9.*t10.*t12.*t13.*2.0;
t252 = qd1_2.*t6.*t8.*t9.*t11.*t15.*2.0;
t253 = qd1_3.*t4.*t5.*t6.*t8.*t11.*2.0;
t254 = qd2_3.*t4.*t6.*t8.*t11.*t13.*2.0;
t255 = qd2_3.*t4.*t10.*t12.*t14.*t15.*2.0;
t256 = t3.*t4.*t5.*t7.*t10.*t13.*3.0;
t257 = t2.*t6.*t9.*t11.*t12.*t13.*3.0;
t258 = t2.*t4.*t11.*t12.*t14.*t15.*3.0;
t259 = qd2_2.*t2.*t9.*t11.*t12.*t14.*t15.*2.0;
t260 = qd2_3.*t2.*t4.*t7.*t11.*t14.*t15.*2.0;
t261 = qd2_2.*t3.*t5.*t7.*t9.*t10.*t13.*2.0;
t262 = qd1_3.*t2.*t5.*t6.*t7.*t9.*t11.*2.0;
t263 = qd2_3.*t3.*t5.*t8.*t9.*t11.*t13.*2.0;
t264 = qd2_3.*t2.*t6.*t7.*t9.*t11.*t13.*2.0;
t95 = t76+t77+t78+t79+t80+t81+t82+t83+t84+t85+t86+t87+t88+t89+t90+t91+t92+t93+t94-t248-t249-t250-t251-t252-t253-t254-t255-t256-t257-t258-t259-t260-t261-t262-t263-t264;
t96 = qd1_3.*t2.*t3.*t9.*2.0;
t97 = qd1_2.*t2.*t6.*t9.*t15.*2.0;
t98 = qd1_3.*t2.*t4.*t5.*t6.*2.0;
t99 = qd1_3.*t3.*t4.*t7.*t8.*2.0;
t100 = qd2_3.*t2.*t4.*t6.*t13.*2.0;
t101 = qd1_2.*t7.*t8.*t9.*t13.*t14.*2.0;
t102 = qd1_2.*t4.*t6.*t7.*t8.*t15.*2.0;
t103 = qd2_3.*t2.*t3.*t5.*t9.*t13.*2.0;
t104 = qd2_2.*t4.*t6.*t8.*t12.*t13.*2.0;
t105 = t3.*t4.*t5.*t8.*t12.*t13.*3.0;
t106 = qd1_2.*t3.*t5.*t7.*t8.*t9.*t15.*2.0;
t107 = qd2_2.*t3.*t5.*t8.*t9.*t12.*t13.*2.0;
t108 = qd2_3.*t3.*t4.*t5.*t7.*t8.*t13.*2.0;
t220 = qd1_2.*t2.*t4.*t13.*t14.*2.0;
t221 = qd2_3.*t2.*t9.*t14.*t15.*2.0;
t222 = t6.*t8.*t9.*t12.*t13.*3.0;
t223 = t4.*t8.*t12.*t14.*t15.*3.0;
t224 = qd2_2.*t8.*t9.*t12.*t14.*t15.*2.0;
t225 = qd2_3.*t4.*t7.*t8.*t14.*t15.*2.0;
t226 = qd1_2.*t2.*t3.*t4.*t5.*t15.*2.0;
t227 = qd1_3.*t5.*t6.*t7.*t8.*t9.*2.0;
t228 = qd2_3.*t6.*t7.*t8.*t9.*t13.*2.0;
t109 = t96+t97+t98+t99+t100+t101+t102+t103+t104+t105+t106+t107+t108-t220-t221-t222-t223-t224-t225-t226-t227-t228;
t110 = t3.*t7.*t9.*t11.*3.0;
t111 = qd2_3.*t3.*t4.*t8.*t10.*2.0;
t112 = t4.*t5.*t6.*t7.*t11.*3.0;
t113 = t2.*t3.*t9.*t10.*t12.*3.0;
t114 = qd2_2.*t5.*t6.*t7.*t9.*t11.*2.0;
t115 = qd2_3.*t2.*t3.*t7.*t9.*t10.*2.0;
t116 = qd1_3.*t9.*t11.*t12.*t14.*t15.*2.0;
t117 = t2.*t4.*t5.*t6.*t10.*t12.*3.0;
t118 = qd1_3.*t3.*t4.*t5.*t8.*t10.*t13.*2.0;
t119 = qd1_3.*t2.*t4.*t6.*t7.*t10.*t13.*2.0;
t120 = qd2_2.*t2.*t5.*t6.*t9.*t10.*t12.*2.0;
t121 = qd2_3.*t2.*t4.*t5.*t6.*t7.*t10.*2.0;
t122 = qd1_3.*t2.*t3.*t5.*t7.*t9.*t10.*t13.*2.0;
t279 = qd2_2.*t3.*t4.*t7.*t11.*2.0;
t280 = qd2_3.*t3.*t9.*t11.*t12.*2.0;
t281 = qd1_3.*t4.*t8.*t10.*t14.*t15.*2.0;
t282 = qd1_3.*t4.*t6.*t11.*t12.*t13.*2.0;
t283 = qd2_3.*t4.*t5.*t6.*t11.*t12.*2.0;
t284 = qd2_2.*t2.*t3.*t4.*t10.*t12.*2.0;
t285 = qd1_3.*t6.*t8.*t9.*t10.*t13.*2.0;
t286 = qd2_3.*t5.*t6.*t8.*t9.*t10.*2.0;
t287 = qd1_3.*t3.*t5.*t9.*t11.*t12.*t13.*2.0;
t288 = qd1_3.*t2.*t7.*t9.*t10.*t14.*t15.*2.0;
t123 = t110+t111+t112+t113+t114+t115+t116+t117+t118+t119+t120+t121+t122-t279-t280-t281-t282-t283-t284-t285-t286-t287-t288;
t124 = qd2_2.*t3.*t4.*t7.*t10.*2.0;
t125 = qd2_3.*t3.*t9.*t10.*t12.*2.0;
t126 = qd2_3.*t3.*t4.*t8.*t11.*2.0;
t127 = t2.*t3.*t9.*t11.*t12.*3.0;
t128 = qd1_3.*t4.*t6.*t10.*t12.*t13.*2.0;
t129 = qd2_3.*t4.*t5.*t6.*t10.*t12.*2.0;
t130 = qd2_3.*t2.*t3.*t7.*t9.*t11.*2.0;
t131 = t2.*t4.*t5.*t6.*t11.*t12.*3.0;
t132 = qd1_3.*t3.*t5.*t9.*t10.*t12.*t13.*2.0;
t133 = qd1_3.*t3.*t4.*t5.*t8.*t11.*t13.*2.0;
t134 = qd1_3.*t2.*t4.*t6.*t7.*t11.*t13.*2.0;
t135 = qd2_2.*t2.*t5.*t6.*t9.*t11.*t12.*2.0;
t136 = qd2_3.*t2.*t4.*t5.*t6.*t7.*t11.*2.0;
t137 = qd1_3.*t2.*t3.*t5.*t7.*t9.*t11.*t13.*2.0;
t265 = t3.*t7.*t9.*t10.*3.0;
t266 = t4.*t5.*t6.*t7.*t10.*3.0;
t267 = qd1_3.*t6.*t8.*t9.*t11.*t13.*2.0;
t268 = qd2_3.*t5.*t6.*t8.*t9.*t11.*2.0;
t269 = qd1_3.*t9.*t10.*t12.*t14.*t15.*2.0;
t270 = qd1_3.*t4.*t8.*t11.*t14.*t15.*2.0;
t271 = qd2_2.*t2.*t3.*t4.*t11.*t12.*2.0;
t272 = qd2_2.*t5.*t6.*t7.*t9.*t10.*2.0;
t273 = qd1_3.*t2.*t7.*t9.*t11.*t14.*t15.*2.0;
t138 = t124+t125+t126+t127+t128+t129+t130+t131+t132+t133+t134+t135+t136+t137-t265-t266-t267-t268-t269-t270-t271-t272-t273;
t139 = t3.*t8.*t9.*t12.*3.0;
t140 = qd1_3.*t2.*t6.*t9.*t13.*2.0;
t141 = qd2_3.*t2.*t5.*t6.*t9.*2.0;
t142 = qd2_3.*t3.*t7.*t8.*t9.*2.0;
t143 = qd1_3.*t2.*t4.*t14.*t15.*2.0;
t144 = t4.*t5.*t6.*t8.*t12.*3.0;
t145 = qd1_3.*t4.*t6.*t7.*t8.*t13.*2.0;
t146 = qd2_2.*t5.*t6.*t8.*t9.*t12.*2.0;
t147 = qd2_3.*t4.*t5.*t6.*t7.*t8.*2.0;
t148 = qd1_3.*t3.*t5.*t7.*t8.*t9.*t13.*2.0;
t289 = qd2_3.*t2.*t3.*t4.*2.0;
t290 = qd2_2.*t3.*t4.*t8.*t12.*2.0;
t291 = qd1_3.*t7.*t8.*t9.*t14.*t15.*2.0;
t292 = qd1_3.*t2.*t3.*t4.*t5.*t13.*2.0;
t149 = t139+t140+t141+t142+t143+t144+t145+t146+t147+t148-t289-t290-t291-t292;
t150 = t4.*t7.*t11.*t13.*t14.*3.0;
t151 = qd2_2.*t7.*t9.*t11.*t13.*t14.*2.0;
t152 = qd1_2.*t4.*t8.*t10.*t14.*t15.*2.0;
t153 = qd1_2.*t4.*t6.*t11.*t12.*t13.*2.0;
t154 = qd2_2.*t4.*t6.*t7.*t11.*t15.*2.0;
t155 = qd2_3.*t6.*t9.*t11.*t12.*t15.*2.0;
t156 = t3.*t4.*t5.*t7.*t11.*t15.*3.0;
t157 = t2.*t4.*t10.*t12.*t13.*t14.*3.0;
t158 = qd1_2.*t6.*t8.*t9.*t10.*t13.*2.0;
t159 = qd2_2.*t2.*t4.*t6.*t10.*t12.*t15.*2.0;
t160 = t2.*t3.*t4.*t5.*t10.*t12.*t15.*3.0;
t161 = qd1_2.*t3.*t5.*t9.*t11.*t12.*t13.*2.0;
t162 = qd1_2.*t2.*t7.*t9.*t10.*t14.*t15.*2.0;
t163 = qd2_2.*t3.*t5.*t7.*t9.*t11.*t15.*2.0;
t164 = qd2_2.*t2.*t9.*t10.*t12.*t13.*t14.*2.0;
t165 = qd2_3.*t2.*t4.*t7.*t10.*t13.*t14.*2.0;
t166 = qd2_2.*t2.*t3.*t5.*t9.*t10.*t12.*t15.*2.0;
t167 = qd2_3.*t2.*t3.*t4.*t5.*t7.*t10.*t15.*2.0;
t314 = t6.*t7.*t9.*t11.*t15.*3.0;
t315 = qd2_3.*t8.*t9.*t10.*t13.*t14.*2.0;
t316 = qd2_3.*t4.*t6.*t8.*t10.*t15.*2.0;
t317 = qd1_2.*t9.*t11.*t12.*t14.*t15.*2.0;
t318 = qd2_3.*t4.*t11.*t12.*t13.*t14.*2.0;
t319 = t2.*t6.*t9.*t10.*t12.*t15.*3.0;
t320 = qd2_3.*t3.*t4.*t5.*t11.*t12.*t15.*2.0;
t321 = qd1_2.*t3.*t4.*t5.*t8.*t10.*t13.*2.0;
t322 = qd1_2.*t2.*t4.*t6.*t7.*t10.*t13.*2.0;
t323 = qd2_3.*t3.*t5.*t8.*t9.*t10.*t15.*2.0;
t324 = qd2_3.*t2.*t6.*t7.*t9.*t10.*t15.*2.0;
t325 = qd1_2.*t2.*t3.*t5.*t7.*t9.*t10.*t13.*2.0;
t168 = t150+t151+t152+t153+t154+t155+t156+t157+t158+t159+t160+t161+t162+t163+t164+t165+t166+t167-t314-t315-t316-t317-t318-t319-t320-t321-t322-t323-t324-t325;
t169 = t4.*t7.*t10.*t13.*t14.*3.0;
t170 = qd2_2.*t7.*t9.*t10.*t13.*t14.*2.0;
t171 = qd1_2.*t4.*t6.*t10.*t12.*t13.*2.0;
t172 = qd2_2.*t4.*t6.*t7.*t10.*t15.*2.0;
t173 = qd2_3.*t8.*t9.*t11.*t13.*t14.*2.0;
t174 = qd2_3.*t6.*t9.*t10.*t12.*t15.*2.0;
t175 = qd2_3.*t4.*t6.*t8.*t11.*t15.*2.0;
t176 = t3.*t4.*t5.*t7.*t10.*t15.*3.0;
t177 = t2.*t6.*t9.*t11.*t12.*t15.*3.0;
t178 = qd1_2.*t3.*t5.*t9.*t10.*t12.*t13.*2.0;
t179 = qd2_2.*t3.*t5.*t7.*t9.*t10.*t15.*2.0;
t180 = qd1_2.*t3.*t4.*t5.*t8.*t11.*t13.*2.0;
t181 = qd1_2.*t2.*t4.*t6.*t7.*t11.*t13.*2.0;
t182 = qd2_3.*t3.*t5.*t8.*t9.*t11.*t15.*2.0;
t183 = qd2_3.*t2.*t6.*t7.*t9.*t11.*t15.*2.0;
t184 = qd1_2.*t2.*t3.*t5.*t7.*t9.*t11.*t13.*2.0;
t299 = t6.*t7.*t9.*t10.*t15.*3.0;
t300 = qd1_2.*t6.*t8.*t9.*t11.*t13.*2.0;
t301 = qd1_2.*t9.*t10.*t12.*t14.*t15.*2.0;
t302 = qd2_3.*t4.*t10.*t12.*t13.*t14.*2.0;
t303 = qd1_2.*t4.*t8.*t11.*t14.*t15.*2.0;
t304 = t2.*t4.*t11.*t12.*t13.*t14.*3.0;
t305 = qd1_2.*t2.*t7.*t9.*t11.*t14.*t15.*2.0;
t306 = qd2_3.*t3.*t4.*t5.*t10.*t12.*t15.*2.0;
t307 = qd2_2.*t2.*t9.*t11.*t12.*t13.*t14.*2.0;
t308 = qd2_3.*t2.*t4.*t7.*t11.*t13.*t14.*2.0;
t309 = qd2_2.*t2.*t4.*t6.*t11.*t12.*t15.*2.0;
t310 = t2.*t3.*t4.*t5.*t11.*t12.*t15.*3.0;
t311 = qd2_2.*t2.*t3.*t5.*t9.*t11.*t12.*t15.*2.0;
t312 = qd2_3.*t2.*t3.*t4.*t5.*t7.*t11.*t15.*2.0;
t185 = t169+t170+t171+t172+t173+t174+t175+t176+t177+t178+t179+t180+t181+t182+t183+t184-t299-t300-t301-t302-t303-t304-t305-t306-t307-t308-t309-t310-t311-t312;
t186 = qd2_3.*t2.*t9.*t13.*t14.*2.0;
t187 = qd2_3.*t2.*t4.*t6.*t15.*2.0;
t188 = t4.*t8.*t12.*t13.*t14.*3.0;
t189 = qd1_2.*t7.*t8.*t9.*t14.*t15.*2.0;
t190 = qd2_2.*t8.*t9.*t12.*t13.*t14.*2.0;
t191 = qd2_3.*t4.*t7.*t8.*t13.*t14.*2.0;
t192 = qd1_2.*t2.*t3.*t4.*t5.*t13.*2.0;
t193 = qd2_3.*t2.*t3.*t5.*t9.*t15.*2.0;
t194 = qd2_2.*t4.*t6.*t8.*t12.*t15.*2.0;
t195 = t3.*t4.*t5.*t8.*t12.*t15.*3.0;
t196 = qd2_2.*t3.*t5.*t8.*t9.*t12.*t15.*2.0;
t197 = qd2_3.*t3.*t4.*t5.*t7.*t8.*t15.*2.0;
t293 = qd1_2.*t2.*t6.*t9.*t13.*2.0;
t294 = qd1_2.*t2.*t4.*t14.*t15.*2.0;
t295 = t6.*t8.*t9.*t12.*t15.*3.0;
t296 = qd1_2.*t4.*t6.*t7.*t8.*t13.*2.0;
t297 = qd2_3.*t6.*t7.*t8.*t9.*t15.*2.0;
t298 = qd1_2.*t3.*t5.*t7.*t8.*t9.*t13.*2.0;
t198 = t186+t187+t188+t189+t190+t191+t192+t193+t194+t195+t196+t197-t293-t294-t295-t296-t297-t298;
t204 = qd2_2.*t7.*t10;
t205 = qd2_3.*t8.*t11;
t247 = qd2_2.*t2.*t11.*t12;
t206 = t204+t205-t247;
t208 = qd2_3.*t2;
t209 = qd2_2.*t8.*t12;
t210 = t208+t209;
t213 = qd2_2.*t7.*t11;
t214 = qd2_2.*t2.*t10.*t12;
t229 = qd2_3.*t8.*t10;
t215 = t213+t214-t229;
t274 = t11.*t12;
t313 = t2.*t7.*t10;
t275 = t274-t313;
t276 = t10.*t12;
t277 = t2.*t7.*t11;
t278 = t276+t277;
t326 = t2.*t6.*t9.*t13.*(1.0./1.0e1);
t327 = t2.*t4.*t14.*t15.*(1.0./1.0e1);
t328 = t4.*t6.*t7.*t8.*t13.*(1.0./1.0e1);
t329 = t3.*t5.*t7.*t8.*t9.*t13.*(1.0./1.0e1);
t388 = t7.*t8.*t9.*t14.*t15.*(1.0./1.0e1);
t389 = t2.*t3.*t4.*t5.*t13.*(1.0./1.0e1);
t330 = t326+t327+t328+t329-t388-t389;
t331 = t6.*t8.*t9.*t10.*t13.*(1.0./1.0e1);
t332 = t4.*t8.*t10.*t14.*t15.*(1.0./1.0e1);
t333 = t4.*t6.*t11.*t12.*t13.*(1.0./1.0e1);
t334 = t3.*t5.*t9.*t11.*t12.*t13.*(1.0./1.0e1);
t335 = t2.*t7.*t9.*t10.*t14.*t15.*(1.0./1.0e1);
t379 = t9.*t11.*t12.*t14.*t15.*(1.0./1.0e1);
t380 = t3.*t4.*t5.*t8.*t10.*t13.*(1.0./1.0e1);
t381 = t2.*t4.*t6.*t7.*t10.*t13.*(1.0./1.0e1);
t382 = t2.*t3.*t5.*t7.*t9.*t10.*t13.*(1.0./1.0e1);
t336 = t331+t332+t333+t334+t335-t379-t380-t381-t382;
t337 = t4.*t6.*t10.*t12.*t13.*(1.0./1.0e1);
t338 = t3.*t5.*t9.*t10.*t12.*t13.*(1.0./1.0e1);
t339 = t3.*t4.*t5.*t8.*t11.*t13.*(1.0./1.0e1);
t340 = t2.*t4.*t6.*t7.*t11.*t13.*(1.0./1.0e1);
t341 = t2.*t3.*t5.*t7.*t9.*t11.*t13.*(1.0./1.0e1);
t383 = t6.*t8.*t9.*t11.*t13.*(1.0./1.0e1);
t384 = t9.*t10.*t12.*t14.*t15.*(1.0./1.0e1);
t385 = t4.*t8.*t11.*t14.*t15.*(1.0./1.0e1);
t386 = t2.*t7.*t9.*t11.*t14.*t15.*(1.0./1.0e1);
t342 = t337+t338+t339+t340+t341-t383-t384-t385-t386;
t343 = t2.*t3.*t9.*(1.0./1.0e1);
t344 = t2.*t4.*t5.*t6.*(1.0./1.0e1);
t345 = t3.*t4.*t7.*t8.*(1.0./1.0e1);
t394 = t5.*t6.*t7.*t8.*t9.*(1.0./1.0e1);
t346 = t343+t344+t345-t394;
t347 = t3.*t8.*t9.*t10.*(1.0./1.0e1);
t348 = t3.*t4.*t11.*t12.*(1.0./1.0e1);
t349 = t4.*t5.*t6.*t8.*t10.*(1.0./1.0e1);
t350 = t2.*t5.*t6.*t7.*t9.*t10.*(1.0./1.0e1);
t390 = t5.*t6.*t9.*t11.*t12.*(1.0./1.0e1);
t391 = t2.*t3.*t4.*t7.*t10.*(1.0./1.0e1);
t351 = t347+t348+t349+t350-t390-t391;
t352 = t3.*t8.*t9.*t11.*(1.0./1.0e1);
t353 = t5.*t6.*t9.*t10.*t12.*(1.0./1.0e1);
t354 = t4.*t5.*t6.*t8.*t11.*(1.0./1.0e1);
t355 = t2.*t5.*t6.*t7.*t9.*t11.*(1.0./1.0e1);
t392 = t3.*t4.*t10.*t12.*(1.0./1.0e1);
t393 = t2.*t3.*t4.*t7.*t11.*(1.0./1.0e1);
t356 = t352+t353+t354+t355-t392-t393;
t357 = t2.*t6.*t9.*t15.*(1.0./1.0e1);
t358 = t7.*t8.*t9.*t13.*t14.*(1.0./1.0e1);
t359 = t4.*t6.*t7.*t8.*t15.*(1.0./1.0e1);
t360 = t3.*t5.*t7.*t8.*t9.*t15.*(1.0./1.0e1);
t400 = t2.*t4.*t13.*t14.*(1.0./1.0e1);
t401 = t2.*t3.*t4.*t5.*t15.*(1.0./1.0e1);
t361 = t357+t358+t359+t360-t400-t401;
t362 = t4.*t8.*t10.*t13.*t14.*(1.0./1.0e1);
t363 = t2.*t7.*t9.*t10.*t13.*t14.*(1.0./1.0e1);
t364 = t3.*t4.*t5.*t8.*t10.*t15.*(1.0./1.0e1);
t365 = t2.*t4.*t6.*t7.*t10.*t15.*(1.0./1.0e1);
t366 = t2.*t3.*t5.*t7.*t9.*t10.*t15.*(1.0./1.0e1);
t395 = t6.*t8.*t9.*t10.*t15.*(1.0./1.0e1);
t396 = t9.*t11.*t12.*t13.*t14.*(1.0./1.0e1);
t397 = t4.*t6.*t11.*t12.*t15.*(1.0./1.0e1);
t398 = t3.*t5.*t9.*t11.*t12.*t15.*(1.0./1.0e1);
t367 = t362+t363+t364+t365+t366-t395-t396-t397-t398;
t368 = t9.*t10.*t12.*t13.*t14.*(1.0./1.0e1);
t369 = t4.*t8.*t11.*t13.*t14.*(1.0./1.0e1);
t370 = t4.*t6.*t10.*t12.*t15.*(1.0./1.0e1);
t371 = t3.*t5.*t9.*t10.*t12.*t15.*(1.0./1.0e1);
t372 = t2.*t7.*t9.*t11.*t13.*t14.*(1.0./1.0e1);
t373 = t3.*t4.*t5.*t8.*t11.*t15.*(1.0./1.0e1);
t374 = t2.*t4.*t6.*t7.*t11.*t15.*(1.0./1.0e1);
t375 = t2.*t3.*t5.*t7.*t9.*t11.*t15.*(1.0./1.0e1);
t399 = t6.*t8.*t9.*t11.*t15.*(1.0./1.0e1);
t376 = t368+t369+t370+t371+t372+t373+t374+t375-t399;
t377 = qd2_2.*t11.*t12;
t387 = qd2_2.*t2.*t7.*t10;
t378 = t377-t387;
out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,-dqt2.*(t2.*t29.*(1.0./2.0)+t8.*t10.*t20.*(1.0./2.0)+t8.*t11.*t25.*(1.0./2.0)),dqt2.*(t20.*t32.*(-1.0./2.0)+t25.*t40.*(1.0./2.0)+t8.*t12.*t29.*(1.0./2.0)),dqt2.*(t25.*t206.*(1.0./2.0)-t20.*t215.*(1.0./2.0)+t29.*t210.*(1.0./2.0)),-dqt2.*(qd2_3.*t20.*t275.*(1.0./2.0)-qd2_3.*t25.*t278.*(1.0./2.0)+qd2_3.*t7.*t8.*t29.*(1.0./2.0)),dqt2.*(t20.*t378.*(1.0./2.0)-qd2_2.*t25.*t278.*(1.0./2.0)+qd2_2.*t7.*t8.*t29.*(1.0./2.0)),0.0,-dqt2.*(t2.*t54.*(-1.0./2.0)+t8.*t10.*t38.*(1.0./2.0)+t8.*t11.*t49.*(1.0./2.0)),-dqt2.*(t32.*t38.*(1.0./2.0)-t40.*t49.*(1.0./2.0)+t8.*t12.*t54.*(1.0./2.0)),-dqt2.*(t38.*t215.*(1.0./2.0)-t49.*t206.*(1.0./2.0)+t54.*t210.*(1.0./2.0)),dqt2.*(qd2_3.*t38.*t275.*(-1.0./2.0)+qd2_3.*t49.*t278.*(1.0./2.0)+qd2_3.*t7.*t8.*t54.*(1.0./2.0)),-dqt2.*(t38.*t378.*(-1.0./2.0)+qd2_2.*t49.*t278.*(1.0./2.0)+qd2_2.*t7.*t8.*t54.*(1.0./2.0)),0.0,dqt2.*(t2.*t109.*(1.0./2.0)+t8.*t10.*t74.*(1.0./2.0)-t8.*t11.*t95.*(1.0./2.0)),dqt2.*(t32.*t74.*(1.0./2.0)+t40.*t95.*(1.0./2.0)-t8.*t12.*t109.*(1.0./2.0)),-dqt2.*(t74.*t215.*(-1.0./2.0)-t95.*t206.*(1.0./2.0)+t109.*t210.*(1.0./2.0)+t275.*t336.*(1.0./2.0)+t278.*t342.*(1.0./2.0)+t7.*t8.*t330.*(1.0./2.0)),dqt2.*(t2.*t330.*(-1.0./2.0)+qd2_3.*t74.*t275.*(1.0./2.0)+qd2_3.*t95.*t278.*(1.0./2.0)-t8.*t10.*t336.*(1.0./2.0)+t8.*t11.*t342.*(1.0./2.0)+qd2_3.*t7.*t8.*t109.*(1.0./2.0)),-dqt2.*(t32.*t336.*(1.0./2.0)+t40.*t342.*(1.0./2.0)+t74.*t378.*(1.0./2.0)+qd2_2.*t95.*t278.*(1.0./2.0)-t8.*t12.*t330.*(1.0./2.0)+qd2_2.*t7.*t8.*t109.*(1.0./2.0)),0.0,dqt2.*(t2.*t149.*(-1.0./2.0)+t8.*t10.*t123.*(1.0./2.0)+t8.*t11.*t138.*(1.0./2.0)),dqt2.*(t32.*t123.*(1.0./2.0)-t40.*t138.*(1.0./2.0)+t8.*t12.*t149.*(1.0./2.0)),dqt2.*(t123.*t215.*(1.0./2.0)-t138.*t206.*(1.0./2.0)+t149.*t210.*(1.0./2.0)-t275.*t351.*(1.0./2.0)+t278.*t356.*(1.0./2.0)-t7.*t8.*t346.*(1.0./2.0)),-dqt2.*(t2.*t346.*(1.0./2.0)-qd2_3.*t123.*t275.*(1.0./2.0)+qd2_3.*t138.*t278.*(1.0./2.0)+t8.*t10.*t351.*(1.0./2.0)+t8.*t11.*t356.*(1.0./2.0)+qd2_3.*t7.*t8.*t149.*(1.0./2.0)),dqt2.*(t32.*t351.*(-1.0./2.0)+t40.*t356.*(1.0./2.0)-t123.*t378.*(1.0./2.0)+qd2_2.*t138.*t278.*(1.0./2.0)+t8.*t12.*t346.*(1.0./2.0)+qd2_2.*t7.*t8.*t149.*(1.0./2.0)),0.0,-dqt2.*(t2.*t198.*(1.0./2.0)-t8.*t10.*t168.*(1.0./2.0)+t8.*t11.*t185.*(1.0./2.0)),dqt2.*(t32.*t168.*(1.0./2.0)+t40.*t185.*(1.0./2.0)+t8.*t12.*t198.*(1.0./2.0)),dqt2.*(t168.*t215.*(1.0./2.0)+t185.*t206.*(1.0./2.0)+t198.*t210.*(1.0./2.0)-t275.*t367.*(1.0./2.0)+t278.*t376.*(1.0./2.0)+t7.*t8.*t361.*(1.0./2.0)),dqt2.*(t2.*t361.*(1.0./2.0)+qd2_3.*t168.*t275.*(1.0./2.0)+qd2_3.*t185.*t278.*(1.0./2.0)-t8.*t10.*t367.*(1.0./2.0)-t8.*t11.*t376.*(1.0./2.0)-qd2_3.*t7.*t8.*t198.*(1.0./2.0)),-dqt2.*(t32.*t367.*(1.0./2.0)-t40.*t376.*(1.0./2.0)+t168.*t378.*(1.0./2.0)+qd2_2.*t185.*t278.*(1.0./2.0)+t8.*t12.*t361.*(1.0./2.0)-qd2_2.*t7.*t8.*t198.*(1.0./2.0))],[6,6]);