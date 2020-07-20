function out1 = Coriolis_fun_21(in1,in2,in3,in4)
%CORIOLIS_FUN_21
%    OUT1 = CORIOLIS_FUN_21(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    14-May-2020 18:17:59

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
t3 = cos(qt2);
t4 = sin(qd1_6);
t5 = cos(qd1_6);
t6 = cos(qd2_6);
t7 = cos(qd1_4);
t8 = cos(qd2_5);
t9 = sin(qt2);
t10 = cos(qd1_5);
t11 = sin(qd1_4);
t12 = sin(qd1_5);
t13 = sin(qd2_5);
t14 = sin(qd2_4);
t15 = cos(qd2_4);
t16 = t3.*t4;
t17 = t5.*t7.*t9;
t18 = t16+t17;
t19 = t3.*t5.*t7.*(3.0./2.0);
t20 = qd1_2.*t5.*t9.*t10;
t21 = qd1_2.*t3.*t4.*t7.*t10;
t22 = t9.*t10.*t11.*(3.0./4.0);
t23 = t4.*t7.*t9.*t12.*(3.0./4.0);
t169 = t3.*t5.*t12.*(3.0./4.0);
t24 = t22+t23-t169;
t25 = t7.*t10.*(3.0./4.0);
t167 = t4.*t9.*(3.0./2.0);
t168 = qd1_2.*t3.*t11.*t12;
t185 = t4.*t11.*t12.*(3.0./4.0);
t26 = t19+t20+t21+t25-t167-t168-t185;
t27 = t5.*t8.*t11.*t14.*3.0;
t28 = t4.*t6.*t9.*t15.*3.0;
t29 = qd1_2.*t7.*t8.*t12.*t14.*2.0;
t30 = qd2_2.*t7.*t10.*t13.*t14.*2.0;
t31 = t3.*t8.*t10.*t11.*t14.*3.0;
t32 = t2.*t3.*t4.*t8.*t15.*3.0;
t33 = t4.*t6.*t11.*t12.*t15.*3.0;
t34 = t2.*t5.*t11.*t13.*t15.*3.0;
t35 = t5.*t8.*t9.*t12.*t14.*3.0;
t36 = qd1_2.*t3.*t6.*t11.*t12.*t15.*2.0;
t37 = qd1_2.*t3.*t5.*t10.*t13.*t14.*2.0;
t38 = qd2_2.*t3.*t5.*t8.*t12.*t14.*2.0;
t39 = qd2_3.*t6.*t9.*t10.*t11.*t15.*2.0;
t40 = qd1_2.*t4.*t8.*t10.*t11.*t14.*2.0;
t41 = qd1_2.*t2.*t7.*t12.*t13.*t15.*2.0;
t42 = qd1_2.*t9.*t11.*t12.*t13.*t14.*2.0;
t43 = t2.*t5.*t7.*t8.*t9.*t15.*3.0;
t44 = t3.*t4.*t7.*t8.*t12.*t14.*3.0;
t45 = t2.*t3.*t10.*t11.*t13.*t15.*3.0;
t46 = t2.*t5.*t9.*t12.*t13.*t15.*3.0;
t47 = qd1_2.*t2.*t4.*t10.*t11.*t13.*t15.*2.0;
t48 = qd2_2.*t2.*t4.*t8.*t11.*t12.*t15.*2.0;
t49 = t2.*t3.*t4.*t7.*t12.*t13.*t15.*3.0;
t50 = qd2_3.*t2.*t3.*t8.*t10.*t11.*t15.*2.0;
t51 = qd2_3.*t4.*t6.*t7.*t9.*t12.*t15.*2.0;
t52 = qd2_2.*t2.*t3.*t5.*t12.*t13.*t15.*2.0;
t53 = qd2_3.*t2.*t5.*t8.*t9.*t12.*t15.*2.0;
t54 = qd1_2.*t2.*t4.*t7.*t8.*t9.*t10.*t15.*2.0;
t55 = qd2_3.*t2.*t3.*t4.*t7.*t8.*t12.*t15.*2.0;
t210 = t6.*t7.*t10.*t15.*3.0;
t211 = t3.*t4.*t13.*t14.*3.0;
t212 = t3.*t5.*t6.*t7.*t15.*3.0;
t213 = t5.*t7.*t9.*t13.*t14.*3.0;
t214 = qd2_2.*t2.*t7.*t8.*t10.*t15.*2.0;
t215 = qd1_2.*t5.*t6.*t9.*t10.*t15.*2.0;
t216 = qd2_3.*t3.*t5.*t6.*t12.*t15.*2.0;
t217 = qd2_2.*t8.*t9.*t10.*t11.*t14.*2.0;
t218 = qd2_3.*t3.*t10.*t11.*t13.*t14.*2.0;
t219 = qd2_3.*t5.*t9.*t12.*t13.*t14.*2.0;
t220 = qd2_2.*t4.*t11.*t12.*t13.*t14.*2.0;
t221 = qd1_2.*t4.*t7.*t9.*t10.*t13.*t14.*2.0;
t222 = qd1_2.*t2.*t8.*t9.*t11.*t12.*t15.*2.0;
t223 = qd2_2.*t4.*t7.*t8.*t9.*t12.*t14.*2.0;
t224 = qd2_3.*t3.*t4.*t7.*t12.*t13.*t14.*2.0;
t225 = qd2_2.*t2.*t9.*t10.*t11.*t13.*t15.*2.0;
t226 = qd1_2.*t3.*t4.*t6.*t7.*t10.*t15.*2.0;
t227 = qd1_2.*t2.*t3.*t5.*t8.*t10.*t15.*2.0;
t228 = qd2_2.*t2.*t4.*t7.*t9.*t12.*t13.*t15.*2.0;
t56 = t27+t28+t29+t30+t31+t32+t33+t34+t35+t36+t37+t38+t39+t40+t41+t42+t43+t44+t45+t46+t47+t48+t49+t50+t51+t52+t53+t54+t55-t210-t211-t212-t213-t214-t215-t216-t217-t218-t219-t220-t221-t222-t223-t224-t225-t226-t227-t228;
t57 = t3.*t4.*t6.*t14.*3.0;
t58 = t3.*t5.*t7.*t13.*t15.*3.0;
t59 = t8.*t9.*t10.*t11.*t15.*3.0;
t60 = t5.*t6.*t7.*t9.*t14.*3.0;
t61 = qd2_2.*t3.*t8.*t10.*t11.*t15.*2.0;
t62 = qd1_2.*t5.*t9.*t10.*t13.*t15.*2.0;
t63 = qd2_2.*t5.*t8.*t9.*t12.*t15.*2.0;
t64 = qd2_3.*t3.*t6.*t10.*t11.*t14.*2.0;
t65 = qd2_3.*t3.*t5.*t12.*t13.*t15.*2.0;
t66 = qd2_3.*t5.*t6.*t9.*t12.*t14.*2.0;
t67 = t2.*t3.*t5.*t7.*t8.*t14.*3.0;
t68 = t4.*t7.*t8.*t9.*t12.*t15.*3.0;
t69 = t2.*t3.*t5.*t12.*t13.*t14.*3.0;
t70 = qd1_2.*t3.*t4.*t7.*t10.*t13.*t15.*2.0;
t71 = qd2_2.*t3.*t4.*t7.*t8.*t12.*t15.*2.0;
t72 = qd1_2.*t4.*t6.*t7.*t9.*t10.*t14.*2.0;
t73 = qd1_2.*t2.*t5.*t8.*t9.*t10.*t14.*2.0;
t74 = qd2_3.*t3.*t4.*t6.*t7.*t12.*t14.*2.0;
t75 = qd2_3.*t2.*t3.*t5.*t8.*t12.*t14.*2.0;
t76 = qd1_2.*t2.*t3.*t4.*t7.*t8.*t10.*t14.*2.0;
t229 = t4.*t9.*t13.*t15.*3.0;
t230 = t3.*t5.*t8.*t12.*t15.*3.0;
t231 = t2.*t4.*t8.*t9.*t14.*3.0;
t232 = qd1_2.*t3.*t5.*t6.*t10.*t14.*2.0;
t233 = qd1_2.*t3.*t11.*t12.*t13.*t15.*2.0;
t234 = qd2_3.*t9.*t10.*t11.*t13.*t15.*2.0;
t235 = qd1_2.*t6.*t9.*t11.*t12.*t14.*2.0;
t236 = t2.*t9.*t10.*t11.*t13.*t14.*3.0;
t237 = qd1_2.*t2.*t3.*t8.*t11.*t12.*t14.*2.0;
t238 = qd2_3.*t4.*t7.*t9.*t12.*t13.*t15.*2.0;
t239 = qd2_2.*t2.*t3.*t10.*t11.*t13.*t14.*2.0;
t240 = qd2_3.*t2.*t8.*t9.*t10.*t11.*t14.*2.0;
t241 = qd2_2.*t2.*t5.*t9.*t12.*t13.*t14.*2.0;
t242 = t2.*t4.*t7.*t9.*t12.*t13.*t14.*3.0;
t243 = qd2_2.*t2.*t3.*t4.*t7.*t12.*t13.*t14.*2.0;
t244 = qd2_3.*t2.*t4.*t7.*t8.*t9.*t12.*t14.*2.0;
t77 = t57+t58+t59+t60+t61+t62+t63+t64+t65+t66+t67+t68+t69+t70+t71+t72+t73+t74+t75+t76-t229-t230-t231-t232-t233-t234-t235-t236-t237-t238-t239-t240-t241-t242-t243-t244;
t78 = t3.*t5.*t7.*t13.*t14.*3.0;
t79 = t8.*t9.*t10.*t11.*t14.*3.0;
t80 = t2.*t4.*t8.*t9.*t15.*3.0;
t81 = qd2_2.*t3.*t8.*t10.*t11.*t14.*2.0;
t82 = qd1_2.*t6.*t9.*t11.*t12.*t15.*2.0;
t83 = qd1_2.*t5.*t9.*t10.*t13.*t14.*2.0;
t84 = qd2_2.*t5.*t8.*t9.*t12.*t14.*2.0;
t85 = qd2_3.*t3.*t5.*t12.*t13.*t14.*2.0;
t86 = t4.*t7.*t8.*t9.*t12.*t14.*3.0;
t87 = t2.*t9.*t10.*t11.*t13.*t15.*3.0;
t88 = qd1_2.*t3.*t5.*t6.*t10.*t15.*2.0;
t89 = qd2_2.*t2.*t5.*t9.*t12.*t13.*t15.*2.0;
t90 = t2.*t4.*t7.*t9.*t12.*t13.*t15.*3.0;
t91 = qd1_2.*t3.*t4.*t7.*t10.*t13.*t14.*2.0;
t92 = qd1_2.*t2.*t3.*t8.*t11.*t12.*t15.*2.0;
t93 = qd2_2.*t3.*t4.*t7.*t8.*t12.*t14.*2.0;
t94 = qd2_2.*t2.*t3.*t10.*t11.*t13.*t15.*2.0;
t95 = qd2_3.*t2.*t8.*t9.*t10.*t11.*t15.*2.0;
t96 = qd2_2.*t2.*t3.*t4.*t7.*t12.*t13.*t15.*2.0;
t97 = qd2_3.*t2.*t4.*t7.*t8.*t9.*t12.*t15.*2.0;
t194 = t3.*t4.*t6.*t15.*3.0;
t195 = t4.*t9.*t13.*t14.*3.0;
t196 = t5.*t6.*t7.*t9.*t15.*3.0;
t197 = t3.*t5.*t8.*t12.*t14.*3.0;
t198 = qd2_3.*t3.*t6.*t10.*t11.*t15.*2.0;
t199 = qd2_3.*t5.*t6.*t9.*t12.*t15.*2.0;
t200 = qd1_2.*t3.*t11.*t12.*t13.*t14.*2.0;
t201 = qd2_3.*t9.*t10.*t11.*t13.*t14.*2.0;
t202 = t2.*t3.*t5.*t7.*t8.*t15.*3.0;
t203 = t2.*t3.*t5.*t12.*t13.*t15.*3.0;
t204 = qd2_3.*t4.*t7.*t9.*t12.*t13.*t14.*2.0;
t205 = qd1_2.*t4.*t6.*t7.*t9.*t10.*t15.*2.0;
t206 = qd1_2.*t2.*t5.*t8.*t9.*t10.*t15.*2.0;
t207 = qd2_3.*t3.*t4.*t6.*t7.*t12.*t15.*2.0;
t208 = qd2_3.*t2.*t3.*t5.*t8.*t12.*t15.*2.0;
t209 = qd1_2.*t2.*t3.*t4.*t7.*t8.*t10.*t15.*2.0;
t98 = t78+t79+t80+t81+t82+t83+t84+t85+t86+t87+t88+t89+t90+t91+t92+t93+t94+t95+t96+t97-t194-t195-t196-t197-t198-t199-t200-t201-t202-t203-t204-t205-t206-t207-t208-t209;
t99 = t2.*t7.*t10.*3.0;
t100 = t2.*t3.*t5.*t7.*3.0;
t101 = t3.*t4.*t6.*t8.*3.0;
t102 = t5.*t6.*t11.*t13.*3.0;
t103 = qd1_2.*t6.*t7.*t12.*t13.*2.0;
t104 = qd1_2.*t2.*t5.*t9.*t10.*2.0;
t105 = qd2_3.*t2.*t3.*t5.*t12.*2.0;
t106 = t5.*t6.*t7.*t8.*t9.*3.0;
t107 = t3.*t6.*t10.*t11.*t13.*3.0;
t108 = t5.*t6.*t9.*t12.*t13.*3.0;
t109 = qd2_3.*t3.*t6.*t8.*t10.*t11.*2.0;
t110 = qd1_2.*t2.*t3.*t4.*t7.*t10.*2.0;
t111 = qd2_2.*t3.*t5.*t6.*t12.*t13.*2.0;
t112 = qd2_3.*t5.*t6.*t8.*t9.*t12.*2.0;
t113 = qd1_2.*t4.*t6.*t10.*t11.*t13.*2.0;
t114 = qd2_2.*t4.*t6.*t8.*t11.*t12.*2.0;
t115 = t3.*t4.*t6.*t7.*t12.*t13.*3.0;
t116 = qd1_2.*t4.*t6.*t7.*t8.*t9.*t10.*2.0;
t117 = qd2_3.*t3.*t4.*t6.*t7.*t8.*t12.*2.0;
t171 = t2.*t4.*t9.*3.0;
t172 = t2.*t4.*t11.*t12.*3.0;
t173 = qd2_2.*t6.*t7.*t8.*t10.*2.0;
t174 = qd1_2.*t2.*t3.*t11.*t12.*2.0;
t175 = qd2_3.*t2.*t9.*t10.*t11.*2.0;
t176 = qd1_2.*t6.*t8.*t9.*t11.*t12.*2.0;
t177 = qd2_2.*t6.*t9.*t10.*t11.*t13.*2.0;
t178 = qd2_3.*t2.*t4.*t7.*t9.*t12.*2.0;
t179 = qd1_2.*t3.*t5.*t6.*t8.*t10.*2.0;
t180 = qd2_2.*t4.*t6.*t7.*t9.*t12.*t13.*2.0;
t118 = t99+t100+t101+t102+t103+t104+t105+t106+t107+t108+t109+t110+t111+t112+t113+t114+t115+t116+t117-t171-t172-t173-t174-t175-t176-t177-t178-t179-t180;
t119 = t3.*t4.*t13.*t15.*3.0;
t120 = t4.*t6.*t9.*t14.*3.0;
t121 = t5.*t7.*t9.*t13.*t15.*3.0;
t122 = t2.*t3.*t4.*t8.*t14.*3.0;
t123 = t4.*t6.*t11.*t12.*t14.*3.0;
t124 = t2.*t5.*t11.*t13.*t14.*3.0;
t125 = qd2_2.*t8.*t9.*t10.*t11.*t15.*2.0;
t126 = qd2_3.*t3.*t10.*t11.*t13.*t15.*2.0;
t127 = qd1_2.*t3.*t6.*t11.*t12.*t14.*2.0;
t128 = qd2_3.*t6.*t9.*t10.*t11.*t14.*2.0;
t129 = qd2_3.*t5.*t9.*t12.*t13.*t15.*2.0;
t130 = qd1_2.*t2.*t7.*t12.*t13.*t14.*2.0;
t131 = qd2_2.*t4.*t11.*t12.*t13.*t15.*2.0;
t132 = t2.*t5.*t7.*t8.*t9.*t14.*3.0;
t133 = t2.*t3.*t10.*t11.*t13.*t14.*3.0;
t134 = t2.*t5.*t9.*t12.*t13.*t14.*3.0;
t135 = qd2_3.*t4.*t6.*t7.*t9.*t12.*t14.*2.0;
t136 = qd2_2.*t2.*t3.*t5.*t12.*t13.*t14.*2.0;
t137 = qd2_3.*t2.*t5.*t8.*t9.*t12.*t14.*2.0;
t138 = qd1_2.*t2.*t4.*t10.*t11.*t13.*t14.*2.0;
t139 = qd2_2.*t2.*t4.*t8.*t11.*t12.*t14.*2.0;
t140 = t2.*t3.*t4.*t7.*t12.*t13.*t14.*3.0;
t141 = qd1_2.*t4.*t7.*t9.*t10.*t13.*t15.*2.0;
t142 = qd2_2.*t4.*t7.*t8.*t9.*t12.*t15.*2.0;
t143 = qd2_3.*t3.*t4.*t7.*t12.*t13.*t15.*2.0;
t144 = qd2_3.*t2.*t3.*t8.*t10.*t11.*t14.*2.0;
t145 = qd1_2.*t2.*t4.*t7.*t8.*t9.*t10.*t14.*2.0;
t146 = qd2_3.*t2.*t3.*t4.*t7.*t8.*t12.*t14.*2.0;
t245 = t6.*t7.*t10.*t14.*3.0;
t246 = t5.*t8.*t11.*t15.*3.0;
t247 = qd1_2.*t7.*t8.*t12.*t15.*2.0;
t248 = qd2_2.*t7.*t10.*t13.*t15.*2.0;
t249 = t3.*t8.*t10.*t11.*t15.*3.0;
t250 = t3.*t5.*t6.*t7.*t14.*3.0;
t251 = t5.*t8.*t9.*t12.*t15.*3.0;
t252 = qd1_2.*t3.*t5.*t10.*t13.*t15.*2.0;
t253 = qd2_2.*t3.*t5.*t8.*t12.*t15.*2.0;
t254 = qd1_2.*t4.*t8.*t10.*t11.*t15.*2.0;
t255 = qd2_2.*t2.*t7.*t8.*t10.*t14.*2.0;
t256 = qd1_2.*t5.*t6.*t9.*t10.*t14.*2.0;
t257 = qd2_3.*t3.*t5.*t6.*t12.*t14.*2.0;
t258 = qd1_2.*t9.*t11.*t12.*t13.*t15.*2.0;
t259 = t3.*t4.*t7.*t8.*t12.*t15.*3.0;
t260 = qd1_2.*t2.*t8.*t9.*t11.*t12.*t14.*2.0;
t261 = qd2_2.*t2.*t9.*t10.*t11.*t13.*t14.*2.0;
t262 = qd1_2.*t3.*t4.*t6.*t7.*t10.*t14.*2.0;
t263 = qd1_2.*t2.*t3.*t5.*t8.*t10.*t14.*2.0;
t264 = qd2_2.*t2.*t4.*t7.*t9.*t12.*t13.*t14.*2.0;
t147 = t119+t120+t121+t122+t123+t124+t125+t126+t127+t128+t129+t130+t131+t132+t133+t134+t135+t136+t137+t138+t139+t140+t141+t142+t143+t144+t145+t146-t245-t246-t247-t248-t249-t250-t251-t252-t253-t254-t255-t256-t257-t258-t259-t260-t261-t262-t263-t264;
t148 = t2.*t3.*t4.*3.0;
t149 = t2.*t5.*t7.*t9.*3.0;
t150 = t4.*t6.*t8.*t9.*3.0;
t151 = qd2_3.*t2.*t3.*t10.*t11.*2.0;
t152 = qd2_3.*t2.*t5.*t9.*t12.*2.0;
t153 = t6.*t9.*t10.*t11.*t13.*3.0;
t154 = qd1_2.*t3.*t6.*t8.*t11.*t12.*2.0;
t155 = qd2_2.*t3.*t6.*t10.*t11.*t13.*2.0;
t156 = qd2_3.*t6.*t8.*t9.*t10.*t11.*2.0;
t157 = qd1_2.*t2.*t4.*t7.*t9.*t10.*2.0;
t158 = qd2_3.*t2.*t3.*t4.*t7.*t12.*2.0;
t159 = qd2_2.*t5.*t6.*t9.*t12.*t13.*2.0;
t160 = t4.*t6.*t7.*t9.*t12.*t13.*3.0;
t161 = qd2_2.*t3.*t4.*t6.*t7.*t12.*t13.*2.0;
t162 = qd2_3.*t4.*t6.*t7.*t8.*t9.*t12.*2.0;
t265 = qd1_2.*t2.*t3.*t5.*t10.*2.0;
t266 = qd1_2.*t2.*t9.*t11.*t12.*2.0;
t267 = t3.*t5.*t6.*t7.*t8.*3.0;
t268 = t3.*t5.*t6.*t12.*t13.*3.0;
t269 = qd1_2.*t5.*t6.*t8.*t9.*t10.*2.0;
t270 = qd2_3.*t3.*t5.*t6.*t8.*t12.*2.0;
t271 = qd1_2.*t3.*t4.*t6.*t7.*t8.*t10.*2.0;
t163 = t148+t149+t150+t151+t152+t153+t154+t155+t156+t157+t158+t159+t160+t161+t162-t265-t266-t267-t268-t269-t270-t271;
t164 = t9.*t10.*t11;
t165 = t4.*t7.*t9.*t12;
t166 = t164+t165-t3.*t5.*t12;
t170 = t19+t20+t21-t167-t168;
t181 = qd1_3.*t3.*t4;
t182 = qd1_2.*t3.*t5.*t12;
t183 = qd1_3.*t5.*t7.*t9;
t184 = t181+t182+t183-qd1_2.*t9.*t10.*t11-qd1_2.*t4.*t7.*t9.*t12;
t186 = t3.*t4.*(3.0./2.0);
t187 = t5.*t7.*t9.*(3.0./2.0);
t188 = qd1_2.*t4.*t7.*t9.*t10;
t295 = qd1_2.*t3.*t5.*t10;
t296 = qd1_2.*t9.*t11.*t12;
t189 = t186+t187+t188-t295-t296;
t190 = qd1_3.*t4.*t9;
t191 = qd1_2.*t3.*t10.*t11;
t192 = qd1_2.*t5.*t9.*t12;
t193 = qd1_2.*t3.*t4.*t7.*t12;
t272 = t2.*t3.*t4;
t273 = t2.*t5.*t7.*t9;
t274 = t4.*t6.*t8.*t9;
t275 = t272+t273+t274-t3.*t5.*t6.*t7.*t8;
t276 = t2.*t3.*t10.*t11.*(1.0./1.0e1);
t277 = t2.*t5.*t9.*t12.*(1.0./1.0e1);
t278 = t6.*t8.*t9.*t10.*t11.*(1.0./1.0e1);
t279 = t2.*t3.*t4.*t7.*t12.*(1.0./1.0e1);
t280 = t4.*t6.*t11.*t12.*t13.*(1.0./1.0e1);
t281 = t4.*t6.*t7.*t8.*t9.*t12.*(1.0./1.0e1);
t395 = t6.*t7.*t10.*t13.*(1.0./1.0e1);
t396 = t3.*t5.*t6.*t8.*t12.*(1.0./1.0e1);
t282 = t276+t277+t278+t279+t280+t281-t395-t396;
t283 = t5.*t8.*t11.*t14;
t284 = t4.*t6.*t9.*t15;
t285 = t2.*t3.*t4.*t8.*t15;
t286 = t2.*t5.*t11.*t13.*t15;
t287 = t2.*t5.*t7.*t8.*t9.*t15;
t288 = t283+t284+t285+t286+t287-t3.*t4.*t13.*t14-t3.*t5.*t6.*t7.*t15-t5.*t7.*t9.*t13.*t14;
t289 = t6.*t9.*t10.*t11.*t15.*(1.0./1.0e1);
t290 = t2.*t3.*t8.*t10.*t11.*t15.*(1.0./1.0e1);
t291 = t4.*t6.*t7.*t9.*t12.*t15.*(1.0./1.0e1);
t292 = t2.*t5.*t8.*t9.*t12.*t15.*(1.0./1.0e1);
t293 = t2.*t3.*t4.*t7.*t8.*t12.*t15.*(1.0./1.0e1);
t397 = t3.*t5.*t6.*t12.*t15.*(1.0./1.0e1);
t398 = t3.*t10.*t11.*t13.*t14.*(1.0./1.0e1);
t399 = t5.*t9.*t12.*t13.*t14.*(1.0./1.0e1);
t400 = t3.*t4.*t7.*t12.*t13.*t14.*(1.0./1.0e1);
t294 = t289+t290+t291+t292+t293-t397-t398-t399-t400;
t297 = t2.*t3.*t5.*t7;
t298 = t3.*t4.*t6.*t8;
t299 = t5.*t6.*t11.*t13;
t300 = t5.*t6.*t7.*t8.*t9;
t301 = t297+t298+t299+t300-t2.*t4.*t9;
t302 = t2.*t3.*t5.*t12.*(1.0./1.0e1);
t303 = t3.*t6.*t8.*t10.*t11.*(1.0./1.0e1);
t304 = t5.*t6.*t8.*t9.*t12.*(1.0./1.0e1);
t305 = t3.*t4.*t6.*t7.*t8.*t12.*(1.0./1.0e1);
t381 = t2.*t9.*t10.*t11.*(1.0./1.0e1);
t382 = t2.*t4.*t7.*t9.*t12.*(1.0./1.0e1);
t306 = t302+t303+t304+t305-t381-t382;
t307 = t3.*t4.*t13.*t15;
t308 = t4.*t6.*t9.*t14;
t309 = t5.*t7.*t9.*t13.*t15;
t310 = t2.*t3.*t4.*t8.*t14;
t311 = t2.*t5.*t11.*t13.*t14;
t312 = t2.*t5.*t7.*t8.*t9.*t14;
t313 = t307+t308+t309+t310+t311+t312-t5.*t8.*t11.*t15-t3.*t5.*t6.*t7.*t14;
t314 = t3.*t10.*t11.*t13.*t15.*(1.0./1.0e1);
t315 = t6.*t9.*t10.*t11.*t14.*(1.0./1.0e1);
t316 = t5.*t9.*t12.*t13.*t15.*(1.0./1.0e1);
t317 = t3.*t4.*t7.*t12.*t13.*t15.*(1.0./1.0e1);
t318 = t2.*t3.*t8.*t10.*t11.*t14.*(1.0./1.0e1);
t319 = t4.*t6.*t7.*t9.*t12.*t14.*(1.0./1.0e1);
t320 = t2.*t5.*t8.*t9.*t12.*t14.*(1.0./1.0e1);
t321 = t2.*t3.*t4.*t7.*t8.*t12.*t14.*(1.0./1.0e1);
t357 = t3.*t5.*t6.*t12.*t14.*(1.0./1.0e1);
t322 = t314+t315+t316+t317+t318+t319+t320+t321-t357;
t323 = t3.*t4.*t6.*t15;
t324 = t4.*t9.*t13.*t14;
t325 = t5.*t6.*t7.*t9.*t15;
t326 = t2.*t3.*t5.*t7.*t8.*t15;
t327 = t323+t324+t325+t326-t2.*t4.*t8.*t9.*t15-t3.*t5.*t7.*t13.*t14;
t328 = t7.*t8.*t10.*t14.*(1.0./1.0e1);
t329 = t3.*t6.*t10.*t11.*t15.*(1.0./1.0e1);
t330 = t2.*t7.*t10.*t13.*t15.*(1.0./1.0e1);
t331 = t5.*t6.*t9.*t12.*t15.*(1.0./1.0e1);
t332 = t9.*t10.*t11.*t13.*t14.*(1.0./1.0e1);
t333 = t3.*t4.*t6.*t7.*t12.*t15.*(1.0./1.0e1);
t334 = t2.*t3.*t5.*t8.*t12.*t15.*(1.0./1.0e1);
t335 = t4.*t7.*t9.*t12.*t13.*t14.*(1.0./1.0e1);
t416 = t3.*t5.*t12.*t13.*t14.*(1.0./1.0e1);
t417 = t4.*t8.*t11.*t12.*t14.*(1.0./1.0e1);
t418 = t2.*t8.*t9.*t10.*t11.*t15.*(1.0./1.0e1);
t419 = t2.*t4.*t11.*t12.*t13.*t15.*(1.0./1.0e1);
t420 = t2.*t4.*t7.*t8.*t9.*t12.*t15.*(1.0./1.0e1);
t336 = t328+t329+t330+t331+t332+t333+t334+t335-t416-t417-t418-t419-t420;
t337 = t3.*t4.*t6.*t14;
t338 = t3.*t5.*t7.*t13.*t15;
t339 = t5.*t6.*t7.*t9.*t14;
t340 = t2.*t3.*t5.*t7.*t8.*t14;
t341 = t337+t338+t339+t340-t4.*t9.*t13.*t15-t2.*t4.*t8.*t9.*t14;
t342 = t3.*t6.*t10.*t11.*t14.*(1.0./1.0e1);
t343 = t3.*t5.*t12.*t13.*t15.*(1.0./1.0e1);
t344 = t4.*t8.*t11.*t12.*t15.*(1.0./1.0e1);
t345 = t2.*t7.*t10.*t13.*t14.*(1.0./1.0e1);
t346 = t5.*t6.*t9.*t12.*t14.*(1.0./1.0e1);
t347 = t3.*t4.*t6.*t7.*t12.*t14.*(1.0./1.0e1);
t348 = t2.*t3.*t5.*t8.*t12.*t14.*(1.0./1.0e1);
t375 = t7.*t8.*t10.*t15.*(1.0./1.0e1);
t376 = t9.*t10.*t11.*t13.*t15.*(1.0./1.0e1);
t377 = t4.*t7.*t9.*t12.*t13.*t15.*(1.0./1.0e1);
t378 = t2.*t8.*t9.*t10.*t11.*t14.*(1.0./1.0e1);
t379 = t2.*t4.*t11.*t12.*t13.*t14.*(1.0./1.0e1);
t380 = t2.*t4.*t7.*t8.*t9.*t12.*t14.*(1.0./1.0e1);
t349 = t342+t343+t344+t345+t346+t347+t348-t375-t376-t377-t378-t379-t380;
t350 = t5.*t11.*(3.0./2.0);
t351 = qd1_2.*t7.*t12;
t352 = t3.*t10.*t11.*(3.0./4.0);
t353 = t5.*t9.*t12.*(3.0./4.0);
t354 = qd1_2.*t4.*t10.*t11;
t355 = t3.*t4.*t7.*t12.*(3.0./4.0);
t356 = t350+t351+t352+t353+t354+t355;
t358 = t3.*t6.*t10.*t11.*t14;
t359 = t3.*t5.*t12.*t13.*t15;
t360 = t4.*t8.*t11.*t12.*t15;
t361 = t2.*t7.*t10.*t13.*t14;
t362 = t5.*t6.*t9.*t12.*t14;
t363 = t3.*t4.*t6.*t7.*t12.*t14;
t364 = t2.*t3.*t5.*t8.*t12.*t14;
t365 = t358+t359+t360+t361+t362+t363+t364-t7.*t8.*t10.*t15-t9.*t10.*t11.*t13.*t15-t2.*t8.*t9.*t10.*t11.*t14-t2.*t4.*t11.*t12.*t13.*t14-t4.*t7.*t9.*t12.*t13.*t15-t2.*t4.*t7.*t8.*t9.*t12.*t14;
t366 = t3.*t10.*t11.*t13.*t15;
t367 = t6.*t9.*t10.*t11.*t14;
t368 = t5.*t9.*t12.*t13.*t15;
t369 = t3.*t4.*t7.*t12.*t13.*t15;
t370 = t2.*t3.*t8.*t10.*t11.*t14;
t371 = t4.*t6.*t7.*t9.*t12.*t14;
t372 = t2.*t5.*t8.*t9.*t12.*t14;
t373 = t2.*t3.*t4.*t7.*t8.*t12.*t14;
t374 = t366+t367+t368+t369+t370+t371+t372+t373-t3.*t5.*t6.*t12.*t14;
t383 = t2.*t3.*t10.*t11;
t384 = t2.*t5.*t9.*t12;
t385 = t6.*t8.*t9.*t10.*t11;
t386 = t2.*t3.*t4.*t7.*t12;
t387 = t4.*t6.*t11.*t12.*t13;
t388 = t4.*t6.*t7.*t8.*t9.*t12;
t389 = t383+t384+t385+t386+t387+t388-t6.*t7.*t10.*t13-t3.*t5.*t6.*t8.*t12;
t390 = t2.*t3.*t5.*t12;
t391 = t3.*t6.*t8.*t10.*t11;
t392 = t5.*t6.*t8.*t9.*t12;
t393 = t3.*t4.*t6.*t7.*t8.*t12;
t394 = t390+t391+t392+t393-t2.*t9.*t10.*t11-t2.*t4.*t7.*t9.*t12;
t401 = t7.*t8.*t10.*t14;
t402 = t3.*t6.*t10.*t11.*t15;
t403 = t2.*t7.*t10.*t13.*t15;
t404 = t5.*t6.*t9.*t12.*t15;
t405 = t9.*t10.*t11.*t13.*t14;
t406 = t3.*t4.*t6.*t7.*t12.*t15;
t407 = t2.*t3.*t5.*t8.*t12.*t15;
t408 = t4.*t7.*t9.*t12.*t13.*t14;
t409 = t401+t402+t403+t404+t405+t406+t407+t408-t3.*t5.*t12.*t13.*t14-t4.*t8.*t11.*t12.*t14-t2.*t8.*t9.*t10.*t11.*t15-t2.*t4.*t11.*t12.*t13.*t15-t2.*t4.*t7.*t8.*t9.*t12.*t15;
t410 = t6.*t9.*t10.*t11.*t15;
t411 = t2.*t3.*t8.*t10.*t11.*t15;
t412 = t4.*t6.*t7.*t9.*t12.*t15;
t413 = t2.*t5.*t8.*t9.*t12.*t15;
t414 = t2.*t3.*t4.*t7.*t8.*t12.*t15;
t415 = t410+t411+t412+t413+t414-t3.*t5.*t6.*t12.*t15-t3.*t10.*t11.*t13.*t14-t5.*t9.*t12.*t13.*t14-t3.*t4.*t7.*t12.*t13.*t14;
t421 = t9.*t11.*t12;
t422 = t3.*t5.*t10;
t423 = t421+t422-t4.*t7.*t9.*t10;
t424 = t5.*t9.*t10;
t425 = t3.*t4.*t7.*t10;
t426 = t424+t425-t3.*t11.*t12;
out1 = reshape([0.0,-dqt2.*(t18.*(t19+t20+t21-t4.*t9.*(3.0./2.0)-qd1_2.*t3.*t11.*t12).*(-1.0./2.0)+t18.*t26.*(1.0./2.0)-t56.*t327.*(1.0./2.0)+t98.*t288.*(1.0./2.0)-t77.*t313.*(1.0./2.0)+t118.*t275.*(1.0./2.0)+t163.*t301.*(1.0./2.0)-t147.*t341.*(1.0./2.0)+t5.*t11.*t24.*(1.0./2.0)),dqt2.*(t24.*(t7.*t10-t4.*t11.*t12).*(1.0./2.0)-t26.*t166.*(1.0./2.0)+t166.*t170.*(1.0./2.0)-t77.*t365.*(1.0./2.0)+t56.*t415.*(1.0./2.0)+t98.*t409.*(1.0./2.0)+t118.*t394.*(1.0./2.0)+t147.*t374.*(1.0./2.0)-t163.*t389.*(1.0./2.0)),-dqt2.*(t77.*(t6.*t7.*t12.*t14.*(3.0./2.0)+qd1_3.*t4.*t6.*t9.*t14+qd1_3.*t3.*t4.*t13.*t15-qd1_2.*t7.*t8.*t10.*t15-qd1_3.*t5.*t8.*t11.*t15+qd2_2.*t7.*t12.*t13.*t15+t4.*t6.*t10.*t11.*t14.*(3.0./2.0)-t5.*t8.*t9.*t10.*t15.*(3.0./2.0)+t3.*t8.*t11.*t12.*t15.*(3.0./2.0)+qd1_3.*t2.*t3.*t4.*t8.*t14-qd1_3.*t3.*t5.*t6.*t7.*t14+qd1_2.*t3.*t6.*t10.*t11.*t14+qd1_2.*t2.*t7.*t10.*t13.*t14+qd1_2.*t5.*t6.*t9.*t12.*t14+qd1_3.*t2.*t5.*t11.*t13.*t14+qd1_2.*t3.*t5.*t12.*t13.*t15-qd2_3.*t3.*t5.*t6.*t10.*t14+qd1_2.*t4.*t8.*t11.*t12.*t15+qd1_3.*t5.*t7.*t9.*t13.*t15-qd2_2.*t3.*t5.*t8.*t10.*t15+qd2_2.*t2.*t7.*t8.*t12.*t14-qd1_2.*t9.*t10.*t11.*t13.*t15+qd2_2.*t4.*t10.*t11.*t13.*t15+qd2_3.*t5.*t9.*t10.*t13.*t15-qd2_3.*t6.*t9.*t11.*t12.*t14-qd2_2.*t8.*t9.*t11.*t12.*t15-qd2_3.*t3.*t11.*t12.*t13.*t15-t3.*t4.*t7.*t8.*t10.*t15.*(3.0./2.0)+t2.*t5.*t9.*t10.*t13.*t14.*(3.0./2.0)-t2.*t3.*t11.*t12.*t13.*t14.*(3.0./2.0)+qd1_2.*t2.*t3.*t5.*t8.*t12.*t14+qd1_2.*t3.*t4.*t6.*t7.*t12.*t14+qd1_3.*t2.*t5.*t7.*t8.*t9.*t14-qd1_2.*t2.*t8.*t9.*t10.*t11.*t14-qd1_2.*t2.*t4.*t11.*t12.*t13.*t14+qd2_2.*t2.*t3.*t5.*t10.*t13.*t14+qd2_2.*t2.*t4.*t8.*t10.*t11.*t14+qd2_3.*t2.*t5.*t8.*t9.*t10.*t14-qd1_2.*t4.*t7.*t9.*t12.*t13.*t15-qd2_3.*t2.*t3.*t8.*t11.*t12.*t14+qd2_3.*t4.*t6.*t7.*t9.*t10.*t14+qd2_2.*t4.*t7.*t8.*t9.*t10.*t15+qd2_3.*t3.*t4.*t7.*t10.*t13.*t15+qd2_2.*t2.*t9.*t11.*t12.*t13.*t14+t2.*t3.*t4.*t7.*t10.*t13.*t14.*(3.0./2.0)-qd1_2.*t2.*t4.*t7.*t8.*t9.*t12.*t14+qd2_3.*t2.*t3.*t4.*t7.*t8.*t10.*t14-qd2_2.*t2.*t4.*t7.*t9.*t10.*t13.*t14).*(1.0./2.0)-t98.*(t6.*t7.*t12.*t15.*(3.0./2.0)-qd1_3.*t3.*t4.*t13.*t14+qd1_3.*t4.*t6.*t9.*t15+qd1_2.*t7.*t8.*t10.*t14+qd1_3.*t5.*t8.*t11.*t14-qd2_2.*t7.*t12.*t13.*t14+t4.*t6.*t10.*t11.*t15.*(3.0./2.0)+t5.*t8.*t9.*t10.*t14.*(3.0./2.0)-t3.*t8.*t11.*t12.*t14.*(3.0./2.0)+qd1_3.*t2.*t3.*t4.*t8.*t15-qd1_3.*t3.*t5.*t6.*t7.*t15+qd1_2.*t3.*t6.*t10.*t11.*t15+qd1_2.*t2.*t7.*t10.*t13.*t15-qd1_2.*t3.*t5.*t12.*t13.*t14+qd1_2.*t5.*t6.*t9.*t12.*t15+qd1_3.*t2.*t5.*t11.*t13.*t15-qd1_2.*t4.*t8.*t11.*t12.*t14-qd1_3.*t5.*t7.*t9.*t13.*t14+qd2_2.*t3.*t5.*t8.*t10.*t14-qd2_3.*t3.*t5.*t6.*t10.*t15+qd2_2.*t2.*t7.*t8.*t12.*t15+qd1_2.*t9.*t10.*t11.*t13.*t14-qd2_2.*t4.*t10.*t11.*t13.*t14-qd2_3.*t5.*t9.*t10.*t13.*t14+qd2_2.*t8.*t9.*t11.*t12.*t14+qd2_3.*t3.*t11.*t12.*t13.*t14-qd2_3.*t6.*t9.*t11.*t12.*t15+t3.*t4.*t7.*t8.*t10.*t14.*(3.0./2.0)+t2.*t5.*t9.*t10.*t13.*t15.*(3.0./2.0)-t2.*t3.*t11.*t12.*t13.*t15.*(3.0./2.0)+qd1_2.*t2.*t3.*t5.*t8.*t12.*t15+qd1_2.*t3.*t4.*t6.*t7.*t12.*t15+qd1_3.*t2.*t5.*t7.*t8.*t9.*t15-qd1_2.*t2.*t8.*t9.*t10.*t11.*t15-qd1_2.*t2.*t4.*t11.*t12.*t13.*t15+qd2_2.*t2.*t3.*t5.*t10.*t13.*t15+qd1_2.*t4.*t7.*t9.*t12.*t13.*t14+qd2_2.*t2.*t4.*t8.*t10.*t11.*t15+qd2_3.*t2.*t5.*t8.*t9.*t10.*t15-qd2_2.*t4.*t7.*t8.*t9.*t10.*t14-qd2_3.*t2.*t3.*t8.*t11.*t12.*t15-qd2_3.*t3.*t4.*t7.*t10.*t13.*t14+qd2_3.*t4.*t6.*t7.*t9.*t10.*t15+qd2_2.*t2.*t9.*t11.*t12.*t13.*t15+t2.*t3.*t4.*t7.*t10.*t13.*t15.*(3.0./2.0)-qd1_2.*t2.*t4.*t7.*t8.*t9.*t12.*t15+qd2_3.*t2.*t3.*t4.*t7.*t8.*t10.*t15-qd2_2.*t2.*t4.*t7.*t9.*t10.*t13.*t15).*(1.0./2.0)-t189.*(t190+t191+t192+t193-qd1_3.*t3.*t5.*t7).*(1.0./2.0)+t189.*(t190+t191+t192+t193+t7.*t12.*(3.0./4.0)+t4.*t10.*t11.*(3.0./4.0)-qd1_3.*t3.*t5.*t7).*(1.0./2.0)+t294.*(t7.*t8.*t12.*t14+t3.*t5.*t10.*t13.*t14-t5.*t6.*t9.*t10.*t15+t3.*t6.*t11.*t12.*t15+t4.*t8.*t10.*t11.*t14+t2.*t7.*t12.*t13.*t15+t9.*t11.*t12.*t13.*t14-t2.*t3.*t5.*t8.*t10.*t15-t3.*t4.*t6.*t7.*t10.*t15+t2.*t4.*t10.*t11.*t13.*t15-t2.*t8.*t9.*t11.*t12.*t15-t4.*t7.*t9.*t10.*t13.*t14+t2.*t4.*t7.*t8.*t9.*t10.*t15).*(1.0./2.0)-t322.*(t7.*t8.*t12.*t15+t5.*t6.*t9.*t10.*t14+t3.*t5.*t10.*t13.*t15-t3.*t6.*t11.*t12.*t14-t2.*t7.*t12.*t13.*t14+t4.*t8.*t10.*t11.*t15+t9.*t11.*t12.*t13.*t15+t2.*t3.*t5.*t8.*t10.*t14+t3.*t4.*t6.*t7.*t10.*t14-t2.*t4.*t10.*t11.*t13.*t14+t2.*t8.*t9.*t11.*t12.*t14-t4.*t7.*t9.*t10.*t13.*t15-t2.*t4.*t7.*t8.*t9.*t10.*t14).*(1.0./2.0)-t26.*t184.*(1.0./2.0)+t170.*t184.*(1.0./2.0)+t306.*(t2.*t5.*t9.*t10-t2.*t3.*t11.*t12+t6.*t7.*t12.*t13+t2.*t3.*t4.*t7.*t10-t3.*t5.*t6.*t8.*t10+t4.*t6.*t10.*t11.*t13-t6.*t8.*t9.*t11.*t12+t4.*t6.*t7.*t8.*t9.*t10).*(1.0./2.0)-t24.*(qd1_2.*t7.*t10+qd1_3.*t5.*t11+t5.*t9.*t10.*(3.0./4.0)-t3.*t11.*t12.*(3.0./4.0)-qd1_2.*t4.*t11.*t12+t3.*t4.*t7.*t10.*(3.0./4.0)).*(1.0./2.0)+t282.*(t2.*t3.*t5.*t10+t2.*t9.*t11.*t12-t2.*t4.*t7.*t9.*t10+t5.*t6.*t8.*t9.*t10-t3.*t6.*t8.*t11.*t12+t3.*t4.*t6.*t7.*t8.*t10).*(1.0./2.0)+t336.*(t3.*t5.*t6.*t10.*t15+t5.*t9.*t10.*t13.*t14-t3.*t11.*t12.*t13.*t14+t6.*t9.*t11.*t12.*t15-t2.*t5.*t8.*t9.*t10.*t15+t2.*t3.*t8.*t11.*t12.*t15+t3.*t4.*t7.*t10.*t13.*t14-t4.*t6.*t7.*t9.*t10.*t15-t2.*t3.*t4.*t7.*t8.*t10.*t15).*(1.0./2.0)-t349.*(-t3.*t5.*t6.*t10.*t14+t5.*t9.*t10.*t13.*t15-t6.*t9.*t11.*t12.*t14-t3.*t11.*t12.*t13.*t15+t2.*t5.*t8.*t9.*t10.*t14-t2.*t3.*t8.*t11.*t12.*t14+t4.*t6.*t7.*t9.*t10.*t14+t3.*t4.*t7.*t10.*t13.*t15+t2.*t3.*t4.*t7.*t8.*t10.*t14).*(1.0./2.0)-t163.*(t2.*t7.*t12.*(-3.0./2.0)-qd1_3.*t2.*t4.*t9-t2.*t4.*t10.*t11.*(3.0./2.0)+qd1_3.*t2.*t3.*t5.*t7+qd1_3.*t3.*t4.*t6.*t8-qd1_2.*t2.*t3.*t10.*t11-qd1_2.*t2.*t5.*t9.*t12+qd2_3.*t2.*t3.*t5.*t10+qd1_2.*t6.*t7.*t10.*t13+qd1_3.*t5.*t6.*t11.*t13+qd2_2.*t6.*t7.*t8.*t12+qd2_3.*t2.*t9.*t11.*t12+t5.*t6.*t9.*t10.*t13.*(3.0./2.0)-t3.*t6.*t11.*t12.*t13.*(3.0./2.0)-qd1_2.*t2.*t3.*t4.*t7.*t12+qd1_2.*t3.*t5.*t6.*t8.*t12+qd1_3.*t5.*t6.*t7.*t8.*t9-qd2_3.*t2.*t4.*t7.*t9.*t10-qd1_2.*t6.*t8.*t9.*t10.*t11-qd1_2.*t4.*t6.*t11.*t12.*t13+qd2_2.*t3.*t5.*t6.*t10.*t13+qd2_2.*t4.*t6.*t8.*t10.*t11+qd2_3.*t5.*t6.*t8.*t9.*t10-qd2_3.*t3.*t6.*t8.*t11.*t12+qd2_2.*t6.*t9.*t11.*t12.*t13+t3.*t4.*t6.*t7.*t10.*t13.*(3.0./2.0)-qd1_2.*t4.*t6.*t7.*t8.*t9.*t12+qd2_3.*t3.*t4.*t6.*t7.*t8.*t10-qd2_2.*t4.*t6.*t7.*t9.*t10.*t13).*(1.0./2.0)-t118.*(qd1_3.*t2.*t3.*t4+qd1_2.*t2.*t3.*t5.*t12+qd1_3.*t2.*t5.*t7.*t9+qd1_3.*t4.*t6.*t8.*t9-qd1_2.*t2.*t9.*t10.*t11+qd2_3.*t2.*t5.*t9.*t10-qd2_3.*t2.*t3.*t11.*t12-t3.*t5.*t6.*t10.*t13.*(3.0./2.0)-t6.*t9.*t11.*t12.*t13.*(3.0./2.0)-qd1_3.*t3.*t5.*t6.*t7.*t8-qd1_2.*t2.*t4.*t7.*t9.*t12+qd2_3.*t2.*t3.*t4.*t7.*t10+qd1_2.*t3.*t6.*t8.*t10.*t11+qd1_2.*t5.*t6.*t8.*t9.*t12-qd2_3.*t3.*t5.*t6.*t8.*t10+qd2_2.*t5.*t6.*t9.*t10.*t13-qd2_2.*t3.*t6.*t11.*t12.*t13-qd2_3.*t6.*t8.*t9.*t11.*t12+t4.*t6.*t7.*t9.*t10.*t13.*(3.0./2.0)+qd1_2.*t3.*t4.*t6.*t7.*t8.*t12+qd2_2.*t3.*t4.*t6.*t7.*t10.*t13+qd2_3.*t4.*t6.*t7.*t8.*t9.*t10).*(1.0./2.0)+t356.*(t3.*t5.*t10.*(3.0./4.0)+t9.*t11.*t12.*(3.0./4.0)-t4.*t7.*t9.*t10.*(3.0./4.0)).*(1.0./2.0)+t56.*(qd1_3.*t3.*t4.*t6.*t15+qd1_3.*t4.*t9.*t13.*t14+t3.*t5.*t8.*t10.*t14.*(3.0./2.0)+t8.*t9.*t11.*t12.*t14.*(3.0./2.0)-qd1_3.*t2.*t4.*t8.*t9.*t15+qd1_2.*t3.*t5.*t6.*t12.*t15-qd1_3.*t3.*t5.*t7.*t13.*t14+qd1_3.*t5.*t6.*t7.*t9.*t15+qd1_2.*t3.*t10.*t11.*t13.*t14-qd1_2.*t6.*t9.*t10.*t11.*t15+qd1_2.*t5.*t9.*t12.*t13.*t14-qd2_2.*t5.*t8.*t9.*t10.*t14-qd2_3.*t3.*t5.*t10.*t13.*t14+qd2_3.*t5.*t6.*t9.*t10.*t15+qd2_2.*t3.*t8.*t11.*t12.*t14-qd2_3.*t3.*t6.*t11.*t12.*t15-qd2_3.*t9.*t11.*t12.*t13.*t14+t2.*t3.*t5.*t10.*t13.*t15.*(3.0./2.0)-t4.*t7.*t8.*t9.*t10.*t14.*(3.0./2.0)+t2.*t9.*t11.*t12.*t13.*t15.*(3.0./2.0)+qd1_3.*t2.*t3.*t5.*t7.*t8.*t15-qd1_2.*t2.*t3.*t8.*t10.*t11.*t15-qd1_2.*t2.*t5.*t8.*t9.*t12.*t15+qd1_2.*t3.*t4.*t7.*t12.*t13.*t14-qd1_2.*t4.*t6.*t7.*t9.*t12.*t15+qd2_3.*t2.*t3.*t5.*t8.*t10.*t15-qd2_2.*t3.*t4.*t7.*t8.*t10.*t14+qd2_3.*t3.*t4.*t6.*t7.*t10.*t15-qd2_2.*t2.*t5.*t9.*t10.*t13.*t15+qd2_2.*t2.*t3.*t11.*t12.*t13.*t15+qd2_3.*t2.*t8.*t9.*t11.*t12.*t15+qd2_3.*t4.*t7.*t9.*t10.*t13.*t14-t2.*t4.*t7.*t9.*t10.*t13.*t15.*(3.0./2.0)-qd1_2.*t2.*t3.*t4.*t7.*t8.*t12.*t15-qd2_2.*t2.*t3.*t4.*t7.*t10.*t13.*t15-qd2_3.*t2.*t4.*t7.*t8.*t9.*t10.*t15).*(1.0./2.0)-t147.*(-qd1_3.*t3.*t4.*t6.*t14+qd1_3.*t4.*t9.*t13.*t15+t3.*t5.*t8.*t10.*t15.*(3.0./2.0)+t8.*t9.*t11.*t12.*t15.*(3.0./2.0)+qd1_3.*t2.*t4.*t8.*t9.*t14-qd1_2.*t3.*t5.*t6.*t12.*t14-qd1_3.*t5.*t6.*t7.*t9.*t14-qd1_3.*t3.*t5.*t7.*t13.*t15+qd1_2.*t6.*t9.*t10.*t11.*t14+qd1_2.*t3.*t10.*t11.*t13.*t15+qd1_2.*t5.*t9.*t12.*t13.*t15-qd2_3.*t5.*t6.*t9.*t10.*t14-qd2_2.*t5.*t8.*t9.*t10.*t15-qd2_3.*t3.*t5.*t10.*t13.*t15+qd2_3.*t3.*t6.*t11.*t12.*t14+qd2_2.*t3.*t8.*t11.*t12.*t15-qd2_3.*t9.*t11.*t12.*t13.*t15-t2.*t3.*t5.*t10.*t13.*t14.*(3.0./2.0)-t4.*t7.*t8.*t9.*t10.*t15.*(3.0./2.0)-t2.*t9.*t11.*t12.*t13.*t14.*(3.0./2.0)-qd1_3.*t2.*t3.*t5.*t7.*t8.*t14+qd1_2.*t2.*t3.*t8.*t10.*t11.*t14+qd1_2.*t2.*t5.*t8.*t9.*t12.*t14+qd1_2.*t4.*t6.*t7.*t9.*t12.*t14-qd2_3.*t2.*t3.*t5.*t8.*t10.*t14+qd1_2.*t3.*t4.*t7.*t12.*t13.*t15-qd2_3.*t3.*t4.*t6.*t7.*t10.*t14-qd2_2.*t3.*t4.*t7.*t8.*t10.*t15+qd2_2.*t2.*t5.*t9.*t10.*t13.*t14-qd2_2.*t2.*t3.*t11.*t12.*t13.*t14-qd2_3.*t2.*t8.*t9.*t11.*t12.*t14+qd2_3.*t4.*t7.*t9.*t10.*t13.*t15+t2.*t4.*t7.*t9.*t10.*t13.*t14.*(3.0./2.0)+qd1_2.*t2.*t3.*t4.*t7.*t8.*t12.*t14+qd2_2.*t2.*t3.*t4.*t7.*t10.*t13.*t14+qd2_3.*t2.*t4.*t7.*t8.*t9.*t10.*t14).*(1.0./2.0)),dqt2.*(t24.*(t4.*t9.*(3.0./4.0)+qd1_3.*t7.*t12-t3.*t5.*t7.*(3.0./4.0)+qd1_3.*t4.*t10.*t11).*(1.0./2.0)-t26.*(qd1_3.*t3.*t5.*t10+qd1_3.*t9.*t11.*t12-qd1_3.*t4.*t7.*t9.*t10).*(1.0./2.0)-t356.*(t3.*t4.*(3.0./4.0)+t5.*t7.*t9.*(3.0./4.0)).*(1.0./2.0)-t275.*t282.*(1.0./2.0)+t288.*t294.*(1.0./2.0)+t301.*t306.*(1.0./2.0)+t313.*t322.*(1.0./2.0)-t327.*t336.*(1.0./2.0)-t341.*t349.*(1.0./2.0)-t118.*(-qd2_3.*t2.*t4.*t9+t3.*t4.*t6.*t13.*(3.0./2.0)+qd1_3.*t2.*t3.*t5.*t10+qd2_3.*t2.*t3.*t5.*t7+qd2_3.*t3.*t4.*t6.*t8+qd1_3.*t2.*t9.*t11.*t12-qd2_2.*t4.*t6.*t9.*t13+t5.*t6.*t7.*t9.*t13.*(3.0./2.0)-qd1_3.*t2.*t4.*t7.*t9.*t10+qd1_3.*t5.*t6.*t8.*t9.*t10-qd1_3.*t3.*t6.*t8.*t11.*t12+qd2_2.*t3.*t5.*t6.*t7.*t13+qd2_3.*t5.*t6.*t7.*t8.*t9+qd1_3.*t3.*t4.*t6.*t7.*t8.*t10).*(1.0./2.0)+t163.*(t2.*t5.*t11.*(3.0./2.0)+qd2_3.*t2.*t3.*t4+t4.*t6.*t9.*t13.*(3.0./2.0)+qd1_3.*t2.*t5.*t9.*t10-qd1_3.*t2.*t3.*t11.*t12+qd2_3.*t2.*t5.*t7.*t9+qd2_2.*t3.*t4.*t6.*t13+qd2_3.*t4.*t6.*t8.*t9+qd1_3.*t6.*t7.*t12.*t13-qd2_2.*t5.*t6.*t8.*t11-t3.*t5.*t6.*t7.*t13.*(3.0./2.0)+qd1_3.*t2.*t3.*t4.*t7.*t10-qd1_3.*t3.*t5.*t6.*t8.*t10-qd2_3.*t3.*t5.*t6.*t7.*t8+qd1_3.*t4.*t6.*t10.*t11.*t13-qd1_3.*t6.*t8.*t9.*t11.*t12+qd2_2.*t5.*t6.*t7.*t9.*t13+qd1_3.*t4.*t6.*t7.*t8.*t9.*t10).*(1.0./2.0)+t189.*(t5.*t11.*(3.0./4.0)+qd1_3.*t5.*t9.*t10-qd1_3.*t3.*t11.*t12+qd1_3.*t3.*t4.*t7.*t10).*(1.0./2.0)-t56.*(t3.*t4.*t8.*t14.*(3.0./2.0)-qd2_2.*t4.*t8.*t9.*t14-qd2_3.*t3.*t4.*t13.*t14+qd2_3.*t4.*t6.*t9.*t15+t2.*t3.*t4.*t13.*t15.*(3.0./2.0)+t5.*t7.*t8.*t9.*t14.*(3.0./2.0)-qd1_3.*t3.*t5.*t6.*t10.*t15+qd2_3.*t2.*t3.*t4.*t8.*t15+qd2_2.*t3.*t5.*t7.*t8.*t14-qd2_3.*t3.*t5.*t6.*t7.*t15-qd1_3.*t5.*t9.*t10.*t13.*t14-qd2_2.*t2.*t4.*t9.*t13.*t15+qd1_3.*t3.*t11.*t12.*t13.*t14-qd1_3.*t6.*t9.*t11.*t12.*t15-qd2_3.*t5.*t7.*t9.*t13.*t14+t2.*t5.*t7.*t9.*t13.*t15.*(3.0./2.0)+qd1_3.*t2.*t5.*t8.*t9.*t10.*t15-qd1_3.*t2.*t3.*t8.*t11.*t12.*t15-qd1_3.*t3.*t4.*t7.*t10.*t13.*t14+qd1_3.*t4.*t6.*t7.*t9.*t10.*t15+qd2_2.*t2.*t3.*t5.*t7.*t13.*t15+qd2_3.*t2.*t5.*t7.*t8.*t9.*t15+qd1_3.*t2.*t3.*t4.*t7.*t8.*t10.*t15).*(1.0./2.0)-t147.*(t3.*t4.*t8.*t15.*(-3.0./2.0)+qd2_3.*t4.*t6.*t9.*t14+qd2_2.*t4.*t8.*t9.*t15+qd2_3.*t3.*t4.*t13.*t15+t2.*t3.*t4.*t13.*t14.*(3.0./2.0)-t5.*t7.*t8.*t9.*t15.*(3.0./2.0)-qd1_3.*t3.*t5.*t6.*t10.*t14+qd2_3.*t2.*t3.*t4.*t8.*t14-qd2_3.*t3.*t5.*t6.*t7.*t14-qd2_2.*t3.*t5.*t7.*t8.*t15-qd2_2.*t2.*t4.*t9.*t13.*t14+qd1_3.*t5.*t9.*t10.*t13.*t15-qd1_3.*t6.*t9.*t11.*t12.*t14-qd1_3.*t3.*t11.*t12.*t13.*t15+qd2_3.*t5.*t7.*t9.*t13.*t15+t2.*t5.*t7.*t9.*t13.*t14.*(3.0./2.0)+qd1_3.*t2.*t5.*t8.*t9.*t10.*t14-qd1_3.*t2.*t3.*t8.*t11.*t12.*t14+qd1_3.*t4.*t6.*t7.*t9.*t10.*t14+qd1_3.*t3.*t4.*t7.*t10.*t13.*t15+qd2_2.*t2.*t3.*t5.*t7.*t13.*t14+qd2_3.*t2.*t5.*t7.*t8.*t9.*t14+qd1_3.*t2.*t3.*t4.*t7.*t8.*t10.*t14).*(1.0./2.0)+t77.*(t4.*t8.*t9.*t15.*(3.0./2.0)+t5.*t6.*t11.*t14.*(3.0./2.0)+qd2_3.*t3.*t4.*t6.*t14+qd2_2.*t3.*t4.*t8.*t15+qd1_3.*t7.*t8.*t12.*t15-qd2_3.*t4.*t9.*t13.*t15+qd2_2.*t5.*t11.*t13.*t15-t3.*t5.*t7.*t8.*t15.*(3.0./2.0)-t2.*t4.*t9.*t13.*t14.*(3.0./2.0)+qd1_3.*t5.*t6.*t9.*t10.*t14-qd2_2.*t2.*t3.*t4.*t13.*t14+qd1_3.*t3.*t5.*t10.*t13.*t15-qd1_3.*t3.*t6.*t11.*t12.*t14-qd2_3.*t2.*t4.*t8.*t9.*t14-qd1_3.*t2.*t7.*t12.*t13.*t14+qd1_3.*t4.*t8.*t10.*t11.*t15+qd2_2.*t2.*t5.*t8.*t11.*t14+qd2_3.*t5.*t6.*t7.*t9.*t14+qd2_2.*t5.*t7.*t8.*t9.*t15+qd2_3.*t3.*t5.*t7.*t13.*t15+qd1_3.*t9.*t11.*t12.*t13.*t15+t2.*t3.*t5.*t7.*t13.*t14.*(3.0./2.0)+qd1_3.*t2.*t3.*t5.*t8.*t10.*t14+qd1_3.*t3.*t4.*t6.*t7.*t10.*t14+qd2_3.*t2.*t3.*t5.*t7.*t8.*t14-qd1_3.*t2.*t4.*t10.*t11.*t13.*t14+qd1_3.*t2.*t8.*t9.*t11.*t12.*t14-qd1_3.*t4.*t7.*t9.*t10.*t13.*t15-qd2_2.*t2.*t5.*t7.*t9.*t13.*t14-qd1_3.*t2.*t4.*t7.*t8.*t9.*t10.*t14).*(1.0./2.0)+t98.*(t4.*t8.*t9.*t14.*(3.0./2.0)-t5.*t6.*t11.*t15.*(3.0./2.0)+qd2_2.*t3.*t4.*t8.*t14-qd2_3.*t3.*t4.*t6.*t15+qd1_3.*t7.*t8.*t12.*t14-qd2_3.*t4.*t9.*t13.*t14+qd2_2.*t5.*t11.*t13.*t14-t3.*t5.*t7.*t8.*t14.*(3.0./2.0)+t2.*t4.*t9.*t13.*t15.*(3.0./2.0)+qd1_3.*t3.*t5.*t10.*t13.*t14-qd1_3.*t5.*t6.*t9.*t10.*t15+qd2_2.*t2.*t3.*t4.*t13.*t15+qd1_3.*t3.*t6.*t11.*t12.*t15+qd1_3.*t4.*t8.*t10.*t11.*t14+qd2_3.*t2.*t4.*t8.*t9.*t15+qd1_3.*t2.*t7.*t12.*t13.*t15-qd2_2.*t2.*t5.*t8.*t11.*t15+qd2_2.*t5.*t7.*t8.*t9.*t14+qd2_3.*t3.*t5.*t7.*t13.*t14-qd2_3.*t5.*t6.*t7.*t9.*t15+qd1_3.*t9.*t11.*t12.*t13.*t14-t2.*t3.*t5.*t7.*t13.*t15.*(3.0./2.0)-qd1_3.*t2.*t3.*t5.*t8.*t10.*t15-qd1_3.*t3.*t4.*t6.*t7.*t10.*t15-qd2_3.*t2.*t3.*t5.*t7.*t8.*t15+qd1_3.*t2.*t4.*t10.*t11.*t13.*t15-qd1_3.*t2.*t8.*t9.*t11.*t12.*t15-qd1_3.*t4.*t7.*t9.*t10.*t13.*t14+qd2_2.*t2.*t5.*t7.*t9.*t13.*t15+qd1_3.*t2.*t4.*t7.*t8.*t9.*t10.*t15).*(1.0./2.0)+qd1_3.*t170.*t423.*(1.0./2.0)-qd1_3.*t189.*t426.*(1.0./2.0)),-dqt2.*(t189.*(t20+t21+t25-t168-t185).*(1.0./2.0)+t24.*t356.*(1.0./2.0)-t282.*t394.*(1.0./2.0)+t322.*t365.*(1.0./2.0)-t306.*t389.*(1.0./2.0)+t294.*t409.*(1.0./2.0)+t349.*t374.*(1.0./2.0)+t336.*t415.*(1.0./2.0)+t118.*(-qd1_2.*t2.*t3.*t5.*t10-qd1_2.*t2.*t9.*t11.*t12+qd2_3.*t2.*t3.*t10.*t11+qd2_3.*t2.*t5.*t9.*t12-t3.*t5.*t6.*t12.*t13.*(3.0./2.0)+t6.*t9.*t10.*t11.*t13.*(3.0./2.0)+qd1_2.*t2.*t4.*t7.*t9.*t10-qd1_2.*t5.*t6.*t8.*t9.*t10+qd2_3.*t2.*t3.*t4.*t7.*t12+qd1_2.*t3.*t6.*t8.*t11.*t12-qd2_3.*t3.*t5.*t6.*t8.*t12+qd2_2.*t3.*t6.*t10.*t11.*t13+qd2_2.*t5.*t6.*t9.*t12.*t13+qd2_3.*t6.*t8.*t9.*t10.*t11+t4.*t6.*t7.*t9.*t12.*t13.*(3.0./2.0)-qd1_2.*t3.*t4.*t6.*t7.*t8.*t10+qd2_2.*t3.*t4.*t6.*t7.*t12.*t13+qd2_3.*t4.*t6.*t7.*t8.*t9.*t12).*(1.0./2.0)-t77.*(t6.*t7.*t10.*t14.*(-3.0./2.0)-qd1_2.*t7.*t8.*t12.*t15-qd2_2.*t7.*t10.*t13.*t15-t3.*t8.*t10.*t11.*t15.*(3.0./2.0)+t4.*t6.*t11.*t12.*t14.*(3.0./2.0)-t5.*t8.*t9.*t12.*t15.*(3.0./2.0)-qd1_2.*t5.*t6.*t9.*t10.*t14-qd1_2.*t3.*t5.*t10.*t13.*t15+qd1_2.*t3.*t6.*t11.*t12.*t14+qd1_2.*t2.*t7.*t12.*t13.*t14-qd1_2.*t4.*t8.*t10.*t11.*t15-qd2_2.*t2.*t7.*t8.*t10.*t14-qd2_3.*t3.*t5.*t6.*t12.*t14-qd2_2.*t3.*t5.*t8.*t12.*t15-qd1_2.*t9.*t11.*t12.*t13.*t15+qd2_3.*t6.*t9.*t10.*t11.*t14+qd2_2.*t8.*t9.*t10.*t11.*t15+qd2_3.*t3.*t10.*t11.*t13.*t15+qd2_2.*t4.*t11.*t12.*t13.*t15+qd2_3.*t5.*t9.*t12.*t13.*t15-t3.*t4.*t7.*t8.*t12.*t15.*(3.0./2.0)+t2.*t3.*t10.*t11.*t13.*t14.*(3.0./2.0)+t2.*t5.*t9.*t12.*t13.*t14.*(3.0./2.0)-qd1_2.*t2.*t3.*t5.*t8.*t10.*t14-qd1_2.*t3.*t4.*t6.*t7.*t10.*t14+qd1_2.*t2.*t4.*t10.*t11.*t13.*t14-qd1_2.*t2.*t8.*t9.*t11.*t12.*t14+qd1_2.*t4.*t7.*t9.*t10.*t13.*t15+qd2_2.*t2.*t3.*t5.*t12.*t13.*t14+qd2_3.*t2.*t3.*t8.*t10.*t11.*t14+qd2_2.*t2.*t4.*t8.*t11.*t12.*t14+qd2_3.*t2.*t5.*t8.*t9.*t12.*t14+qd2_3.*t4.*t6.*t7.*t9.*t12.*t14+qd2_2.*t4.*t7.*t8.*t9.*t12.*t15+qd2_3.*t3.*t4.*t7.*t12.*t13.*t15-qd2_2.*t2.*t9.*t10.*t11.*t13.*t14+t2.*t3.*t4.*t7.*t12.*t13.*t14.*(3.0./2.0)+qd1_2.*t2.*t4.*t7.*t8.*t9.*t10.*t14+qd2_3.*t2.*t3.*t4.*t7.*t8.*t12.*t14-qd2_2.*t2.*t4.*t7.*t9.*t12.*t13.*t14).*(1.0./2.0)+t98.*(t6.*t7.*t10.*t15.*(-3.0./2.0)+qd1_2.*t7.*t8.*t12.*t14+qd2_2.*t7.*t10.*t13.*t14+t3.*t8.*t10.*t11.*t14.*(3.0./2.0)+t4.*t6.*t11.*t12.*t15.*(3.0./2.0)+t5.*t8.*t9.*t12.*t14.*(3.0./2.0)+qd1_2.*t3.*t5.*t10.*t13.*t14-qd1_2.*t5.*t6.*t9.*t10.*t15+qd1_2.*t3.*t6.*t11.*t12.*t15+qd1_2.*t4.*t8.*t10.*t11.*t14+qd1_2.*t2.*t7.*t12.*t13.*t15-qd2_2.*t2.*t7.*t8.*t10.*t15+qd2_2.*t3.*t5.*t8.*t12.*t14-qd2_3.*t3.*t5.*t6.*t12.*t15+qd1_2.*t9.*t11.*t12.*t13.*t14-qd2_2.*t8.*t9.*t10.*t11.*t14-qd2_3.*t3.*t10.*t11.*t13.*t14+qd2_3.*t6.*t9.*t10.*t11.*t15-qd2_2.*t4.*t11.*t12.*t13.*t14-qd2_3.*t5.*t9.*t12.*t13.*t14+t3.*t4.*t7.*t8.*t12.*t14.*(3.0./2.0)+t2.*t3.*t10.*t11.*t13.*t15.*(3.0./2.0)+t2.*t5.*t9.*t12.*t13.*t15.*(3.0./2.0)-qd1_2.*t2.*t3.*t5.*t8.*t10.*t15-qd1_2.*t3.*t4.*t6.*t7.*t10.*t15+qd1_2.*t2.*t4.*t10.*t11.*t13.*t15-qd1_2.*t2.*t8.*t9.*t11.*t12.*t15-qd1_2.*t4.*t7.*t9.*t10.*t13.*t14+qd2_2.*t2.*t3.*t5.*t12.*t13.*t15+qd2_3.*t2.*t3.*t8.*t10.*t11.*t15+qd2_2.*t2.*t4.*t8.*t11.*t12.*t15+qd2_3.*t2.*t5.*t8.*t9.*t12.*t15-qd2_2.*t4.*t7.*t8.*t9.*t12.*t14-qd2_3.*t3.*t4.*t7.*t12.*t13.*t14+qd2_3.*t4.*t6.*t7.*t9.*t12.*t15-qd2_2.*t2.*t9.*t10.*t11.*t13.*t15+t2.*t3.*t4.*t7.*t12.*t13.*t15.*(3.0./2.0)+qd1_2.*t2.*t4.*t7.*t8.*t9.*t10.*t15+qd2_3.*t2.*t3.*t4.*t7.*t8.*t12.*t15-qd2_2.*t2.*t4.*t7.*t9.*t12.*t13.*t15).*(1.0./2.0)+t24.*(t351+t352+t353+t354+t355).*(1.0./2.0)+t163.*(t2.*t7.*t10.*(3.0./2.0)-t2.*t4.*t11.*t12.*(3.0./2.0)+qd1_2.*t2.*t5.*t9.*t10-qd1_2.*t2.*t3.*t11.*t12+qd2_3.*t2.*t3.*t5.*t12+qd1_2.*t6.*t7.*t12.*t13-qd2_2.*t6.*t7.*t8.*t10-qd2_3.*t2.*t9.*t10.*t11+t3.*t6.*t10.*t11.*t13.*(3.0./2.0)+t5.*t6.*t9.*t12.*t13.*(3.0./2.0)+qd1_2.*t2.*t3.*t4.*t7.*t10-qd1_2.*t3.*t5.*t6.*t8.*t10+qd1_2.*t4.*t6.*t10.*t11.*t13-qd2_3.*t2.*t4.*t7.*t9.*t12-qd1_2.*t6.*t8.*t9.*t11.*t12+qd2_2.*t3.*t5.*t6.*t12.*t13+qd2_3.*t3.*t6.*t8.*t10.*t11+qd2_2.*t4.*t6.*t8.*t11.*t12+qd2_3.*t5.*t6.*t8.*t9.*t12-qd2_2.*t6.*t9.*t10.*t11.*t13+t3.*t4.*t6.*t7.*t12.*t13.*(3.0./2.0)+qd1_2.*t4.*t6.*t7.*t8.*t9.*t10+qd2_3.*t3.*t4.*t6.*t7.*t8.*t12-qd2_2.*t4.*t6.*t7.*t9.*t12.*t13).*(1.0./2.0)-t26.*(-t188+t295+t296).*(1.0./2.0)+t56.*(t3.*t5.*t8.*t12.*t14.*(-3.0./2.0)+t8.*t9.*t10.*t11.*t14.*(3.0./2.0)+qd1_2.*t3.*t5.*t6.*t10.*t15+qd1_2.*t5.*t9.*t10.*t13.*t14-qd1_2.*t3.*t11.*t12.*t13.*t14+qd1_2.*t6.*t9.*t11.*t12.*t15+qd2_2.*t3.*t8.*t10.*t11.*t14-qd2_3.*t3.*t6.*t10.*t11.*t15+qd2_2.*t5.*t8.*t9.*t12.*t14+qd2_3.*t3.*t5.*t12.*t13.*t14-qd2_3.*t5.*t6.*t9.*t12.*t15-qd2_3.*t9.*t10.*t11.*t13.*t14-t2.*t3.*t5.*t12.*t13.*t15.*(3.0./2.0)+t4.*t7.*t8.*t9.*t12.*t14.*(3.0./2.0)+t2.*t9.*t10.*t11.*t13.*t15.*(3.0./2.0)-qd1_2.*t2.*t5.*t8.*t9.*t10.*t15+qd1_2.*t2.*t3.*t8.*t11.*t12.*t15+qd1_2.*t3.*t4.*t7.*t10.*t13.*t14-qd1_2.*t4.*t6.*t7.*t9.*t10.*t15-qd2_3.*t2.*t3.*t5.*t8.*t12.*t15+qd2_2.*t3.*t4.*t7.*t8.*t12.*t14-qd2_3.*t3.*t4.*t6.*t7.*t12.*t15+qd2_2.*t2.*t3.*t10.*t11.*t13.*t15+qd2_2.*t2.*t5.*t9.*t12.*t13.*t15+qd2_3.*t2.*t8.*t9.*t10.*t11.*t15-qd2_3.*t4.*t7.*t9.*t12.*t13.*t14+t2.*t4.*t7.*t9.*t12.*t13.*t15.*(3.0./2.0)-qd1_2.*t2.*t3.*t4.*t7.*t8.*t10.*t15+qd2_2.*t2.*t3.*t4.*t7.*t12.*t13.*t15+qd2_3.*t2.*t4.*t7.*t8.*t9.*t12.*t15).*(1.0./2.0)-t147.*(t3.*t5.*t8.*t12.*t15.*(-3.0./2.0)+t8.*t9.*t10.*t11.*t15.*(3.0./2.0)-qd1_2.*t3.*t5.*t6.*t10.*t14+qd1_2.*t5.*t9.*t10.*t13.*t15-qd1_2.*t6.*t9.*t11.*t12.*t14-qd1_2.*t3.*t11.*t12.*t13.*t15+qd2_3.*t3.*t6.*t10.*t11.*t14+qd2_2.*t3.*t8.*t10.*t11.*t15+qd2_3.*t5.*t6.*t9.*t12.*t14+qd2_2.*t5.*t8.*t9.*t12.*t15+qd2_3.*t3.*t5.*t12.*t13.*t15-qd2_3.*t9.*t10.*t11.*t13.*t15+t2.*t3.*t5.*t12.*t13.*t14.*(3.0./2.0)+t4.*t7.*t8.*t9.*t12.*t15.*(3.0./2.0)-t2.*t9.*t10.*t11.*t13.*t14.*(3.0./2.0)+qd1_2.*t2.*t5.*t8.*t9.*t10.*t14-qd1_2.*t2.*t3.*t8.*t11.*t12.*t14+qd1_2.*t4.*t6.*t7.*t9.*t10.*t14+qd1_2.*t3.*t4.*t7.*t10.*t13.*t15+qd2_3.*t2.*t3.*t5.*t8.*t12.*t14+qd2_3.*t3.*t4.*t6.*t7.*t12.*t14+qd2_2.*t3.*t4.*t7.*t8.*t12.*t15-qd2_2.*t2.*t3.*t10.*t11.*t13.*t14-qd2_2.*t2.*t5.*t9.*t12.*t13.*t14-qd2_3.*t2.*t8.*t9.*t10.*t11.*t14-qd2_3.*t4.*t7.*t9.*t12.*t13.*t15-t2.*t4.*t7.*t9.*t12.*t13.*t14.*(3.0./2.0)+qd1_2.*t2.*t3.*t4.*t7.*t8.*t10.*t14-qd2_2.*t2.*t3.*t4.*t7.*t12.*t13.*t14-qd2_3.*t2.*t4.*t7.*t8.*t9.*t12.*t14).*(1.0./2.0)+qd1_2.*t170.*t423.*(1.0./2.0)-qd1_2.*t189.*t426.*(1.0./2.0)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[6,6]);
