function T = Force_mat(in1)
%FORCE_MAT
%    T = FORCE_MAT(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    07-May-2020 17:35:11

Tt1 = in1(1,:);
Tt2 = in1(2,:);
T = reshape([Tt1,0.0,Tt2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[4,6]);