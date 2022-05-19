clc

A = build_lhs(xs,ys);
b = build_rhs(xs,ys,alpha);
gam = A/b;