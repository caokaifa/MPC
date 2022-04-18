  a=[1,0;0,1;-1,0;0,-1];
  b=[2;2;2;2];
   c = [2;1];
   x = sdpvar(2,1);
   diagnostics = optimize(a*x<b,c'*x);
   value(x)