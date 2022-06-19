H=[1 -1;-1 2];
f=[-2;-6];
Aeq=[1 1];
beq=0;
[x,fval,exitflag]=quadprog(H,f,[],[],Aeq,beq)
