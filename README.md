# MPC

关于MPC解可以参考资料：
1.上面给出的MPC.pdf文档
2.关于OSQP求解器：https://osqp.org/docs/examples/mpc.html#matlab，这个python采用的 Sparse QP Formulation
3.https://zhuanlan.zhihu.com/p/141871796  Condensed QP Formulation

关于Hierarchical-MPC的工程是基于ROS开发的，这个仿真可以使用f110_simulator工程下面的仿真环境，使用如下:

1.把这个Hierarchical-MPC与f110_simulator放到工作空间编译
2.启动仿真环境：roslaunch f110_simulator simulator.launch
3.启动控制算法：roslaunch hmpc_auto_race hmpc.launch
