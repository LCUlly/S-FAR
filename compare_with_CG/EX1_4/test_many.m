clc
clear
% 
load("data_EX32.mat")
[x_199_4,num_199_4,Err_199_4] = FAR199(0.05);
[x_15_4,num_15_4,Err_15_4] = FAR15(0.05);
[x_land_4,num_land_4,Err_land_4] = Landweber(0.05);
[x_09_4,num_09_4,Err_09_4] = FAR09(0.05);
[x_CG_4,num_CG_4,Err_CG_4] = CG(0.05);
% 
[x_199_3,num_199_3,Err_199_3] = FAR199(0.01);
[x_15_3,num_15_3,Err_15_3] = FAR15(0.01);
[x_land_3,num_land_3,Err_land_3] = Landweber(0.01);
[x_09_3,num_09_3,Err_09_3] = FAR09(0.01);
[x_CG_3,num_CG_3,Err_CG_3] = CG(0.01);
% 
[x_199_2,num_199_2,Err_199_2] = FAR199(0.001);
[x_15_2,num_15_2,Err_15_2] = FAR15(0.001);
[x_land_2,num_land_2,Err_land_2] = Landweber(0.001);
[x_09_2,num_09_2,Err_09_2] = FAR09(0.001);
[x_CG_2,num_CG_2,Err_CG_2] = CG(0.001);
% 
% % 
[x_199_1,num_199_1,Err_199_1] = FAR199(0.0001);
[x_15_1,num_15_1,Err_15_1] = FAR15(0.0001);
[x_land_1,num_land_1,Err_land_1] = Landweber(0.0001);
[x_09_1,num_09_1,Err_09_1] = FAR09(0.0001);
[x_CG_1,num_CG_1,Err_CG_1] = CG(0.0001);
% 
save('data_EX32_compare_CG')
% 
% 
% 
% 
