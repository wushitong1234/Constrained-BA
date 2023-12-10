clear; clc;

number=1; tol=1e-10;
opts_lap.parp=1; opts_lap.lambda=1;
opts_gas.parp=2; opts_gas.sigma=1;
for i = 1:5
    D = 0.1+(i-1)*0.2;
    [t_gas(i,1),t_gas(i,2),t_gas(i,3),e_gas(i,1),e_gas(i,2),e_gas(i,3),n_gas(i,1),n_gas(i,2),n_gas(i,3),n_gas(i,4),R_gas(i,1),R_gas(i,2)] = RD_CBA(D,number,tol,opts_gas);
    [t_lap(i,1),t_lap(i,2),t_lap(i,3),e_lap(i,1),e_lap(i,2),e_lap(i,3),n_lap(i,1),n_lap(i,2),n_lap(i,3),n_lap(i,4),R_lap(i,1),R_lap(i,2)] = RD_CBA(D,number,tol,opts_lap);
end
save ('data_R_10_100.mat','t_gas','t_lap','e_gas','e_lap','n_gas','n_lap','R_gas','R_lap')