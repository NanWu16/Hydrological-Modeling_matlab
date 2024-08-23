clc;
clearvars -except CLM run_time

global run_time
run_time = 0;

addpath(genpath('../read'))
addpath(genpath('../write'))
addpath(genpath('../initialize'))

CON = read_conf();


% CLM = read_climatic(CON);
disp(strcat("Calibration period: ", num2str(CON.IYB), "��", num2str(CON.IMB), "��", num2str(CON.IDB), "�� --- ", num2str(CON.IYE), "��", num2str(CON.IME), "��", num2str(CON.IDE), "��")) % ��ӡģ�����е���ֹ������


[CON.IY1, CON.IM1, CON.ID1, ~, CON.IY2, CON.IM2, CON.ID2, ~] = data_length(CON.IYB, CON.IMB, CON.IDB, 0, CON.IYE, CON.IME, CON.IDE, 23);
CON.ntime = int32(length(CON.IY1));
CON.ntime1 = int32(length(CON.IY2));

RUN = read_run(CON);
disp(strcat("Total length of observed discharge:", num2str(CON.ntime1)))


UL = read_underlying(CON);


disp('Read lumped parameters')
[PARAM, param_bound] = read_params(CON);

n_iterations    = 20;
n_complexes     = 20;
n_parameters    = size(param_bound, 1);
[best_params, calibration_results] = sceua_N(param_bound, n_iterations, n_complexes, n_parameters, CLM, RUN, UL, CON);

y_sim = GXAJ_Snow17_Frozensoil_func(CLM, best_params, UL, CON);
nash_sutcliffe_efficiency(RUN,y_sim)
plot(1:length(RUN),RUN,1:length(RUN),y_sim)