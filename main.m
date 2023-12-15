clear;clc;
addpath(genpath(pwd));

source_dir = 'H:\lxy\data\1_normalization\9_9\';
dis0_dir = './Results/data0/';
if exist(dis0_dir, 'dir') == 0
    mkdir(dis0_dir)
end
disOcc_dir = './Results/dataOcc/';
if exist(disOcc_dir, 'dir') == 0
    mkdir(disOcc_dir)
end
disMLR_dir = './Results/dataMLR/';
if exist(disMLR_dir, 'dir') == 0
    mkdir(disMLR_dir)
end
disTLM_dir = './Results/dataTLM/';
if exist(disTLM_dir, 'dir') == 0
    mkdir(disTLM_dir)
end
disFinal_dir = './Results/dataFinal/';
if exist(disFinal_dir, 'dir') == 0
    mkdir(disFinal_dir)
end
dataset = 'HCInew';
name = 'sideboard';

paths.source_path = [source_dir, name];
paths.dis0_path = [dis0_dir, name];
% paths.dis0_path = "H:\lxy\data\4_data_withoutAbl\hci\9_9\data0_add\antinous.mat";
paths.disOccCost_path = [disOcc_dir, name, 'Occ'];
paths.disOccFinal_path = [disOcc_dir, name, 'Final'];
paths.disMLR_path = [disMLR_dir, name];
paths.disTLM1_path = [disTLM_dir, name, '1'];
paths.disTLMMask_path = [disTLM_dir, name, 'Mask'];
paths.disTLM_path = [disTLM_dir, name, 'Tlm'];
paths.disFinal_path = [disFinal_dir, name];
%%
tic
InitialDepth
disp([name, ' Initial Depth Estimation Done'])
OcclusionCost
disp(' ')
disp([name, ' Occlusion Cost Calc Done'])
toc
tic
OcclusionRefine
disp(' ')
disp([name, ' Occlusion Refinement Done'])
toc
tic
MLRefine
disp(' ')
disp([name, ' Multi-Line Refinement Done'])
toc
tic
TexturelessRefine
disp(' ')
disp([name, ' Texture-less Refinement Done'])
toc
tic
GlobalOptimization
disp(' ')
disp([name, ' Global Optimization Done'])
toc
disp(['Result has been saved to ', paths.disFinal_path])