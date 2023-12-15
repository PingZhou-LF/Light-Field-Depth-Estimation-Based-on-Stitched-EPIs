%% global optimization
clearvars -except source_dir dis0_dir disOcc_dir disMLR_dir disTLM_dir disFinal_dir dataset name paths

th=0.07;
load(paths.source_path);
load(paths.disTLM_path);
csai=squeeze(LF(5,5,:,:,:));

TLM=zeros(size(csai,1),size(csai,2));
TLM((disTlm-dis1)~=0)=1;
%%
para.c=c0;
para.cocc=c_occ;
para.cTlm=cTlm;
para.csai=csai;
para.dis2=disTlm;
para.dis0=dis0;
para.disOcc=dis_occ;
para.TLM=TLM;

para.alpha=0.0001;
para.occEdgeEnforce= 5;
para.lcwEdgeEnforce= 2;
para.Th_minEstVar = 0.5;
para.Th_lcwEdgEnforce = 0.1;

[disFinal,cwn]= GlobalOptimize1(para);
%%    save data
save(paths.disFinal_path,'c0','cTlm','c_occ','cwn','dis0','dis1','disFinal','disTlm','dis_occ','para');