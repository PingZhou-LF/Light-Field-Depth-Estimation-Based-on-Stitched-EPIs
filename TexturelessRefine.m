clearvars -except source_dir dis0_dir disOcc_dir disMLR_dir disTLM_dir disFinal_dir dataset name paths
load(paths.dis0_path);%info,dis0
dis0=info.dis0;
load(paths.source_path,'LF');%
csai=squeeze(LF(5,5,:,:,:));
% line detection
thresEPI=0.05;
lenThres=10;
thresLineNum=4;

sxLineDetection
tyLineDetection
clearvars -except labelSx sxLine tlmIndSx labelTy tyLine tlmIndTy LF csai dis0 flagEdgeSx flagEdgeTy nameAll name nameDataset stSize experiment ...
    source_dir dis0_dir disOcc_dir disMLR_dir disTLM_dir disFinal_dir dataset name paths
%% 
[imgMerge,imgSx,imgTy] = linMerge(LF,sxLine,tyLine);

save(paths.disTLM1_path,'labelSx','sxLine','tlmIndSx','labelTy','tyLine','tlmIndTy','flagEdgeSx','flagEdgeTy','imgMerge','imgSx','imgTy');

%%
clearvars -except source_dir dis0_dir disOcc_dir disMLR_dir disTLM_dir disFinal_dir dataset name paths

load(paths.dis0_path);%info,dis0
dis0=info.dis0;
load(paths.source_path,'LF');

load(paths.disTLM1_path);
load(paths.disMLR_path);
csai=squeeze(LF(5,5,:,:,:));

disp([name, ' Area Segmentation'])

type=3;%1:边界明显2：边界不明显;3:边界极不明显
areaSeg

%% Depth Refine
thresArea=1000;
[labelFinal,~] = unique(tlmRefine);
labelFinal=labelFinal(2:end);
% dis1=dis_occ;
disTlm=dis1;
[Ny,Nx]=size(disTlm);
labelRefine=labelFinal;
del=[];
for i=1:1:max(max(labelFinal))
    maskTemp=zeros(Ny,Nx);
    maskTemp(tlmRefine==i)=1;
    if(sum(sum(maskTemp))<=thresArea)
        del=[del;i];
    end
end
labelRefine(del,:)=[];

disp([name, ' Depth Refine'])
for i=2:1:size(labelRefine,1)
    maskTemp=zeros(Ny,Nx);
    maskTemp(tlmRefine==labelRefine(i,1))=1;
    [re_depth,~] = AreaOptimize( maskTemp,disTlm );
    disTlm=re_depth;
end

save(paths.disTLMMask_path, 'labelSx','sxLine','tlmIndSx','labelTy','tyLine', ...
    'tlmIndTy','flagEdgeSx','flagEdgeTy','imgMerge','imgSx','imgTy', ...
    'tlmMaskFinal','tlmRefine','paraTLM');
cTlm=c_occ;
cTlm(dis1~=disTlm)=0.7;

save(paths.disTLM_path,'dis0','dis_occ','dis1','c0','c_occ','disTlm','cTlm');