%%
% Create a matrix where we will store the coordinates. It has 9 rows
% because there are 3 replicate images x 3 mapping iterations
XYZ_Mqn_3DMD_7160LM = zeros(9,21480);

%%
%Change directory to first mapping iteration
cd 'C:\Users\Julie White\Box\PSU_KU_Collaboration\Vectra_vs_3DMD\ImagesAndRegistration\3DMD\WithPowder\Mapping_1\IMAGES\41 MAT MAPPED';

%Read in each file and assign to the first three rows.
files = dir('*.mat');
for i = 1:length(files)
    tmp = load(files(i).name);
    A = tmp.in;
    XYZ_Mqn_3DMD_7160LM_2(i,:) = reshape(A.Vertices.',1,21480);
end

%%
%Change directory to second mapping iteration
cd 'C:\Users\Julie White\Box\PSU_KU_Collaboration\Vectra_vs_3DMD\ImagesAndRegistration\3DMD\WithPowder\Mapping_2\IMAGES\41 MAT MAPPED';

%Read in each file and assign to the second three rows.
files = dir('*.mat');
for i = 1:length(files)
    tmp = load(files(i).name);
    A = tmp.in;
    XYZ_Mqn_3DMD_7160LM_2(i+3,:) = reshape(A.Vertices.',1,21480);
end

%%
%Change directory to third mapping iteration
cd 'C:\Users\Julie White\Box\PSU_KU_Collaboration\Vectra_vs_3DMD\ImagesAndRegistration\3DMD\WithPowder\Mapping_3\IMAGES\41 MAT MAPPED';

%Read in each file and assign to the third three rows.
files = dir('*.mat');
for i = 1:length(files)
    tmp = load(files(i).name);
    A = tmp.in;
    XYZ_Mqn_3DMD_7160LM_2(i+6,:) = reshape(A.Vertices.',1,21480);
end

%%
%Change directory to 3DMD folder
cd 'C:\Users\Julie White\Box\PSU_KU_Collaboration\Vectra_vs_3DMD\ImagesAndRegistration\3DMD';

%Write matrix to file 
dlmwrite('XYZ_Mqn_3DMD_7160LM_2.txt',XYZ_Mqn_3DMD_7160LM_2);

%%
% Create a matrix where we will store the coordinates. It has 315 rows
% because there are 35 people x 3 replicate images x 3 mapping iterations
XYZ_35ind_3DMD_Up_7160LM = zeros(315,21480);

%%
%Change directory to first mapping iteration
cd 'C:\Users\Julie White\Box\PSU_KU_Collaboration\Vectra_vs_3dMD\ImagesAndRegistration\3DMD_UpsampledMesh\20190828\41 MAT MAPPED_RUN1';

%Read in each file and assign to the first 105 rows.
files = dir('*.mat');
for i = 1:length(files)
    tmp = load(files(i).name);
    A = tmp.in;
    XYZ_35ind_3DMD_Up_7160LM(i,:) = reshape(A.Vertices.',1,21480);
end

%%
%Change directory to second mapping iteration
cd 'C:\Users\Julie White\Box\PSU_KU_Collaboration\Vectra_vs_3DMD\ImagesAndRegistration\3DMD_UpsampledMesh\20190828\41 MAT MAPPED_RUN2';

%Read in each file and assign to the second 105 rows.
files = dir('*.mat');
for i = 1:length(files)
    tmp = load(files(i).name);
    A = tmp.in;
    XYZ_35ind_3DMD_Up_7160LM(i+105,:) = reshape(A.Vertices.',1,21480);
end

%%
%Change directory to third mapping iteration
cd 'C:\Users\Julie White\Box\PSU_KU_Collaboration\Vectra_vs_3DMD\ImagesAndRegistration\3DMD_UpsampledMesh\20190828\41 MAT MAPPED_RUN3';

%Read in each file and assign to the second 105 rows.
files = dir('*.mat');
for i = 1:length(files)
    tmp = load(files(i).name);
    A = tmp.in;
    XYZ_35ind_3DMD_Up_7160LM(i+210,:) = reshape(A.Vertices.',1,21480);
end

%%
%Change directory to 3DMD folder
cd 'C:\Users\Julie White\Box\PSU_KU_Collaboration\Vectra_vs_3DMD\ImagesAndRegistration\3DMD_UpsampledMesh';

%Write matrix to file 
dlmwrite('XYZ_35ind_3DMD_Up_7160LM.txt',XYZ_35ind_3DMD_Up_7160LM);

%% Now for the Vectra%%

% Create a matrix where we will store the coordinates. It has 9 rows
% because there are 3 replicate images x 3 mapping iterations
XYZ_Mqn_Vectra_7160LM = zeros(9,21480);

%%
%Change directory to first mapping iteration
cd 'C:\Users\Julie White\Box\PSU_KU_Collaboration\Vectra_vs_3DMD\Vectra\Mannequin_Mapping_1\IMAGES\41 MAT MAPPED';

%Read in each file and assign to the first three rows.
files = dir('*.mat');
for i = 1:length(files)
    tmp = load(files(i).name);
    A = tmp.in;
    XYZ_Mqn_Vectra_7160LM(i,:) = reshape(A.Vertices.',1,21480);
end

%%
%Change directory to second mapping iteration
cd 'C:\Users\Julie White\Box\PSU_KU_Collaboration\Vectra_vs_3DMD\Vectra\Mannequin_Mapping_2\IMAGES\41 MAT MAPPED';

%Read in each file and assign to the second three rows.
files = dir('*.mat');
for i = 1:length(files)
    tmp = load(files(i).name);
    A = tmp.in;
    XYZ_Mqn_Vectra_7160LM(i+3,:) = reshape(A.Vertices.',1,21480);
end

%%
%Change directory to third mapping iteration
cd 'C:\Users\Julie White\Box\PSU_KU_Collaboration\Vectra_vs_3DMD\Vectra\Mannequin_Mapping_3\IMAGES\41 MAT MAPPED';

%Read in each file and assign to the third three rows.
files = dir('*.mat');
for i = 1:length(files)
    tmp = load(files(i).name);
    A = tmp.in;
    XYZ_Mqn_Vectra_7160LM(i+6,:) = reshape(A.Vertices.',1,21480);
end

%%
%Change directory to Vectra folder
cd 'C:\Users\Julie White\Box\PSU_KU_Collaboration\Vectra_vs_3DMD\Vectra';

%Write matrix to file 
dlmwrite('XYZ_Mqn_Vectra_7160LM.txt',XYZ_Mqn_Vectra_7160LM);

%%
% Create a matrix where we will store the coordinates. It has 315 rows
% because there are 35 people x 3 replicate images x 3 mapping iterations
XYZ_35ind_Vectra_7160LM = zeros(315,21480);

%%
%Change directory to first mapping iteration
cd 'C:\Users\Julie White\Box\PSU_KU_Collaboration\Vectra_vs_3DMD\Vectra\Mapping_1\IMAGES\41 MAT MAPPED';

%Read in each file and assign to the first 105 rows.
files = dir('*.mat');
for i = 1:length(files)
    tmp = load(files(i).name);
    A = tmp.in;
    XYZ_35ind_Vectra_7160LM(i,:) = reshape(A.Vertices.',1,21480);
end

%%
%Change directory to second mapping iteration
cd 'C:\Users\Julie White\Box\PSU_KU_Collaboration\Vectra_vs_3DMD\Vectra\Mapping_2\IMAGES\41 MAT MAPPED';

%Read in each file and assign to the second 105 rows.
files = dir('*.mat');
for i = 1:length(files)
    tmp = load(files(i).name);
    A = tmp.in;
    XYZ_35ind_Vectra_7160LM(i+105,:) = reshape(A.Vertices.',1,21480);
end

%%
%Change directory to third mapping iteration
cd 'C:\Users\Julie White\Box\PSU_KU_Collaboration\Vectra_vs_3DMD\Vectra\Mapping_3\IMAGES\41 MAT MAPPED';

%Read in each file and assign to the second 105 rows.
files = dir('*.mat');
for i = 1:length(files)
    tmp = load(files(i).name);
    A = tmp.in;
    XYZ_35ind_Vectra_7160LM(i+210,:) = reshape(A.Vertices.',1,21480);
end

%%
%Change directory to Vectra folder
cd 'C:\Users\Julie White\Box\PSU_KU_Collaboration\Vectra_vs_3DMD\Vectra';

%Write matrix to file 
dlmwrite('XYZ_35ind_Vectra_7160LM.txt',XYZ_35ind_Vectra_7160LM);
