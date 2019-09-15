%% Combine 3dMD mannequin coordinates
% Create a matrix where we will store the coordinates. It has 9 rows
% because there are 3 replicate images x 3 mapping iterations
XYZ_Mqn_3DMD_7160LM = zeros(9,21480);

%Read in each file from first mapping iteration and assign to the first three rows.
files = dir('3DMD/Mannequin_WithPowder_Mapped/Mapping_1/*.mat');
for i = 1:length(files)
    tmp = load(fullfile(files(i).folder,files(i).name));
    A = tmp.in;
    XYZ_Mqn_3DMD_7160LM(i,:) = reshape(A.Vertices.',1,21480);
end

%Read in each file from second mapping iteration and assign to the second three rows.
files = dir('3DMD/Mannequin_WithPowder_Mapped/Mapping_2/*.mat');
for i = 1:length(files)
    tmp = load(fullfile(files(i).folder,files(i).name));
    A = tmp.in;
    XYZ_Mqn_3DMD_7160LM(i+3,:) = reshape(A.Vertices.',1,21480);
end

%Read in each file from third mapping iteration and assign to the third three rows.
files = dir('3DMD/Mannequin_WithPowder_Mapped/Mapping_3/*.mat');
for i = 1:length(files)
    tmp = load(fullfile(files(i).folder,files(i).name));
    A = tmp.in;
    XYZ_Mqn_3DMD_7160LM(i+6,:) = reshape(A.Vertices.',1,21480);
end

%Write matrix to file 
dlmwrite('3DMD/Mannequin_WithPowder_Mapped/XYZ_Mqn_3DMD_7160LM.txt',XYZ_Mqn_3DMD_7160LM);

%% Combine Vectra mannequin coordinates
% Create a matrix where we will store the coordinates. It has 9 rows
% because there are 3 replicate images x 3 mapping iterations
XYZ_Mqn_Vectra_7160LM = zeros(9,21480);

%Read in each file from first mapping iteration and assign to the first three rows.
files = dir('Vectra/Mannequin_WithPowder_Mapped/Mapping_1/*.mat');
for i = 1:length(files)
    tmp = load(fullfile(files(i).folder,files(i).name));
    A = tmp.in;
    XYZ_Mqn_Vectra_7160LM(i,:) = reshape(A.Vertices.',1,21480);
end

%Read in each file from second mapping iteration and assign to the second three rows.
files = dir('Vectra/Mannequin_WithPowder_Mapped/Mapping_2/*.mat');
for i = 1:length(files)
    tmp = load(fullfile(files(i).folder,files(i).name));
    A = tmp.in;
    XYZ_Mqn_Vectra_7160LM(i+3,:) = reshape(A.Vertices.',1,21480);
end

%Read in each file from third mapping iteration and assign to the third three rows.
files = dir('Vectra/Mannequin_WithPowder_Mapped/Mapping_3/*.mat');
for i = 1:length(files)
    tmp = load(fullfile(files(i).folder,files(i).name));
    A = tmp.in;
    XYZ_Mqn_Vectra_7160LM(i+6,:) = reshape(A.Vertices.',1,21480);
end

%Write matrix to file 
dlmwrite('Vectra/Mannequin_WithPowder_Mapped/XYZ_Mqn_Vectra_7160LM.txt',XYZ_Mqn_Vectra_7160LM);

%%