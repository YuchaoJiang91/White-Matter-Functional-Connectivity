% This matlab script performs K-means clustering of white matter voxels
% based on BOLD-fMRI within white matter voxels.
% 
% Revised by Yuchao Jiang, March 2021. 
% Email: yuchaojiang@fudan.edu.cn
%
% If you use this code, please cite: 
% [1] Jiang Y, et al. White-matter functional networks changes in patients with schizophrenia[J]. Neuroimage, 2019, 190: 172-181.
% [2] Jiang Y, et al. Dysfunctional white?matter networks in medicated and unmedicated benign epilepsy with centrotemporal spikes[J]. Human brain mapping, 2019, 40(10): 3113-3124.

clear;clc
%% input paramters
% the directory (func_dir) should contain a directory for each subject
fMRI_dir = '/Users/jyc_mac/MRI/WMFN_testdata/BOLD-fMRI-preprocessed'; 
WMmask_file = '/Users/jyc_mac/MRI/WMFN_testdata/WMmask.nii';
outputfiles_dir = '/Users/jyc_mac/MRI/WMFN_testdata/output';

Cluster_K = [2:3]; % Cluster number K = 2 to 5;

%% Load white matter mask
[WMmask Header_WMmask] = y_Read(WMmask_file);
WMvoxels = find(WMmask~=0); 

%% Computing the correlation matrix for clustering from each subject
% correlation matrix is a NxM matrix:
% - rows of the matrix correspond to all WM voxels
% - columns of the matrix correspond to a subsampled WM matrix (~1/4 of all WM voxels)
% - each matrix cell represents the correlation r value between the corresponding voxels BOLD-fMRI signals
% - averaged correlation matrix across subjects will be used for clustering.

% a sub-sampling of all voxels in white matter mask
subsamping = zeros(size(WMmask));    
subsamping(2:2:end, 2:2:end, 2:2:end) = 1; subsamping(1:2:end, 1:2:end, 1:2:end) = 1; % defining a grid of 1s and 0s across the image
subsamping = subsamping(WMvoxels);   % choosing only locations of white-matter voxels
subsamping = find(subsamping);    % getting the indices of subsampled voxels in the whole mask
subsamping = subsamping(randperm(length(subsamping)));     % randomly mixing the voxels' indices

% getting the data from each subject - correlation between all WM voxels and the subsampled WM voxels
FCM_for_clustering = zeros(length(WMvoxels),length(subsamping));
num_subjs_notnan = zeros(length(WMvoxels),length(subsamping));

subjects = dir(fMRI_dir); subjects = subjects(3:end); 
num_subjs = length(subjects);

for subj=1:num_subjs
    disp(['calculating correlation matrix for subject n = ', num2str(subj)])
    
    % loading current subject's functional data
    [func_data, VoxelSize, FileList, Header] = y_ReadAll(fullfile(fMRI_dir, subjects(subj).name));   % load fMRI 4D signals for each subject using y_ReadAll.m
    allvoxels_timecourses = reshape(func_data, [size(func_data,1)*size(func_data,2)*size(func_data,3), size(func_data,4)]); % 4D signals to 2D, voxels x timepoints 
    WMvoxels_timecourses = allvoxels_timecourses(WMvoxels,:);  
    clear allvoxels_timecourses; clear func_data;

    % finding voxels that have zero values or have NaN values
    WM_voxels_with_data = find(var(WMvoxels_timecourses')~=0 & ~isnan(var(WMvoxels_timecourses')));
    WM_voxels_p_with_data = find(var(WMvoxels_timecourses(subsamping,:)')~=0 & ~isnan(var(WMvoxels_timecourses(subsamping,:)')));
    
    % Calculating correlation matrix of each WM voxel to the subsampled voxels
    corr_matrix = corr(WMvoxels_timecourses(WM_voxels_with_data,:)', WMvoxels_timecourses(subsamping(WM_voxels_p_with_data),:)');
    
    % adding the current connectivity matrix to the sum of all subjects, for later averaging
    FCM_for_clustering(WM_voxels_with_data, WM_voxels_p_with_data) = FCM_for_clustering(WM_voxels_with_data, WM_voxels_p_with_data) + corr_matrix;
    num_subjs_notnan(WM_voxels_with_data, WM_voxels_p_with_data) = num_subjs_notnan(WM_voxels_with_data, WM_voxels_p_with_data) + 1;
    
    clear WMvoxels_timecourses; clear corr_matrix;
end

% calculating mean connectivity matrix from all subjects
FCM_for_clustering = FCM_for_clustering ./ num_subjs_notnan;
FCM_for_clustering(isnan(FCM_for_clustering)) = 0;  % set NaN to 0
missing_voxels = find(std(FCM_for_clustering')==0);   % finding voxels with no data
WMvoxels(missing_voxels) = []; FCM_for_clustering(missing_voxels,:)=[];

clear num_subjs_notnan;


%% K-means clustering of average connectivity matrix (i.e., FCM_for_clustering)

for K= Cluster_K
    disp(['clustering when K = ' num2str(K)]);
    IDX_allsubjs = kmeans(FCM_for_clustering, K,'distance','correlation','replicates',10);          
    clustering_results = zeros(size(WMmask)); clustering_results(WMvoxels) = IDX_allsubjs;   
    y_Write(clustering_results,Header_WMmask,[outputfiles_dir,filesep,'WM_clustering_K',num2str(K),'.nii']); % saving the results to nifit
    clear clustering_results;
end


%% Selecting the best K on the group-level data, by measuring stability of clustering solutions
% Here we separate the correlation matrix columns into 4 groups
% (cross-validation folds), therefore selecting a subset of the features.
% We perform the clustering on each fold, and then measure the similarity
% between clustering solutions for all pairs of folds. This is repeated for
% each K (number of clusters). 


CV_repeated_times = 10; % 4fold-CV will be repeated 10 times 

for i = 1:CV_repeated_times
    CV_perm{i} = randperm(size(FCM_for_clustering,2));
end

for CV_i = 1:CV_repeated_times
    temp_r = CV_perm{CV_i};
    FCM_for_clustering_CV = FCM_for_clustering(:,temp_r);
    % separating the data into subsets of features
    num_CV_folds = 4;
    IDX_folds_new = cell(1,num_CV_folds);
    for K=Cluster_K
        disp(['CV is repeated time = ' num2str(CV_i),'; K = ',num2str(K)]);
        IDX_folds = zeros(size(FCM_for_clustering_CV,1),num_CV_folds); IDX_folds_new{K} = zeros(size(IDX_folds));
        size_fold = size(FCM_for_clustering_CV,2)/num_CV_folds;
        for c=1:num_CV_folds        % going over folds (sub-matrices)
            %  disp(c)
            mat_corr_current = FCM_for_clustering_CV(:,round((c-1)*size_fold+1):round(c*size_fold));     % the sub-correlation-matrix
            IDX_folds(:,c) = kmeans(mat_corr_current, K,'distance','correlation','replicates',10);  % calculating the clustering result for this K
        end
       % computing the difference between adjacency matrices for each fold
        size_chunk = 100;   % (if voxels>1000, recommend size_chunk=100)£¬ otherwise it takes too much memory (~300 million numbers per matrix)
        num_chunks = floor(size(IDX_folds,1) / size_chunk);
            sum_diff_adjmats_folds = zeros(num_CV_folds); sum_common_connections_adjmats = zeros(num_CV_folds); sum_all_connections_adjmats = zeros(num_CV_folds);
            for ch1=1:num_chunks
                for ch2=1:num_chunks    % iterating over all adjacency matrix parts combinations
                current_clustering_adjmats = zeros(size_chunk,size_chunk,num_CV_folds);
                current_chunk1 = (ch1-1)*size_chunk+1 : ch1*size_chunk;
                current_chunk2 = (ch2-1)*size_chunk+1 : ch2*size_chunk;

                % creating the current adjacency matrix part, for all folds
                for c=1:num_CV_folds
                    for i=1:size_chunk
                        for j=1:size_chunk 
                            % In the adjacency matrix, cell (i,j) equals 1 if voxels (i,j) belong to the same cluster and 0 otherwise
                            % This allows comparison of clustering results even if the labels of the same clusters in each result are different
                            % (e.g. if the occipital cluster in solution 1 is labeled as cluster number 4, and in solution 2 it's labeled as cluster 7)
                            if IDX_folds(current_chunk1(i),c)==IDX_folds(current_chunk2(j),c)
                                current_clustering_adjmats(i,j,c)=1;
                            end
                        end
                    end
                end

                % Computing the difference between adjacency matrix for different folds, and adding this difference to the sum matrix
                for c1=1:num_CV_folds
                    for c2=1:num_CV_folds
                        sum_diff_adjmats_folds(c1,c2) = sum_diff_adjmats_folds(c1,c2) + (sum(sum(current_clustering_adjmats(:,:,c1)~=current_clustering_adjmats(:,:,c2))));
                        sum_common_connections_adjmats(c1,c2) = sum_common_connections_adjmats(c1,c2) + (sum(sum(current_clustering_adjmats(:,:,c1) & current_clustering_adjmats(:,:,c2))));
                        sum_all_connections_adjmats(c1,c2) = sum_all_connections_adjmats(c1,c2) + (sum(sum(current_clustering_adjmats(:,:,c1) + current_clustering_adjmats(:,:,c2))));                    
                    end
                end
            end
        end

        % Calculating the average difference between adjacency matrices (across all folds pairs) 
        sum_diff_adjmats_folds_all(K) = mean(sum_diff_adjmats_folds(~eye(num_CV_folds)));   % sum of all differences
        Temp_Dice_coefficient = sum_common_connections_adjmats * 2 ./ sum_all_connections_adjmats;
        Temp_Dice_coefficient_folds_all(K) = mean(Temp_Dice_coefficient(~eye(num_CV_folds)));    % Dice's coef is 1 for perfect match, 0 for no commonalities

    end
    Result.Dice_coefficient(:,CV_i) = Temp_Dice_coefficient_folds_all;
    save([outputfiles_dir,filesep,'CV_similarity.mat'],'Result');
end

Result.Dice_coefficient_mean = mean(Result.Dice_coefficient(2:20,:),2);
Result.Dice_coefficient_std = std(Result.Dice_coefficient(2:20,:),1,2);

% plotting the stability results for all K values, to identify peaks
axis_y = Result.Dice_coefficient_mean;
axis_x = (2:20)';
axis_errbar = Result.Dice_coefficient_std;
figure;errorbar(axis_x, axis_y, axis_errbar); 

save([outputfiles_dir,filesep,'CV_similarity.mat'],'Result');







