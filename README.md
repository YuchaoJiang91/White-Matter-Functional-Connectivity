# White-Matter-Functional-Organization

Using BOLD-fMRI, ten large-scale white-matter functional networks were identified by a cluster analysis of voxel-based white-matter functional connectivity and classified into superficial, middle and deep layers of networks. 

![22](https://user-images.githubusercontent.com/102531632/161035646-05291ca0-2ea0-4464-b4a2-ab41ff6acaa3.jpg)
[Figure 1. Ten clusters of White-matter functional networks. 1. Occipital network; 2. Cerebellar network; 3. Anterior corona radiate network; 4. Orbitofrontal network; 5. Pre/post-central network; 6. Posterior callosum network; 7. Tempofrontal network; 8. Deep network; 9. Superior corona radiate network; 10. Superior temporal network.![image](https://user-images.githubusercontent.com/102531632/161035850-aea41c75-2e3d-4c81-bedf-dbcff084384c.png)
]

# Clustering algorithm
![1](https://user-images.githubusercontent.com/102531632/161032619-cb677dca-10f6-4be3-a09f-cdec2e148d7e.jpg)
                    [Figure 1. Clustering algorithm based on correlation matrix of white matter BOLD-fMRI signals]

Script WMFN_clustering.m performs K-means clustering of white matter voxels based on BOLD-fMRI within white matter voxels.
Ps: These scripts require Matlab and SPM (https://www.fil.ion.ucl.ac.uk/spm/) to be installed. 
Further details can be found in our papers.

# Paper
This method has been used in schizophrenia [1] and epilepsy [2].

If you use this code, please cite: 

[1] Jiang Y, et al. White-matter functional networks changes in patients with schizophrenia[J]. Neuroimage, 2019, 190: 172-181. https://doi.org/10.1016/j.neuroimage.2018.04.018

[2] Jiang Y, et al. Dysfunctional white-matter networks in medicated and unmedicated benign epilepsy with centrotemporal spikes[J]. Human brain mapping, 2019, 40(10): 3113-3124.
https://doi.org/10.1002/hbm.24584



# Acknowledgement

Original methods are provided in the paper "Evidence for functional networks within the human white matter" (Journal of Neuroscience, 2017).

Please cite above paper.

# Questions, suggestions and improvements

Email: yuchaojiang@fudan.edu.cn



