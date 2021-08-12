# Multi-view-clustering-via-deep-concept-factorization
Code for the clustering model: MCDCF.

Run 'performMultiViewClusteringVarDeepConceptFactorization.m' in Matlab to perform MCDCF.

Published paper on Knowledge-Based Systems, please cite the paper if you use this code for any research purpose. Main components were contributed by Shuai Chang(Aero Asukara), All Rrights Reserved, 2019-2020.
Cite info:

@article{CHANG2021106807,
title = {Multi-view clustering via deep concept factorization},
journal = {Knowledge-Based Systems},
volume = {217},
pages = {106807},
year = {2021},
issn = {0950-7051},
doi = {https://doi.org/10.1016/j.knosys.2021.106807},
url = {https://www.sciencedirect.com/science/article/pii/S0950705121000708},
author = {Shuai Chang and Jie Hu and Tianrui Li and Hao Wang and Bo Peng},
keywords = {Multi-view clustering, Deep clustering, Concept factorization, Kernel method},
abstract = {Recent studies have shown the satisfactory results of the matrix factorization technique in Multi-view Clustering (MVC). Compared with the single-layer formed clustering models, the deep matrix factorization clustering models can better perceive the hierarchical information of the data, thereby increasing the clustering performance. Nowadays, a particular matrix factorization technique called Concept Factorization (CF) has got extensive attention in clustering research. Nevertheless, state-of-the-art CF-based clustering methods fail to simultaneously integrate MVC and deep CF into a unified framework. In this paper, a novel MVC model is proposed to tackle this challenge, which brings deep CF to MVC for learning the hierarchical information through performing multi-layer CF and derives a common consensus representation matrix to fetch the shared features among different views. Besides, we impose manifold regularization on each view for locally geometrical structure retention and employ Gaussian kernel to map the original space to a higher Hilbert space for effectively distinguishing the data points. Finally, an efficient optimization algorithm with theoretically guaranteed convergence is developed to solve the proposed model. Experiment results on several open datasets demonstrate the superior performance of the proposed model compared with baseline methods.}
}
