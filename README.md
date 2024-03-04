# Description of the project
This work reduces the pilot overhead by considering the channel subspace. It exploits the existing correlation among RIS elements that shrinks the channel space dimension. It resulted in an accepted paper at the ICC 2023 workshop. 

IEEE
- [M. Haghshenas, P. Ramezani, M. Magarini and E. Björnson, "A New Channel Subspace Characterization for Channel Estimation in RIS-Aided Communications," 2023 IEEE International Conference on Communications Workshops (ICC Workshops), Rome, Italy, 2023, pp. 1523-1528, doi: 10.1109/ICCWorkshops57953.2023.10283646. keywords: {Surface waves;Transmitters;Conferences;Millimeter wave technology;Channel estimation;Line-of-sight propagation;Receivers;Millimeter wave communication;Task analysis;MIMO communication;Holographic MIMO;reconfigurable intelligent surface;channel estimation;channel subspace characterization},](https://doi.org/10.1109/ICCWorkshops57953.2023.10283646)

arXiv (Open Access)
- [Haghshenas, Mehdi, Parisa Ramezani, Maurizio Magarini, and Emil Björnson. "A New Channel Subspace Characterization for Channel Estimation in RIS-Aided Communications." arXiv preprint arXiv:2304.02087 (2023).](https://doi.org/10.48550/arXiv.2304.02087)

# Functions
- AlgRun.m: The primary function that generates Figures 4 and 5.
- AsympRatioFigure.m: The function re-produces Figure 2.
- AzElGraph.m: Function to demonstrate orthogonal beam in Azimuth-Elevation plane as illustrated in Figure 1.
- BuildDFT.m: It generates the DFT codebook used by the Least Squares estimator.
- DFTbookSearch.m: It shows that some columns of the DFT do not correspond to any angle in space.
- UPA_BasisElupnew.m: It iteratively finds a set of azimuth and elevation angle pairs corresponding to orthogonal beams in space.
- UPA_Codebook.m: To generate the RIS configuration from a set of azimuth and elevation angles.
- UPA_Evaluate.m: The array response of a structured RIS.
- correlationUPA.m: produce Figure 3.
  
