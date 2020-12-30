# LowStatisticalInterpretation
This is a simple framework to perform statistical interpretation in the case of very low background contamination.


## in main.C, there are two functions.
1. cal_p0(2000,1) : it is to calculate the p0 value under the background-only hypothesis. The first argument "2000" is the number of toy experiments. The second argument "1" is whether to use "logy" when plotting.
2. cal_upperlimit(1,1000,1): it is to produce the distributions of mu_hat (signal strength) and t_mu (test statistics) for "mu=1". The first argument is "mu=1". The second argument is the number of toy experiments (for each possible number of data events). The third argument "1" is whether to use "logy" when plotting. If we have obtained CLs+b, CLs for many mu values, we can obtain the upper limits of the signal strength at a confidence level (you can use the script "draw.py").
3. Plots from an example are put in the directorys "pic_obs_0/" and "pic_obs_lephad1". A document about this example will appear soon.
4. "input.root" is a root file storing the signal, background and data histograms of observable.


## Features.
1. Very simple codes and easy to understand and modify.
2. Currently, the systematic uncertainties are not implemented. But I will do it in the near future.

