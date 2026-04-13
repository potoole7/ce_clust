# Read Through Notes 

## Introduction
[- Having some trouble with the bibliography (can't get authors to be listed chronologically when using citep and multiple authors), hope to fix soon ...]: #
[- On final page (5), mention "spatial locations", should this also have a mention of random vectors, or are we switching to this more specific language thereafter?]: #
[- Also on page 5, we say the discussion includes future work; do we really discuss future work? Apart from the fact that our method is applicable in many other settings?]: #

## Methods
[- In equation 3, we use $u$, but below in equation 4 we switch to $u_i$, is that okay?]: #
[- At end of 2.1, I say that $-1 < \alpha_{j \mid i} < 1$ is equivalent to $\chi < 1$. Shouldn't we say $0 < \chi < 1$ here? Or $0 \le \chi < 1$? Also what if $-1 < \chi < 0$?] : #
- Should I define Equation 6 as $H_{i, s}(y)$ when it first appears, to match Equation 10, or is it fine to define later?
[- At the end of 2.2, I mention that we choose $q$ with threshold stability plots; do I need to include these in supplementary materials for each location?]: #
[- Also here, should I mention that parameter estimates are stable *across sites*, and not just keep the general statement we have now?]: #
[- In Equation 7, do we want to use $H^*$ instead of $H_2$, to match later equations which compare sites $s$ and $s^*$?]: #
[- In Equation 9, is having ${\Omega^*}^{-1}$ the correct LaTeX syntax for taking the power of a variable with a superscript?]: #
- When I give bounds for the $JS^{G_{\lambda}}$ in terms of the Jeffrey's divergence, should I include the brief derivation (which follows from Nielsen, 2025) in the supplementary materials? Also, do we want to say anything more about this bound? Wanted to say previously that this means our divergence is less susceptible to blowing up than Jeffrey's, but don't think I found a specific reference ..
[- For MC estimation, do we need to say what to do if the 99th empirical percentile across all sites exceeds the maximum observations at any given site?]: #
[- In our algorithm, do we want to have the index below argmin? Also, is "observation" in line 5 the right word to use?]: #
- Where we say "We find in our simulation study (Section~\ref{sec:sim}) that this method correctly identifies $k$.", do we need to include supplementary plots for this?

## Simulation Study
- In 3.1, we say that we estimate $u_i$ as the standard Laplace quantile for $q=0.9$ and $q=0.99$. Shouldn't this just be $u$ then, since it'll be the same for all variables, having come from the standard Laplace distribution?
- In the caption for Figure 1, is it better to have "(top) $\alpha_{\mid i}$" or "$\alpha_{\mid i}$ (top)"? We currently use the first here, but the latter in Figure 7.
- Do we want to move Figure 4 to supplementary materials?
- Is there a need for the line "Having separate $\rhou{t}^s$ values for each variable pair would increase the complexity and signal-to-noise ratio of the clustering task, and likely reduce performance."
- Christian was concerned before that I didn't have enough description to accompany simulation plots (particularly Figure 4), is that better now? I included examples of where clustering performance was best and worst to aid the reader, which may be unnecessary or have a better alternative.

## Application
- Where we say "...data are obtained from ... and ERA5 Reanalysis \citep{Hersbach2020}, respectively.", should we say "ERA5 Reanalysis dataset"? Currently seems a bit abrupt.
- We say (about Figure 6) "Non-zero values are observed at almost all sites, indicating positive asymptotic dependence between precipitation and wind speed." Is it correct to say "asymptotic dependence" here? For $\chi < 0$ don't we have asymptotic independence, but could maybe say "some extremal dependence"?
- Need to refer to Figure 6 when mentioning the elevation of the Wickow Mountains, when this entire paragraph is already all about Figure 6?
- In the last line of 4.1, is there a need to reference notable outliers again?
- In 4.2, should we have CE diagnostic plots and/or threshold stability plots in our supplementary materials? Or is it okay to leave our threshold stability at least, since we show our results are reasonably robust to the choice of $q$ in Figure 10?
- Figure 7 is very large since it has 2 rows of maps, would you like to see if I could get all sub-plots on the same row or something?
- Should we include elbow plots in 4.3 in supplementary materials?
- Some of the discussion of outliers surrounding Figures 8 and 9 is rehashed from Figure 6 ($\chi$), so I should probably cut down on this where possible, or link better.
- Should we include plots of clustering solutions for $k = 2$ and $k = 4$ clusters in supplementary materials?

## Discussion
- Added paragraph on non-spatial possible uses for method, which hasn't been reviewed yet.

## Misc
- Could plots still be smaller? Might squeeze more text in.

## Supplementary Materials
- Would I need one for first submission? Probably if I'm going to be waiting 4 months for a reply.
- Christian's JCGS Bayesian clustering paper had a very "barebones" supplement (and didn't mention it in the main text as far as I can see), while Jordan's Neural Bayes Estimators paper in JCGS (by Matthew Sainsbury-Dale) had a more detailed one which was mentioned frequently in the main paper; which style should I go for?
- Personally, I don't think I'd really need detailed descriptions, just somewhere to dump extra plots.
