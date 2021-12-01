# TorporEcophysiology
Resources for "Ecophysiological models describe biological limits to hibernating bat behavior" by Golas et al. This custom code was developed to compare alternative hypotheses (i.e. alternative mechanistic models) for torpor bout duration, fitting to data from temperature/humidity data loggers attached to big brown bats throughout hibernation. Data is currently under government review and will be made available in a public repository prior to publication.

Initial parameter sensitivity analysis to identify parameters to estimate is performed using the 'TorporBoutPSA_21_04_03.Rmd' file.

Data upload, organization, and Bayesian hierarchical modeling with figure outputs are located in 'TorporEnsembleCode.Rmd'.

Once output has been generated from 'TorporEnsembleCode.Rmd', we predict bat performance with and without Psedugymnoascus destructans in a range of stable microclimates using 'FittingModelOutputs.R'.
