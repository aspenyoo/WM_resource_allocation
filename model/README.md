# WM resource allocation
This folder contains all the functions necessary to estimate model parameters, given data. 

## File descriptions
- **calc_E_EU.m**: calculate the expected cost of the entire experiment, given parameters, for the Maximizing Points model. 
- **calc_EU.m**: calculate the expected utility for a given circle radius, memory precision, and risk preference. for the Maximizing Points model.
- **calc_expectederror_analytical.m**: analytically calculates expected cost of the entire experiment, given parameters, for the Minimizing Error model. 
- **calc_expectederror.m**: calculates expected cost of the entire experiment using sampling, given parameters, for the Minimizing Error model. 
- **calc_nLL.m**: calculates negative log-likelihood of parameters given data and model. 
- **calc_p_Hit.m**: calculates the probability that the saccade landing is within the circle wager. 
- **calc_pdf_r.m**: calculates the pdf of choosing r for a given beta (wager noise parameter) and J (memory precision).
- **calc_pVec_maxpoints.m**: calculates the proportion to allocate to each target that maximizes points.
- **calc_pVec_minerror.m**: calculates the proportion to allocate to each target that minimizes loss
- **check_nonbcon.m**: contains non-bound constrains to the model.
- **fit_parameters.m**: estimate parameters given model and data.
- **gensimtheta.m**: generate data from particular parameter combination.
- **loadconstraints.m**: contains optimization contraints for each model. is loaded when estimating parameters in fit_parameters.m.
- **loadvar.m**: contains variables used for calculating things using sampling. 


