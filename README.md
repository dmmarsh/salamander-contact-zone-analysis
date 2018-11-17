# salamander-contact-zone-analysis
code and data for analysis of salamander demography in contact zones between Plethodon hubrichti and P. cinereus

Uses hierarchical Bayesian models in rstan based on examples in Gelman and Hill (2006) to test for reduction in fitness-related traits
in two species of salamanders in the area where they overlap.  Files/folders indexed with "po" refer to Peaks of Otter salamanders
(Plethodon hubrichti), files/folders indexed with "rb" refer to Redback salamanders (Plethodon cinereus).  For salamander traits,
"bc" is equal to body condition, measured as the residual of the regression of mass on Snout-vent length and tail length, "tl" is
tail loss, and "hat" refers to a 1/0 variable for whether a captured salamander is a hatchling versus and adult.  

Model code adapted from: https://github.com/stan-dev/example-models/wiki/ARM-Models

Use at your own risk.  I am fairly new to rstan, but am posting these to support the concept of making data and code available as part of the 
scientific research/publication process.  Feel free to e-mail me: marshd at wlu.edu if you have any questions. 



