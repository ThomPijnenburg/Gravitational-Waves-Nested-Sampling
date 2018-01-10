Gravitational Waves Nested Sampling
========
## Introduction 
The goal of this project is to write a nested sampling algorithm for a two-dimensional parameter space. This is placed in the context of gravitational wave models hypothesis testing where this algorithm can be applied in a 15 dimensional parameter space. 

### Gravitational Waves Theory

Gravitational Waves are ripples in spacetime that travel through the universe with the speed of light. These can be generated when an extremely strong gravitational field undergoes rapid changes. One recent measurement of these waves was done by LIGO, observing the inspiral and merger of two black holes. This process is so violent and energetic, disrupting the spacetime and producing gravitational waves powerful enough to be measured on earth many, many lightyears away. (more on [LIGO and gravitational waves](https://www.ligo.caltech.edu/page/what-are-gw))

### The experiment
An important part of science, and physics in particular, is to invent (mathematical/numerical) models for processes we observe in nature. The scientist collects a set of input parameters and combines them in a specific way in order to describe and predict the process in question.

When the problem grows in complexity, often the model grows in complexity too. Usually, one can come up with multiple hypotheses for a specific problem, and has to jump through a number of hoops in order to refine them and find the most accurate hypothesis.

In the case of the Gravitational Waves model, it turns out that there are 15 different input parameters that are necessary to describe the measured signal. Among these parameters are masses, spins, distance, sky position and phase at arrival (and more!). In other words, we have to refine our hypothesis in a 15 dimensional parameter space, which is enormous. 

### The Algorithm

The algorithm we are using here is [Nested Sampling](https://en.wikipedia.org/wiki/Nested_sampling_algorithm#cite_note-5), which is developed in order to compare models in Bayesian statistics. Let's go through some terminology:

*Prior* : Prior knowledge of the distribution of the parameter.

*Likelihood*: Probability of the evidence (data) given the parameter.

*Prior mass*: Fraction of volume having a greater likelihood.


The algorithm works as follows:

1. Drop *M* samples across parameter space, sampled from their priors. These are called *Live Points*:
    - Each has a *likelihood* associated to it.
    - Associated with volume such that likelihood is lowest at the surface.
    - Uniformly sampled in prior mass between 0 and 1
2. Discard live point with lowest likelihood (or highes prior mass)
    - Replace by new live point, sampled from the prior, which has a higher likelihood.
    - Statistically assign new value for the prior mass
3. Repeat 2.






