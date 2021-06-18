Population-based computational approaches to investigate cardiac arrhythmia risk.

This package contains the Matlab codes used in our Chapter for the "Handbook of Modern Biophysics" series
to describe the use of population-based computational approaches to investigate cardiac electrophysiology
and proarrythmic mechanisms.

Here, the Morotti et al. model of human atrial myocyte (J Mol Cell Cardiol. 2016 Jul;96:63-71) is used to
create a population of model variants by perturbing the baseline values of some parameters, and to simulate
an electrophysiological protocol that enhances early afterdepolarization proclivity. Linear and logistic
regression analyses are perfomed to the population-level data to quantify how changes in model parameters
affect action potential and calcium transient properties in control conditions, and modulate the probability
of development of arrhythmogenic early afterdepolarizations.

_____________________________________________________________________________________________________
Contents:

readme.txt					this file

morotti_et_al_ham_ina_ran_main.m		loads initial conditions (from yf_ham_ina_ran_1Hz.mat or from 
						yf_ham_ina_ran_ACh0p1_1Hz.mat) and runs the simulation (main file)

morotti_et_al_ham_ina_ran_model_SA.m		excitation-contraction coupling model (ODE file)

morotti_et_al_ham_ina_ran_calcCurrents.m	supporting function for simulation output analysis

SA_00_generate_parameters.m			generation of random scaling factors for perturbing model parameters
						(saved into SA_par_matrix_1000_s0p1.mat)

SA_01_obtain_ICs_control.m			inclusion of parameter perturbations and execution of long 1-Hz pacing
						control simulations to reach steady-state condition for each model variant
						(with final conditions saved into SA_ICs_matrix_1000_s0p1.mat)

SA_02_beat_analysis.m				analysis of action potential and calcium transient properties on each model
						variant in the population (these outputs are determined with the Matlab function
						function_beat_analysis_EAD.m and saved into SA_output_matrix_control_1Hz_300s.mat)

SA_03_linear_regression_analysis.m		performs linear regression analysis (with the function PLS_nipals.m) on
						parameter scaling factors and action potential and calcium transient and
						plots the results		

SA_04_obtain_ICs_Ach0p1				generation of initial conditions used to simulate the pro-arrhythmic protocol
						(with final conditions saved into SA_ICs_matrix_1000_s0p1_Ach0p1.mat)

SA_05_EAD_protocol_analysis.m			cyclic execution of pro-arrhythmic protocol and assessment of presence/absence
						of early afterdepolarizations with function_EAD_occurrence.m
						(with results saved into SA_EAD_outputs_matrix_1000_s0p1.mat)

SA_06_logistic_regression_analysis.m		performs logistic regression analysis and plots the results

rotateXLabels.m					supporting function used for plotting the results
__________________________________________________________________________________________________________


Reference:

Morotti S & Grandi E. Population-based computational approaches to investigate cardiac arrhythmia risk.
In: Jue T. (eds) Molecular Modeling of Ion Channel and Cellular Function in the Heart.
Handbook of Modern Biophysics, vol 7. Springer, New York, NY.
