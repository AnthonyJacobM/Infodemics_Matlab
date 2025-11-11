# Infodemics_Matlab
Infodemics with replicable results using Matlab Programming:

Infodemics_ODE_Sys_V2.m - solve and plot the trajectories of a 5D system with given parameters.
Infodemids_ode_Par_Sweep.m - solve and plot the trajectories of the infodemics system for a given parameter set, but 
iterate a single parameter linearly spaced vector.
matcont_data_extractor.m - loads a matcont continuation file and extracts key information such as state variables and parameters.
combine_matcont_curves.m - combines two curves (structured) from the matcont_data_extractor files.
filter_specific_labels.m - modifies the labels of the existing formats of the labels from the matcont_data_extractor.m file.
plot_nullclines.m - displays the nullclines and phase trajectories corresponding to the time trajectories listed in Infodemics_ODE_Sys.m. 
Figure_Modifications.m - modify existing opened figures to publishable format. 
