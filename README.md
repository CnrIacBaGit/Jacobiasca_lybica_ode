Angela Monti Institute for Applied Mathematics (IAC), CNR, Bari, Italy. Mail: angela.monti@cnr.it

ENEI MATLAB REPOSITORY -- FILE MAP AND RUN INSTRUCTIONS
=======================================================

MAIN FILE
---------
script_solve_ode_Arancio.m
    This is the principal manuscript driver. It loads the forcing archive and
    authoritative parameter CSV, integrates the population model through the
    shared enei_rhs.m routine, and generates the main figures and tables.

MASTER LAUNCHER
---------------
RUN_ALL_ENEI.m
    Runs the complete repository in the intended order. From MATLAB, open the
    MATLAB_code folder and run:

        RUN_ALL_ENEI

    The first operation is verify_repository(), which stops if a required
    project file, MAT variable, CSV column, parameter, or forcing value is bad.

OTHER ENTRY POINTS
------------------
run_enei_analysis.m
    Generates endpoint elasticities, interaction checks, and Protocol A/B
    finite-range robustness outputs.

make_figure_R22.m
    Generates figure_R22_stage_trajectories.png (Figure 5 in the revision).

script_sensitivity_aintro.m
script_sensitivity_temperature.m
script_sensitivity_humidity.m
script_sensitivity_radiation.m
script_sensitivity_wind.m
    Generate the five Appendix-A time-dependent sensitivity figures.

AUTHORITATIVE INPUT FILES
-------------------------
ENEI_parameter_intervals_updated.csv
    The only parameter source used by the current pipeline.

Arancio_dati_2000_2025.mat
    Daily Tvec, Hvec, Rvec and Wvec forcing. The NASA POWER -999 fill value
    has been removed from the published archive.

CORE MODEL FILES
----------------
enei_rhs.m
    Shared right-hand side used by the principal simulation and auxiliary
    analyses.

simulate_enei_model.m
    Single-year Euler integrator used by elasticity/robustness analyses.

enei_sensitivity.m
    Variational equations for Appendix-A sensitivities.



The routine has been implemented and developed by Angela Monti and Fasma Diele. It can be used under the conditions of CC-BY-NC 2.0

A full description of the model is available in: Deborah Lacitignola, Angela Monti, Emanuele Gosamo, Giovanni Attolico, Massimiliano Nitti,
Fasma Diele, Carmela Marangi, Early Nymph Emergence Index (ENEI) for Jacobiasca lybica: An environmentally driven early-warning tool for vineyard risk assessment (submitted, 2026)
