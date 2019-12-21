ENM manuals: R tools
================
Marlon E. Cobos

  - [Description](#description)
  - [ENM manual sections](#enm-manual-sections)
      - [Getting data](#getting-data)
      - [Data cleaning](#data-cleaning)
      - [Delimitation of areas for model calibration
        (M)](#delimitation-of-areas-for-model-calibration-m)
      - [Variables processing](#variables-processing)
      - [Ecological niche modeling with Maxent and
        kuenm](#ecological-niche-modeling-with-maxent-and-kuenm)
      - [Post-modeling analyses](#post-modeling-analyses)
  - [References](#references)

<br>

## Description

This repository was created to store R scripts that help to perform
common procedures in Ecological Niche Modeling exercises. The scripts
are organized according to the main activities that need to be done when
creating ecological niche models (ENMs).

This repository accompanies a series of manuals published in the journal
<a href="https://journals.ku.edu/jbi" target="_blank">Biodiversity
Informatics</a>. These open access publications are available following
their specific links
<a href="https://journals.ku.edu/jbi/article/view/7600" target="_blank">Data
cleaning manual</a> and ENM/SDM manual (in process).

<br>

## ENM manual sections

### Getting data

Currently only one script is available in
<a href="https://github.com/marlonecobos/ENM_manuals/tree/master/Getting_data" target="_blank">this
section</a>. This script allows to download and retain georeferenced
records from the GBIF data
    base.

  - <a href="https://github.com/marlonecobos/ENM_manuals/blob/master/Getting_data/GBIF_data.R" target="_blank">Getting
    GBIF data for one or multiple species</a>

<br>

### Data cleaning

Scripts in this section help to clean multiple types of errors that can
be found in occurrence data sets. These errors decrease the quality of
results obtained when creating ENMs and can lead to misleading
conclusions; hence, data cleaning processes need to be done carefully
and consciously.

All available scripts for performing data cleaning (cleaning of species
occurrences) can be seen
<a href="https://github.com/marlonecobos/ENM_manuals/tree/master/Data_cleaning" target="_blank">here</a>.

List of scripts in this
    section:

  - <a href="https://github.com/marlonecobos/ENM_manuals/blob/master/Data_cleaning/Occurrences_initial_corrections.R" target="_blank">Initial
    corrections to occurrence
    data</a>
  - <a href="https://github.com/marlonecobos/ENM_manuals/blob/master/Data_cleaning/Out_continents_or_M.R" target="_blank">Correction
    of occurrences outside of continents and/or calibration
    areas</a>
  - <a href="https://github.com/marlonecobos/ENM_manuals/blob/master/Data_cleaning/Environmental_outlier_detection.R" target="_blank">Environmental
    outlier detection</a>

<br>

### Delimitation of areas for model calibration (M)

In this section scripts that allow users to create distinct hypotheses
of **M** (areas for model calibration) are stored. See all scripts
<a href="https://github.com/marlonecobos/ENM_manuals/tree/master/M_hypotheses" target="_blank">here</a>.

List of scripts in this
    section:

  - <a href="https://github.com/marlonecobos/ENM_manuals/blob/master/M_hypotheses/Construction_of_simple_Ms.R" target="_blank">Simple
    hypotheses of
    M</a>
  - <a href="https://github.com/marlonecobos/ENM_manuals/blob/master/M_hypotheses/M_from_polygon_intersection.R" target="_blank">M
    based on intersections of multiple hypotheses</a>

<br>

### Variables processing

To see all available scripts for processing environmental variables
(predictors) that can be used in ENMs go to the following
<a href="https://github.com/marlonecobos/ENM_manuals/tree/master/Variables_processing" target="_blank">link</a>.

List of scripts in this
    section:

  - <a href="https://github.com/marlonecobos/ENM_manuals/blob/master/Variables_processing/Masking_variables_with_M.R" target="_blank">Masking
    variables to areas for model
    calibration</a>
  - <a href="https://github.com/marlonecobos/ENM_manuals/blob/master/Variables_processing/Variables_correlation_evaluation.R" target="_blank">Variables
    correlation
    analysis</a>
  - <a href="https://github.com/marlonecobos/ENM_manuals/blob/master/Variables_processing/PCA_raster_and_projections.R" target="_blank">PCA
    with raster layers and
    projections</a>
  - <a href="https://github.com/marlonecobos/ENM_manuals/blob/master/Variables_processing/Variable_sets_from_all_combinations.R" target="_blank">Preparing
    variable sets from all their combinations</a>

<br>

### Ecological niche modeling with Maxent and kuenm

Scripts in this section help to automate critical steps of ecological
niche modeling. All available scripts for performing the analyses can be
seen
<a href="https://github.com/marlonecobos/ENM_manuals/tree/master/ENM_process" target="_blank">here</a>.

List of scripts in this
    section:

  - <a href="https://github.com/marlonecobos/ENM_manuals/blob/master/ENM_process/Model_calibration.R" target="_blank">Model
    calibration</a>
  - <a href="https://github.com/marlonecobos/ENM_manuals/blob/master/ENM_process/Final_models.R" target="_blank">Final
    model creation and
    projections</a>
  - <a href="https://github.com/marlonecobos/ENM_manuals/blob/master/ENM_process/Final_model_evaluation.R" target="_blank">Final
    model evaluation (optional step)</a>

<br>

### Post-modeling analyses

In this section the scripts listed below help to perform important
analyses that are not explicit part of ecological niche modeling
exercises. However, currently, this analyses are some of the common
(good) practices to be done in order to summarize results and make
interpretations more straightforward.

All available scripts from this section can be seen
<a href="https://github.com/marlonecobos/ENM_manuals/tree/master/Post_modeling" target="_blank">here</a>.

List of scripts in this section (under
    construction):

  - <a href="https://github.com/marlonecobos/ENM_manuals/blob/master/Post_modeling/Model_statistics.R" target="_blank">Summarizing
    model outputs (descriptive
    statistics)</a>
  - <a href="https://github.com/marlonecobos/ENM_manuals/blob/master/Post_modeling/Model_variability.R" target="_blank">Layers
    of model variance coming from distinct
    sources</a>
  - <a href="https://github.com/marlonecobos/ENM_manuals/blob/master/Post_modeling/Model_variance_partitioning.R" target="_blank">Hierarchichal
    partitioning of model
    variability</a>
  - <a href="https://github.com/marlonecobos/ENM_manuals/blob/master/Post_modeling/MOP.R" target="_blank">Anlysis
    of extrapolation risks (MOP)</a>

<br>

## References

Important references to be considered are listed below. These references
describe theoretical and methodological bases for many of the analyses
proposed above.

  - Anderson, R.P., Gonzalez, I., 2011. Species-specific tuning
    increases robustness to sampling bias in models of species
    distributions: An implementation with Maxent. Ecol. Modell. 222,
    2796–2811. <https://doi.org/10.1016/j.ecolmodel.2011.04.011>
  - Anderson, R.P., Lew, D., Peterson, A.T., 2003. Evaluating predictive
    models of species’ distributions: Criteria for selecting optimal
    models. Ecol. Modell. 162, 211–232.
    <https://doi.org/10.1016/S0304-3800(02)00349-6>
  - Barve, N., Barve, V., Jiménez-Valverde, A., Lira-Noriega, A., Maher,
    S.P., Peterson, A.T., Soberón, J., Villalobos, F., 2011. The crucial
    role of the accessible area in ecological niche modeling and species
    distribution modeling. Ecol. Modell. 222, 1810–1819.
    <https://doi.org/10.1016/j.ecolmodel.2011.02.011>
  - Cobos, M.E., Osorio-Olvera, L., Peterson, A.T., 2019a. Assessment
    and representation of variability in ecological niche model
    predictions. bioRxiv. <https://doi.org/10.1101/603100>
  - Cobos, M.E., Peterson, A.T., Barve, N., Osorio-Olvera, L., 2019b.
    kuenm: An R package for detailed development of ecological niche
    models using Maxent. PeerJ 7, e6281.
    <https://doi.org/10.7717/peerj.6281>
  - Cobos, M.E., Peterson, A.T., Osorio-Olvera, L., Jiménez-García, D.,
    2019c. An exhaustive analysis of heuristic methods for variable
    selection in ecological niche modeling and species distribution
    modeling. Ecol. Inform. 53, 100983.
    <https://doi.org/10.1016/j.ecoinf.2019.100983>
  - Muscarella, R., Galante, P.J., Soley-Guardia, M., Boria, R.A., Kass,
    J.M., Uriarte, M., Anderson, R.P., 2014. ENMeval: An R package for
    conducting spatially independent evaluations and estimating optimal
    model complexity for Maxent ecological niche models. Methods Ecol.
    Evol. 5, 1198–1205. <https://doi.org/10.1111/2041-210X.12261>
  - Owens, H.L., Campbell, L.P., Dornak, L.L., Saupe, E.E., Barve, N.,
    Soberón, J., Ingenloff, K., Lira-Noriega, A., Hensz, C.M., Myers,
    C.E., Peterson, A.T., 2013. Constraints on interpretation of
    ecological niche models by limited environmental ranges on
    calibration areas. Ecol. Modell. 263, 10–18.
    <https://doi.org/10.1016/j.ecolmodel.2013.04.011>
  - Peterson, A.T., Cobos, M.E., Jiménez‐García, D., 2018. Major
    challenges for correlational ecological niche model projections to
    future climate conditions. Ann. N. Y. Acad. Sci. 1429, 66–77.
    <https://doi.org/10.1111/nyas.13873>
  - Peterson, A.T., Papeş, M., Soberón, J., 2008. Rethinking receiver
    operating characteristic analysis applications in ecological niche
    modeling. Ecol. Modell. 213, 63–72.
    <https://doi.org/10.1016/j.ecolmodel.2007.11.008>
  - Peterson, A.T., Soberón, J., 2012. Species distribution modeling and
    ecological niche modeling: Getting the concepts right. Natureza &
    Conservação 10, 1–6. <https://doi.org/10.4322/natcon.2012.019>
  - Peterson, A.T., Soberón, J., Pearson, R.G., Anderson, R.P.,
    Martínez-Meyer, E., Nakamura, M., Araújo, M.B., 2011. Ecological
    Niches and Geographic Distributions. Princeton University Press,
    Princeton.
  - Peterson, A.T., Soberón, J., Sánchez-Cordero, V., 1999. Conservatism
    of ecological niches in evolutionary time. Science 285, 1265–1267.
    <https://doi.org/10.1126/science.285.5431.1265>
  - Saupe, E.E., Barve, V., Myers, C.E., Soberón, J., Barve, N., Hensz,
    C.M., Peterson, A.T., Owens, H.L., Lira-Noriega, A., 2012. Variation
    in niche and distribution model performance: The need for a priori
    assessment of key causal factors. Ecol. Modell. 237–238, 11–22.
    <https://doi.org/10.1016/j.ecolmodel.2012.04.001>
  - Soberón, J., 2007. Grinnellian and Eltonian niches and geographic
    distributions of species. Ecology Letters 10, 1115–1123.
    <https://doi.org/10.1111/j.1461-0248.2007.01107.x>
  - Soberón, J., Peterson, A.T., 2019. What is the shape of the
    fundamental Grinnellian niche? Theor Ecol.
    <https://doi.org/10.1007/s12080-019-0432-5>
  - Soberón, J., Peterson, A.T., 2005. Interpretation of models of
    fundamental ecological niches and species’ distributional areas.
    Biodiversity Informatics 2, 1–10.
    <https://doi.org/10.17161/bi.v2i0.4>
  - Warren, D.L., Seifert, S.N., 2011. Ecological niche modeling in
    Maxent: The importance of model complexity and the performance of
    model selection criteria. Ecol. Appl. 21, 335–342.
    <https://doi.org/10.1890/10-1171.1>
