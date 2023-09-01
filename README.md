# Hierarchical Survival Models: <br/> Estimation, Prediction, Interpretation

This repository contains the material for the talk titled _Hierarchical Survival Models: Estimation, Prediction, Interpretation_, which was presented by Alessandro Gasparini at the [2023 Northern European Stata Conference](https://www.stata.com/meeting/northern-european23).

Specifically:

* The Stata script [`case-study.do`](https://github.com/RedDoorAnalytics/2023-nesc-hsurv/blob/main/case-study.do) includes the code to replicate the case study;

* The [`dataOvarian.dta`](https://github.com/RedDoorAnalytics/2023-nesc-hsurv/blob/main/dataOvarian.dta) is the case study dataset.
  This dataset combines the data that were collected in double-blind randomized clinical trials in advanced ovarian cancer.
  The dataset includes the following columns:
    + `patientID`, the identification number of a patient;
    + `trialID`, the center in which a patient was treated;
    + `trt`, the treatment indicator;
    + `timeS`, the candidate surrogate (progression-free survival);
    + `statusS`, censoring indicator for for progression-free survival;
    + `timeT`, the true endpoint (survival time);
    + `statusT`, censoring indicator for survival time.

## License

The content of this repository is released under a [Creative Commons Zero v1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/) license.

<a rel="license" href="http://creativecommons.org/publicdomain/zero/1.0/">
    <img src="http://i.creativecommons.org/p/zero/1.0/88x31.png" style="border-style: none;" alt="CC0" />
</a>
