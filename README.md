# MM-diagnosis-model-construction-code
Based on metagenomic sequencing of the fecal samples of newly diagnosed patients with multiple myeloma and healthy controls, the MM-characterizing microorganisms were selected as biomarkers for constructing MM diagnosis model

In this study, an in-house metagenomics sequencing dataset from a cohort of newly diagnosed group of MM patients and healthy controls were first analyzed. The relatively convergent and abundant gut micrbiota composition was discovered in MM. Subsequently, the differential species in MM were identified by DESeq2, and eleven were further verified using qPCR. Furthermore, taking these MM-characterizing species as features, 10 diagnostic models were constructed using distinct machine learning methods, of which 9 were validated with AUC of >0.65 in the independent validation dataset. In particular, when methods rf, adaboost, and nnet were utilized, the robust predictive power was all achieved (AUC >0.8). 

This study provides a novel insight for MM diagnosis. We first characterized gut microbiome as biomarkers and showed the successful establishment of diagnostic models, indicating the potential of gut microbiome analysis towards targeted non-invasive biomarkers for early MM diagnosis. 


