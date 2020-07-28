# NIH_NCBI

Development of Machine Learning-Based Prediction Models for Chemical Modulators of the Glucocorticoid Receptor Signaling Pathway Using Public-Domain Bioactivity Data

Objectives:  We present a model that can predict glucocorticoid receptor (GR) activity based on the structure of its small molecule chemical modulators using a machine learning approach.  The GR signaling pathway varies depending on the small molecules that bind to it, and can act as an agonist, antagonist, or not act at all.  Due to the uncertainty associated with this signaling pathway and the availability of GR high-throughput screening (qHTS) bioassay data on PubChem, the world’s largest freely accessible chemical database, an algorithm can be trained and validated to create a GR behavior model.  We can use this information to determine the chemical substructures that play the largest role in determining GR activity.  We hope this algorithm will allow for a greater understanding of GR pathway dynamics to be used in predictive analytics, intracellular modeling, and drug discovery.  We also pose this model’s development pipeline and use of open-source PubChem data as a framework to predict the behavior of additional receptors.
 
Solution Concept: We used small molecule structure to predict GR activity, building six different machine learning approaches and testing on five different machine-readable chemical structure keys, or molecular fingerprints.  We conducted statistical analysis on the qHTS data to determine the most activity-significant chemical substructures.
 
Measurements and Main Results: Six machine learning approaches, Naïve Bayes (NB), Decision Trees (DT), Random Forest (RF), K Nearest Neighbors (KNN), Support Vector Machine (SVM), and Neural Networks (NN), were all built using Tox21 qHTS data.  Each model took five different molecular fingerprint types as input and predicted GR activity – “Active” or “Inactive” – on a test set and two external datasets. These predictions had an associated Area under the Curve (AUC), Balanced Accuracy Score (BACC), Sensitivity, and Specificity.  While the RF, KNN, SVM, and NN all had a test set AUC of 0.96, the RF model performed the most consistently across all fingerprint types, with an AUC range of 0.86% – 0.96%.
 
Conclusions: Machine learning models built using PubChem open-source bioassay data are a viable approach to predicting GR receptor behavior.  It is necessary to train the models on a greater number of compounds, to increase the general applicability and external dataset performance of the model predictions. 


File Locations:
Eqv All Jupyter Notebook Scripts: Github/NCBI repo
All linux machine scripts, input files, output files etc: Shreya Scripts Zipped File.zip
All presentations/docs/papers: Shreya.zip
All dr kim’s scripts and tox21 filder: All Files.zip
