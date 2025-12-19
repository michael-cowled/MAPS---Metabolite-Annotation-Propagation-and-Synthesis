# MAPS: Metabolome Annotation Propagation and Synthesis
## Data File Explanations and Pipeline Documentation

This repository contains the outputs and documentation for the **MAPS** metabolomics data processing pipeline. MAPS is an automated pipeline designed for processed untargeted metabolomics data.

---

## 1. Main Output Files

### final-annotation-df.csv
The primary output containing a comprehensive list of annotations for all chemical features.

| Column | Description |
| :--- | :--- |
| **feature.ID** | The most intense ion for a particular compound, including related ions matched by ion identity networking. |
| **rt** | Retention time in minutes. |
| **mz** | Mass-to-charge ratio ($m/z$). |
| **compound.name** | Highest confidence annotation as determined by MAPS, standardized to PubChem. |
| **smiles** | Computational interpretation of the molecule, standardized to PubChem. |
| **confidence.score** | A non-probabilistic value that differs based on the tool that provided the annotation. |
| **confidence.level** | Ranked 1 (highest) to 5 (lowest) based on evidence. |
| **NPC.pathway** | Natural Product pathway (NPClassifier). *Note: Currently mirrored from canopus.NPC.pathway due to a bug*. |
| **NPC.superclass** | Natural Product superclass (NPClassifier). *Note: Currently mirrored from canopus.NPC.superclass due to a bug*. |
| **gnps.cluster.ID** | ID for clusters of similar chemical features (Cosine > 0.7) used for annotation propagation. |
| **feature.usi** | Universal Spectrum Identifier for the specific spectrum. |

### samples-df.csv
Lists all data files (samples) and their associated abundances for each feature.
* **samples**: Individual data files included regardless of whether a specific feature is found.
* **area**: Abundance represented as the area under the curve.

---

## 2. Technical Tools & Scoring Thresholds

MAPS integrates several third-party tools. High-confidence annotations are accepted based on the following tool-specific thresholds:

| Tool | Metric | Acceptance Threshold | Description |
| :--- | :--- | :--- | :--- |
| **GNPS** | Cosine Score | **> 0.7** | Measures spectral similarity of fragments. |
| **CSI:FingerID**| Confidence | **> 0.64** | Equates to a 10% False Discovery Rate (FDR). |
| **CANOPUS** | Prob. Score | **> 0.7** | Probability of compound class assignment. |
| **MS2Query** | Tanimoto | **> 0.63** | Recommended distance for reliable analogues. |


---

## 3. MAPS Propagation Tool

The propagation tool attempts to annotate previously unknown features after initial high-confidence assignments.

1. **Clustering**: Groups unknown features with known compounds using the `gnps.cluster.id`.
2. **Matching**: Identifies the highest similarity feature within a cluster by sorting by cosine score.
3. **Transfer**: If the match has Level 1 or 2 confidence, the identity is propagated to the unknown feature as a purported analogue.

---

## 4. Confidence Levels

Annotations are assigned levels based on the following criteria:

* **Level 1**: Authentic Standard (requires MS/MS and retention time match).
* **Level 2**: Putative MS/MS spectral library match.
* **Level 3**: Spectral analogue or *in silico* match.
* **Level 4**: Compound Class identification.
* **Level 5**: Molecular Formula or HRMS.

---

## 5. Supplemental Data Files

| File | Description |
| :--- | :--- |
| **top-10-features.csv** | The 10 most abundant features from the sample dataset. |
| **data_sirius.mgf** | MS1/MS2 info for formula and class identification in SIRIUS. |
| **data_annotations.csv** | Annotations from MZmine library and lipid tools. |
| **formula_identifications.tsv** | Raw formula prediction data from the ZODIAC tool. |
| **mzmine-batch.mzbatch** | The sequence of processing parameters used in MZmine. |
