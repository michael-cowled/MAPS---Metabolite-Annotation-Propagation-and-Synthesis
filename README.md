# MAPS: Metabolome Annotation Propagation and Synthesis
## Data File Explanations and Pipeline Documentation

This repository contains the outputs and documentation for the **MAPS** metabolomics data processing pipeline[cite: 1, 19]. [cite_start]MAPS is an automated pipeline designed for processed untargeted metabolomics data[cite: 1, 19].

---

## 1. Main Output Files

### final-annotation-df.csv
[cite_start]The primary output containing a comprehensive list of annotations for all chemical features [cite: 4, 34-1].

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
[cite_start]Lists all data files (samples) and their associated abundances for each feature [cite: 11, 34-2].
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
| **MS2Query** | Tanimoto | **0.63** | Recommended distance for reliable analogues. |


---

## 3. MAPS Propagation Tool

The propagation tool attempts to annotate previously unknown features after initial high-confidence assignments.

1. **Clustering**: Groups unknown features with known compounds using the `gnps.cluster.id`.
2. **Matching**: Identifies the highest similarity feature within a cluster by sorting by cosine score.
3. **Transfer**: If the match has Level 1 or 2 confidence, the identity is propagated to the unknown feature as a purported analogue.

---

## 4. Confidence Levels

Annotations are assigned levels based on the following criteria:

* [cite_start]**Level 1**: Authentic Standard (requires MS/MS and retention time match) [cite: 27, 34-6].
* [cite_start]**Level 2**: MS/MS spectral library match[cite: 28].
* [cite_start]**Level 3**: Chemical analogue or *in silico* match[cite: 29].
* [cite_start]**Level 4**: Compound Class identification[cite: 30].
* [cite_start]**Level 5**: Molecular Formula or HRMS[cite: 31].

---

## 5. Supplemental Data Files

| File | Description |
| :--- | :--- |
| **top-10-features.csv** | The 10 most abundant features from the sample dataset. |
| **data_sirius.mgf** | [cite_start]MS1/MS2 info for formula and class identification in SIRIUS [cite: 12, 34-9]. |
| **data_annotations.csv** | [cite_start]Annotations from MZmine library and lipid tools [cite: 13, 34-6]. |
| **formula_identifications.tsv** | [cite_start]Raw formula prediction data from the ZODIAC tool [cite: 16, 34-14]. |
| **mzmine-batch.mzbatch** | The sequence of processing parameters used in MZmine. |
