
# S2 – Explanations of Data Files

This document provides a detailed guide to the various files, columns, and analytical tools used in the metabolomics data processing pipeline **MAPS** (currently unpublished).

---

## Main Annotations

### `final-annotation-df.csv`

Comprehensive list of annotations (without intensity). This file represents the primary MAPS output and contains the highest-confidence annotation per feature, standardised to PubChem where applicable.

| ID | File / Column | Description |
|----|---------------|-------------|
| A | feature.ID | Most intense ion for a particular compound (related ions which have been matched by ion identity networking). |
| B | rt | Retention time in minutes. |
| C | mz | Mass to charge ratio (m/z). |
| D | compound.name | Highest confidence annotation as determined by MAPS – standardised to PubChem. |
| E | smiles | Computational interpretation of molecule – standardised to PubChem. |
| F | annotation.type | Annotation tool from which the annotation was derived. |
| G | confidence.level | Annotation confidence level. |
| H | confidence.score | Confidence score is not a probability (or percentage confidence), and is different depending on the annotation.type that provided it. Accepted thresholds were: GNPS cosine score > 0.7 (level 2 and 3), CANOPUS score > 0.7, CSI:FingerID score > 0.64 (10% FDR), MS2Query Tanimoto distance of 0.63. |
| I | id.prob | Identification probability computed as 1/N (Dorrestein, 2025). |
| J | CID | PubChem identifier. |
| K | HMDB.ID | Human Metabolome Database identifier. |
| L | Formula | Molecular formula for unionised species. |
| M | IUPAC | IUPAC compound name. |
| N | Monoisotopic.Mass | Monoisotopic mass of the compound. |
| O | mz.diff.ppm | m/z difference to the theoretical mass of the annotated adduct (in ppm). |
| P | feature.usi | Universal spectrum identifier for the feature. |
| Q | gnps.library.usi | Unique spectrum identifier to matched library spectrum in GNPS2. |
| R | gnps.cluster.ID | Feature-based molecular networking cluster ID. |
| S | gnps.in.silico.bile.acid.info | Only applicable to level 3 annotations derived from in silico GNPS libraries. |
| T | canopus.NPC.pathway | NPClassifier pathway derived from CANOPUS. |
| U | canopus.NPC.pathway.probability | CANOPUS pathway confidence score. |
| V | canopus.NPC.superclass | NPClassifier superclass derived from CANOPUS. |
| W | canopus.NPC.superclass.probability | CANOPUS superclass confidence score. |
| X | canopus.classyfire.subclass | Classyfire subclass derived from CANOPUS. |
| Y | canopus.classyfire.subclass.probability | Classyfire subclass confidence score. |
| Z | canopus.classyfire.specclass | Classyfire specific class derived from CANOPUS. |
| AA | canopus.classyfire.specclass.probability | Classyfire specific class confidence score. |
| AB | zodiac.formula | Molecular formula predicted by ZODIAC. |
| AC | zodiac.confidence.score | ZODIAC confidence score. |
| AD | Propagated.Feature.ID | Feature ID from which annotation was propagated. |
| AE | Propagated.Annotation.Type | Annotation source of propagated annotation. |
| AF | Propagated.Annotation.Class | NPC superclass of propagated annotation. |
| AG | Samples | Samples in which the feature is detected. |

---

## Abundances

### `samples-df.csv`

List of all samples and associated feature abundances (area under the curve). Columns shared with `final-annotation-df.csv` include feature.ID, feature.USI, compound.name, smiles, Formula, IUPAC, and Monoisotopic.Mass.

| Column | Description |
|-------|-------------|
| samples | One row per sample, regardless of feature presence. |
| area | Area under the curve for the feature in that sample. |

---

## Other Documents

- `top-10-features.csv` – Top 10 most abundant features.
- - `DATASET.ID-counts.csv` - Number of features annotated at each confidence level.
- `cytoscape-v2.csv` – Annotation override file for Cytoscape networks.
- `ms2query.csv` – Processed output from MS2Query.
- `data_annotations.csv` – MZMine-derived annotations and authentic standards.
- `DATA_iimn_gnps.mgf` – MS2 spectra for GNPS and MS2Query.
- `DATA_iimn_gnps_quant.csv` – Sample vs feature matrix (MS2 only).
- `data_sirius.mgf` – MS1 and MS2 spectra for SIRIUS.
- `ms1-and-ms2.csv` – Comprehensive MS1 and MS2 feature table.
- `mzmine-batch.mzbatch` – MZMine processing workflow.
- `canopus_formula_summary.tsv` – CANOPUS class predictions from formulae.
- `canopus_structure_summary.tsv` – CANOPUS class predictions from structures.
- `formula_identifications.tsv` – ZODIAC formula predictions.
- `structure_identifications.tsv` – CSI:FingerID structure predictions.
- `_top-100` suffix – Top 100 SIRIUS hits per feature.

---

## Tools

**MAPS** – Metabolite Annotation Propagation and Synthesis: Automated pipeline for processed untargeted metabolomics data (Michael Cowled, Metabolomics Australia, University of Melbourne).

**Propagation Tool** – Propagates level 1–2 annotations to clustered unknown analogues.

**Ion Identity Networking (IIN)** – Links adducts of the same compound in MZMine.

**MZMine** – Open-source MS data processing and QC software.

**GNPS** – Global Natural Product Social Molecular Networking platform.

**FBMN** – Feature-Based Molecular Networking.

**CANOPUS** – Compound class prediction via SIRIUS.

**CSI:FingerID** – In silico structure identification.

**ZODIAC** – Molecular formula prediction.

**MS2Query** – Spectral library search and analogue prediction.

**Cytoscape** – Network visualisation tool.
