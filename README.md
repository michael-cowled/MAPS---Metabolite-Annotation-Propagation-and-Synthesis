# S2 – Explanations of Data Files

This document provides a detailed guide to the various files, columns, and analytical tools used in the metabolomics data processing pipeline **MAPS** (currently unpublished).

---

## Main Annotations

### `final-annotation-df.csv`

Comprehensive list of annotations (without intensity). This file represents the primary MAPS output and contains the highest-confidence annotation per feature, standardised to PubChem where applicable.

| ID | File / Column    | Description                                                                                                                                                     |
| -- | ---------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| A  | feature.ID       | Most intense ion for a particular compound (related ions which have been matched by ion identity networking).                                                   |
| B  | rt               | Retention time in minutes.                                                                                                                                      |
| C  | mz               | Mass to charge ratio (m/z).                                                                                                                                     |
| D  | compound.name    | Highest confidence annotation as determined by MAPS – standardised to PubChem.                                                                                  |
| E  | smiles           | Computational interpretation of molecule – standardised to PubChem.                                                                                             |
| F  | confidence.score | Confidence score is not a probability (or percentage confidence), and is different depending on the annotation.type that provided it. Accepted thresholds were: |

* GNPS: cosine score > 0.7 (level 2 and 3 annotations)

* CANOPUS: score > 0.7

* CSI:FingerID: score > 0.64 (equating to a 10% FDR)

* MS2Query: Tanimoto distance of 0.63 (as recommended in the original publication) |
  | G | mz.diff.ppm | m/z difference to the theoretical mass of the annotated adduct (in ppm). |
  | H | gnps.shared.peaks | Number of MS/MS fragments that match in GNPS spectral libraries (minimum 6). |
  | I | library.name | Only more defined for GNPS libraries (and can infer the original submitter). |
  | J | library.quality | Only more defined for GNPS libraries. Definitions:

* GOLD – Synthetic, complete structural characterization with NMR, crystallography, or other standard methods as defined in the publication guidelines for *Journal of Natural Products*, privileged users

* SILVER – Isolated or lysate/crude, published data showing presence of molecule in the sample

* BRONZE – Any other putative, complete or partial annotation

For stringent analyses, consider accepting those of SILVER or GOLD. |
| K | NPC.pathway | Natural Product class pathway as described by NPClassifier ontologies. Derived from the annotation.type appended to the annotation and may not match canopus.NPC.pathway.

Note: Currently bugged, so appears identical to canopus.NPC.pathway. |
| L | NPC.superclass | Natural Product superclass as described by NPClassifier ontologies. Derived from the annotation.type appended to the annotation and may not match the canopus.NPC.superclass.

Note: Currently bugged, so appears identical to canopus.NPC.superclass. |
| M | gnps.library.usi | Unique spectrum identifier to matched library spectrum in GNPS2. Use the ‘spectrum resolver’ tool on GNPS2 to view and compare this to the feature.usi. |
| N | gnps.in.silico.bile.acid.info | Only applicable to level 3 annotations derived from in silico libraries matched in GNPS2. Contains more specific identifiers of the matched bile acid. |
| O | annotation.type | The annotation tool from which the ascribed annotation was derived (if present). |
| P | confidence.level | Annotation confidence level (see Confidence Levels section below). |
| Q | CID | PubChem identifier. |
| R | Formula | Molecular formula for unionised species. |
| S | IUPAC | International Union of Pure and Applied Chemistry compound name. |
| T | Monoisotopic.Mass | The sum of the masses of its constituent atoms, using the mass of the most abundant stable isotope for each element. |
| U | id.prob | When there are N possible compounds satisfying the matching conditions of an MS/MS spectrum, the identification probability is computed as 1/N (from Dorrestein, 2025). |
| V | canopus.NPC.pathway | Natural Product class pathway as described by NPClassifier ontologies derived from the CANOPUS tool in SIRIUS. May not match the NPC pathway ascribed by the annotation.type. |
| W | canopus.NPC.pathway.probability | The confidence score ascribed by CANOPUS to the NPC pathway in V (not a percentage). |
| X | canopus.NPC.superclass | Natural Product superclass as described by NPClassifier ontologies derived from the CANOPUS tool in SIRIUS. May not match the NPC superclass ascribed by the annotation.type. |
| Y | canopus.NPC.superclass.probability | The confidence score ascribed by CANOPUS to the NPC superclass in X (not a percentage). |
| Z | zodiac.formula | Molecular formula predicted by the ZODIAC tool in SIRIUS. |
| AA | zodiac.confidence.score | The confidence score ascribed by ZODIAC to the formula in Z (not a percentage). |
| AB | gnps.cluster.ID | Derived from FBMN analysis in GNPS2. A cluster of similar chemical features as determined by MS/MS spectral similarity with a cosine score > 0.7 and at least 6 matching fragments. This parameter is used to propagate annotations to unknown analogues (AD, AE, AF). |
| AC | feature.usi | Universal spectrum identifier specific for a specific spectrum or feature (irrespective of dataset). Use the ‘spectrum resolver’ tool on GNPS2 to view and compare this to the corresponding matched library.usi (if applicable). |
| AD | Propagated.Feature.ID | The feature ID of the propagated annotation onto the unknown analogue using the MAPS Propagation tool. |
| AE | Propagated.Annotation.Type | The annotation source of the propagated annotation onto the unknown analogue using the MAPS Propagation tool. |
| AF | Propagated.Annotation.Class | The NPC.superclass of the propagated annotation (if known) onto the unknown analogue using the MAPS Propagation tool. |
| AG | Samples | All data files that contain the feature of interest at the specified noise levels set in MZMine. |

---

-|---------------|-------------|
| A | feature.ID | Most intense ion for a particular compound (related ions which have been matched by ion identity networking). |
| B | rt | Retention time in minutes. |
| C | mz | Mass to charge ratio (m/z). |
| D | compound.name | Highest confidence annotation as determined by MAPS – standardised to PubChem. |
| E | smiles | Computational interpretation of molecule – standardised to PubChem. |
| F | confidence.score | Confidence score is not a probability (or percentage confidence), and is different depending on the annotation.type that provided it.

For GNPS, a cosine score of > 0.7 was accepted (level 2 and 3 annotations).

For CANOPUS, a score > 0.7 was accepted.

For CSI:FingerID, a score > 0.64 was accepted (equating to a 10% FDR).

For MS2Query, a Tanimoto distance of 0.63 was accepted (as recommended in the original publication). |
| G | mz.diff.ppm | m/z difference to the theoretical mass of the annotated adduct (in ppm). |
| H | gnps.shared.peaks | Number of MS/MS fragments that match in GNPS spectral libraries (min. 6). |
| I | library.name | Only more defined for GNPS libraries (and can infer the original submitter). |
| J | library.quality | Only more defined for GNPS libraries.

GOLD – Synthetic, complete structural characterization with NMR, crystallography or other standard methods as defined in the publication guidelines for Journal of Natural Products, privileged users.

SILVER – Isolated or lysate/crude, published data showing presence of molecule in the sample.

BRONZE – Any other putative, complete or partial annotation.

For stringent analyses, consider accepting those of SILVER or GOLD. |
| K | NPC.pathway | Natural Product class pathway as described by NPClassifier ontologies. Derived from the annotation.type appended to the annotation and may not match canopus.NPC.pathway.

Note: Currently bugged, so appears identical to canopus.NPC.pathway. |
| L | NPC.superclass | Natural Product superclass as described by NPClassifier ontologies. Derived from the annotation.type appended to the annotation and may not match the canopus.NPC.superclass.

Note: Currently bugged, so appears identical to canopus.NPC.superclass. |
| M | gnps.library.usi | Unique spectrum identifier to matched library spectrum in GNPS2. Use the ‘spectrum resolver’ tool on GNPS2 to view and compare this to the feature.usi. |
| N | gnps.in.silico.bile.acid.info | Only applicable to level 3 annotations derived from in silico libraries matched in GNPS2. Contains more specific identifiers of the matched bile acid. |
| O | annotation.type | The annotation tool from which the ascribed annotation was derived (if present). |
| P | confidence.level | See confidence levels below. |
| Q | CID | PubChem identifier. |
| R | Formula | Molecular formula for unionised species. |
| S | IUPAC | International Union of Pure and Applied Chemistry compound name. |
| T | Monoisotopic.Mass | The sum of the masses of its constituent atoms, using the mass of the most abundant stable isotope for each element. |
| U | id.prob | When there are N possible compounds satisfying the matching conditions of an MS/MS spectrum, the identification probability is computed as 1/N (from Dorrestein, 2025). |
| V | canopus.NPC.pathway | Natural Product class pathway as described by NPClassifier ontologies derived from the CANOPUS tool in SIRIUS. May not match the NPC pathway ascribed by the annotation.type. |
| W | canopus.NPC.pathway.probability | The confidence score ascribed by CANOPUS to the NPC pathway in V (not a percentage). |
| X | canopus.NPC.superclass | Natural Product superclass as described by NPClassifier ontologies derived from the CANOPUS tool in SIRIUS. May not match the NPC superclass ascribed by the annotation.type. |
| Y | canopus.NPC.superclass.probability | The confidence score ascribed by CANOPUS to the NPC superclass in X (not a percentage). |
| Z | zodiac.formula | Molecular formula predicted by the ZODIAC tool in SIRIUS. |
| AA | zodiac.confidence.score | The confidence score ascribed by ZODIAC to the formula in Z (not a percentage). |
| AB | gnps.cluster.ID | Derived from FBMN analysis in GNPS2. A cluster of similar chemical features as determined by MS/MS spectral similarity with a cosine score > 0.7 and at least 6 matching fragments. This parameter is used to propagate annotations to unknown analogues (AD, AE, AF). |
| AC | feature.usi | Universal spectrum identifier specific for a specific spectrum or feature (irrespective of dataset). Use the ‘spectrum resolver’ tool on GNPS2 to view and compare this to the corresponding matched library.usi (if applicable). |
| AD | Propagated.Feature.ID | The feature ID of the propagated annotation onto the unknown analogue using the MAPS Propagation tool. |
| AE | Propagated.Annotation.Type | The annotation source of the propagated annotation onto the unknown analogue using the MAPS Propagation tool. |
| AF | Propagated.Annotation.Class | The NPC.superclass of the propagated annotation (if known) onto the unknown analogue using the MAPS Propagation tool. |
| AG | Samples | All data files that contain the feature of interest at the specified noise levels set in MZMine. |

---

## Abundances

### `samples-df.csv`

A list of all data files (samples), related to the annotations in `final-annotation-df.csv` and their associated abundances (area). A matrix format can be found in the unprocessed files: `DATA_iimn_gnps_quant.csv` (MS2 features only) and `ms1-and-ms2.csv` found in the mzmine folder.

The following columns are identical to `final-annotation-df.csv`: feature.ID, feature.USI, compound.name, smiles, Formula, IUPAC, and Monoisotopic.Mass.

| Column  | Description                                                                                              |
| ------- | -------------------------------------------------------------------------------------------------------- |
| samples | All samples are given a separate row and are included regardless of whether a specific feature is found. |
| area    | Area under the curve for the associated feature.ID in that given sample.                                 |

---

## Other Documents

1. `top-10-features.csv`
   Same as `samples-df.csv` but limited to the top 10 most abundant features.

2. `cytoscape-v2.csv`
   After importing the associated `.graphml` file found in the GNPS folder into Cytoscape, this file can be imported to overwrite the original annotations from GNPS (as created in MAPS).

3. `ms2query.csv`
   Processed output file following use of MS2Query. See documentation relevant to the tool for explanations.

4. `data_annotations.csv`
   Derived from MZMine spectral library search and lipid annotation tools. For the purposes of MAPS, only in-house authentic standards with retention time matching are included here (and classed as level 1 annotations). Annotated lipids are designated level 3 due to the in silico nature of its algorithm.

5. `DATA_iimn_gnps.mgf`
   File containing the MS2 spectral information which is used for spectral library matching in GNPS and MS2Query.

6. `DATA_iimn_gnps_quant.csv`
   File containing a matrix of filenames (samples) vs feature ID with z = area for all features containing MS2 information. Also includes retention time and precursor m/z.

7. `data_sirius.mgf`
   File containing both the MS1 and MS2 spectral information which is used for formula identification, compound classification, and spectral library matching to in silico standards in SIRIUS.

8. `ms1-and-ms2.csv`
   File containing filenames (samples) vs feature.ID with several additional columns for all features (inclusive of MS1 and MS2). Includes the retention time range, precursor m/z, height, area, ion identity clusters, and associated adduct information. A more difficult file orientation to manipulate.

9. `mzmine-batch.mzbatch`
   The sequence of processing steps and associated parameters used in MZMine.

10. `canopus_formula_summary.tsv`
    Raw compound class information derived from predicted formulae using the CANOPUS tool in SIRIUS. This file is not used by MAPS (MAPS uses `canopus_structure_summary.tsv`).

11. `canopus_structure_summary.tsv`
    Raw compound class information derived from predicted structure using the CANOPUS tool in SIRIUS.

12. `formula_identifications.tsv`
    Raw formulae prediction information from the ZODIAC tool in SIRIUS.

13. `structure_identifications.tsv`
    Raw spectral library identifications from the CSI:FingerID tool in SIRIUS. These are denoted as level 3 annotations due to the use of in silico libraries.

14. Suffix `_top-100` for SIRIUS files
    A secondary export of summary files for which top K = 100 hits per feature are included.

---

## Tools

### MAPS

**Metabolome Annotation Propagation and Synthesis**
Automated pipeline for processed untargeted metabolomics data (created by Michael Cowled, University of Melbourne).

### Propagation

**MAPS Propagation Tool**
Following the assignment of level 1 and 2 annotations, the propagation tool attempts to annotate previously unknown features. For unknown features that cluster with known compounds by `gnps.cluster.ID`, the highest similarity feature is determined (by sorting by cosine score). If that feature is annotated with level 1 or 2 confidence, the annotation is propagated to the unknown feature as a purported analogue.

### IIN

**Ion Identity Networking**
An automated feature in MZMine which connects related adducts for the same compound (by analysis of retention time and MS1 feature shape correlation). MAPS collapses features in the same IIN to a single compound.

### MZMine

Open-source and platform-independent software used for mass spectrometry (MS) data processing, QC checking, and assignment of level 1 annotations to authentic standards (if applicable).

### GNPS

**Global Natural Product Social Molecular Networking**
Open-source tool for performing spectral similarity searching across a broad range of MS/MS libraries and performing feature networking by comparing the spectral similarity of all features to one another.

### FBMN

**Feature-Based Molecular Networking**
Tool for visually representing the similarities between different chemical features, connecting them in a network based on shared fragmentation patterns. Metadata can be overlaid such as biological origin to perform metabolomic analyses.

### CANOPUS

Predicts compound class based on fragmentation patterns as well as either molecular formulae (as determined by ZODIAC) or chemical structure (as predicted by CSI:FingerID).

### CSI:FingerID

Predicts a compound’s chemical structure by comparing its mass spectrometry fragmentation pattern with a vast database of known structures and their fragmentation patterns.

### ZODIAC

Predicts a compound’s molecular formula by analysing its mass spectrum, taking into account isotopic patterns and the quality of the spectrum.

### MS2Query

Annotates compounds by first performing spectral library searching and, if no match is found, uses a machine learning model to predict what kind of compound the molecule is based on its fragmentation patterns (i.e. determines potential analogues).

### Cytoscape

Open-source network visualisation tool.
