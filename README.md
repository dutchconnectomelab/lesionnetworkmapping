# Lesion Network Mapping (LNM) – Example Scripts

This repository contains example scripts for **Lesion Network Mapping (LNM)** and **symptom Lesion Network Mapping (sLNM)** analyses using voxel-wise Lead-DBS workflows and atlas-based connectivity matrices.
 
--

## Notes

- Scripts were developed and tested on **macOS (iMac, x86_64)**  
- MATLAB version: **R2024b**

---

## Getting Started

1. Download or clone this repository.
2. Set the variable `projectDir` at the beginning of each example script to your local project directory  
   (default value is `projectDir`).

---

## Optional: Running Lead-DBS

If you want to run the Lead-DBS steps included in some examples:

1. Download **Lead-DBS** https://www.lead-dbs.org and place it in: utils_leaddbs/


3. Download the **preprocessed full GSP1000 dataset** and place it in: 
   utils_leaddbs/leaddbs/connectomes/fMRI/vol

4. Download **SPM12** (required for Lead-DBS) and place it in:
   utils_leaddbs/spm

5. Open `setup_project.m` and set:
- `leadDBS_path`
- `spm_path`

---

## Running the Examples

### 1. `Example_LNM_step123.m`
- Performs individual (Steps 1 & 2) and group analysis (Step 3) for voxel-wise LNM  
- Based on Horn et al. (2015): https://www.lead-dbs.org/download-lead-dbs/
- Demonstrates the linear formulation: LNM = sum(M x C)

- Shows equivalence of two LNM implementations
- Includes multiple example lesion sets

---

### 2. `Example_LNM_randomspinLesions.m`
- Computes LNM overlap for a lesion set
- Demonstrates LNM map quality using randomized lesions via a spin model
- Includes several example lesion sets

---

### 3. `Example_LNM_specificity_randomLesions.m`
- Compares disorder-specific LNM maps to maps generated from random lesions
- Includes several example lesion sets

---

### 4. `Example_sLNM_step123.m`
- Demonstrates the 3-step **symptom LNM (sLNM)** procedure  
- Steps 1 & 2 identical to Example 1  
- Step 3 correlates lesion connectivity maps with symptom scores `sv`: sLNM = sv x (M x C)

- Includes 4 example lesion sets

---

### 5. `Example_sLNM_randomSv_randomLesions.m`
- Demonstrates that sLNM-like maps can be generated using:
- random lesions
- random symptom vectors
- Includes several published sLNM maps

---

### 6. `Example_regressionModel.m`
- Shows that most variance in LNM and sLNM maps can be explained by intrinsic
properties of the connectivity matrix `C`
- Includes published LNM and sLNM maps for validation

---

Data and code sources
---
Data sources for the examples scripts are listed in DATA_RESOURCES.txt

The functions in ./utils/toolbox_functions/, ./utils/toolbox_functions_windows/, and ./utils_leaddbs/ reuse code from the following toolboxes: Lead-DBS, Brain Connectivity Toolbox (BCT), BrainSpace, FreeSurfer MATLAB utilities, and SPM12. Please ensure that you have obtained the necessary licenses for these toolboxes before use. 

Lead-DBS — Licensed under the GNU General Public License (GPL v3). Redistribution or integration of its code requires compliance with the terms of the GPL.
Recommended citation: Horn A, Kühn A (2014). “Lead-DBS: A toolbox for deep brain stimulation electrode localizations and visualizations.” NeuroImage, doi: 10.1016/j.neuroimage.2014.12.002.

Brain Connectivity Toolbox (BCT) — Distributed under the MIT License.
Rubinov M, Sporns O (2010). “Complex network measures of brain connectivity: Uses and interpretations.” NeuroImage, doi: 10.1016/j.neuroimage.2009.10.003.

BrainSpace - Distributed under the BSD 3-Clause License.
Vos de Wael R et al. (2020). “BrainSpace: a toolbox for the analysis of macroscale gradients in neuroimaging and connectomics datasets.” Communications Biology, doi: 10.1038/s42003-020-0794-7

FreeSurfer MATLAB utilities — Subject to the FreeSurfer license agreement. Use of these functions requires acceptance of the FreeSurfer License and proper acknowledgment.
Fischl B (2012). “FreeSurfer.” NeuroImage, doi: 10.1016/j.neuroimage.2012.01.021.

SPM12 — Subject to the SPM license terms. Users must acknowledge the Wellcome Trust Centre for Neuroimaging and cite relevant SPM publications.
Penny WD, Friston KJ, Ashburner J, Kiebel SJ, Nichols TE (2007). “Statistical Parametric Mapping: The Analysis of Functional Brain Images.” Academic Press, doi: 10.1016/B978-0-12-372560-8.X5000-1

FSL templates and resources — © Oxford Centre for Functional MRI of the Brain (FMRIB). Redistribution requires acceptance of the FSL license.
Jenkinson M et al. (2012). “FSL.” NeuroImage, doi: 10.1016/j.neuroimage.2011.09.015.

FreeSurfer resources — © Massachusetts General Hospital.

Video abstract
---------
A video abstract of Lesion Network Mapping and compressed implementation can be found at:
https://youtu.be/IIISJPYg7Eo

Citation
---------

If you use any scripts, data, or resources from this repository in academic work, please cite the following paper:
van den Heuvel, M. P., Libedinsky, I., Quiroz Monnens, S., Repple, J., Sommer, I., & Cocchi, L.  Investigating the methodological foundation of lesion network mapping. Nature Neuroscience (2025). https://doi.org/10.1038/s41593-025-02196-7

