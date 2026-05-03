# Data Directory

## ⚠️ Privacy Notice

The raw patient-level dataset used in this study **cannot be publicly shared** due to:
- IRB privacy restrictions (University of Sadat City ethical approval)
- Patient confidentiality agreements
- Declaration of Helsinki compliance requirements

## 📋 Cohort Summary (publicly available)

| Cohort | Features | Total n | Controls | HBV Patients |
|--------|----------|---------|----------|--------------|
| A | SNP genotyping (11 SNPs) | 200 | 104 | 96 |
| B | ELISA cytokines (5) | 143 | 46 | 97 |
| C | miRNA expression (12) | 98 | 48 | 50 |
| D | ELISA + miRNA (primary ML cohort) | 77 | 29 | 48 |
| E | All features (SNP + ELISA + miRNA) | 73 | 29 | 44 |

**Total enrolled:** 237 subjects (119 controls, 118 chronic HBV patients)

## 📊 Feature Description

### ELISA Cytokines (5 features)
| Feature | Description | Units | Transform |
|---------|-------------|-------|-----------|
| IL10_ELISA | Interleukin-10 serum level | pg/mL | log₁₀ |
| TGFb_ELISA | TGF-β serum level | pg/mL | log₁₀ |
| IL6_level | Interleukin-6 serum level | pg/mL | log₁₀ |
| TNFA_ELISA | TNF-α serum level | pg/mL | log₁₀ |
| IFN_ELISA | IFN-γ serum level | pg/mL | log₁₀ |

### miRNA Expression (12 features)
| Feature | miRNA | Transform |
|---------|-------|-----------|
| mir_10 | miR-10 | log₂ (2^-ΔΔCt) |
| mir_17 | miR-17 | log₂ |
| mir_21 | miR-21 | log₂ |
| mir_24 | miR-24 | log₂ |
| mir_26 | miR-26 | log₂ |
| mir_122 | miR-122 | log₂ |
| mir_125 | miR-125 | log₂ |
| mir_145 | miR-145 | log₂ |
| mir_146 | miR-146 | log₂ |
| mir_148 | miR-148 | log₂ |
| mir_155 | miR-155 | log₂ |
| mir_221 | miR-221 | log₂ |

### SNP Genotypes (11 features)
| Feature | Locus | Gene | Coding |
|---------|-------|------|--------|
| IL10_1082_Geno | -1082 A/G | IL-10 | 1=GG, 2=GA, 3=AA |
| IL10_819_Geno | -819 C/T | IL-10 | 1=CC, 2=CT, 3=TT |
| TGFb_800_Geno | -800 G/A | TGF-β | 1=GG, 2=GA, 3=AA |
| TGFb_509_Geno | -509 C/T | TGF-β | 1=CC, 2=CT, 3=TT |
| TGFb_codon10_Geno | Codon 10 | TGF-β | 1=TT, 2=TC, 3=CC |
| TGFb_codon25_Geno | Codon 25 | TGF-β | 1=GG, 2=GC, 3=CC |
| TNFA_863_Geno | -863 C/A | TNF-α | 1=CC, 2=CA, 3=AA |
| TNFA_376_Geno | -376 G/A | TNF-α | 1=GG, 2=GA, 3=AA |
| TNFA_308_Geno | -308 G/A | TNF-α | 1=GG, 2=GA, 3=AA |
| TNFA_857_Geno | -857 C/T | TNF-α | 1=CC, 2=CT, 3=TT |
| TNFA_489_Geno | -489 G/A | TNF-α | 1=GG, 2=GA, 3=AA |

## 📬 Data Access Request

To request access to the de-identified dataset for research purposes, contact:

**Prof. Roba M. Talaat**  
Molecular Biology Department, GEBRI  
University of Sadat City, Egypt  
Email: roba.talaat@gebri.usc.edu.eg
