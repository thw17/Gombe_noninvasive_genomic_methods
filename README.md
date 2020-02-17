# Gombe noninvasive genomics
This repository contains the code used in Ozga, Webster, et al.
to investigate the promise of various sources of noninvasively
collected DNA in primate genomics.

For questions, contact: Tim Webster, Arizona State University

# Citation
If you use any code in this repository, please cite our preprint:

Ozga AT, Webster TH, Gilby IC, Wilson MA, Nockerts RS, Wilson ML, Pusey AE,
Li Y, Hahn BH, Stone AC. 2020. Urine as a high quality source of host genomic
DNA from wild populations.

The manuscript has been submitted for peer-review and we'll update the citation
when it's published.

## Goals
We aim to explicitly test the success of DNA capture and
next-generation sequencing using noninvasive sources of DNA from wild chimpanzees from Gombe National Park, Tanzania.

## Samples
This project has 26 total samples from:
- 10 individuals
- 4 sources (urine, feces, dental calculus, dentin)
- 4 capture methods (human exome, PTS whole genome, human whole genome, shotgun)

For the capture methods:
- "human exome" ("Exome" below) uses iDT XGEN human exome baits
- "PTS_WG" are a custom whole genome MYBaits kit produced by Arbor Biosciences using Pan troglodytes schweinfurthii DNA
- "Human_WG" is a MYcroarray whole genome capture kit (Human Masai male DNA)
- "Shotgun" is shotgun sequencing without capture

| Individual | Capture | Urine | Feces | Dental calculus | Dentin |
| ---------- | ------- | ----- | ----- | --------------- | ------ |
| 7365 | Exome | Y | Y | N | N |
| 7365 | PTS_WG | N | Y | N | N |
| 7365 | Human_WG | N | Y | N | N |
| 7150 | Exome | Y | Y | N | N |
| 7150 | PTS_WG | N | Y | N | N |
| 7150 | Human_WG | N | Y | N | N |
| 7069 | Exome | N | Y | N | N |
| 7069 | PTS_WG | N | Y | N | N |
| 7069 | Human_WG | N | Y | N | N |
| 7057 | Exome | N | N | Y | Y |
| 7057 | Human_WG | N | N | Y | N |
| 7057 | Shotgun | N | N | Y | N |
| 7072 | Exome | Y | N | N | N |
| 7535 | Exome | Y | N | N | N |
| 7433 | Exome | N | N | Y | N |
| 7433 | Human_WG | N | N | Y | N |
| 7433 | PTS_WG | N | N | Y | N |
| 7650 | Exome | Y | N | N | N |
| 7507 | Exome | Y | Y | N | N |
| 7507 | Human_WG | N | Y | N | N |
| 7507 | PTS_WG | N | Y | N | N |
| 7323 | Exome | Y | N | N | N |


## Setup


## Running the pipeline
