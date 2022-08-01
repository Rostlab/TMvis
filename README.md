# **_TMvis_**  :atom::dna:

Welcome to **_TMvis_** - a pipeline for transmembrane protein annotation and 3D visualization.

![alt-text](https://github.com/stelviopas/TMPipeline/blob/main/P45880_plddt.png)               
                                                    
**_TMvis_** combines AlphaFold 2 [1] structures from the AlphaFold DB [2] with predicted transmembrane protein (TMP) annotations into interactive 3D visualizations of protein structures embedded into membranes. The TMPs are predicted by TMbed [3], a method based on the protein language model ProtT5 [4], which provides per-residue alpha-helical and beta-barrel transmembrane segment predictions. The respective AlphaFold 2 TMP structures are then enhanced by adding the predicted TMbed topology to the 3D visualization. Further, **_TMvis_** allows to add membrane embeddings predicted by ANVIL [5], or PPM3 [6]. 

![alt-text](https://zenodo.org/api/iiif/v2/3c577630-8539-4e9e-a248-5bc7fd71ec35:7e67b846-f6c8-4df8-a3b1-c5438ea954c6:TMvis_project_overview.png/full/750,/0/default.png)

## Quick access

As an example, we provide a subset of 496 predicted TMPs. TMbed predicted 4.967 TMP for the human proteome (20,375 proteins, UniProt [7] version April 2022; excluding TITIN_HUMAN due to length). We extracted AlphaFold 2 structures with an average per-residue confidence score (pLDDT) of more than 90%, which lead to the subset. Please download the set [here](https://doi.org/10.5281/zenodo.6816083). 

Once you have the dataset ready, you can use the Jupyter notebook `TMvis.ipynb` in the `TMvis` folder for 3D-visualization of every protein structure in the dataset with predictions of ANVIL, PPM3, and TMbed. Additionally, you can visualize the per-residue confidence scores (pLDDT) of AlphaFold.


## Running **_TMvis_** :atom::dna:

### Step 0: Clone repository
`git clone https://github.com/Rostlab/TMvis`

### Step 1: Install dependencies 
1. Python and Conda 
```
conda env create -n TMvis --file TMvis.yml    
conda activate TMvis
```       
2. Nextflow
3. Docker

### Step 3: Generate dataset
See `data/` folder for details.
Make sure that after this step, you have a folder `data/current/` containing one folder with your AlphaFold 2 structures and a text file with TMbed predictions.

### Step 4: Preprocessing and subset extraction
Run ```python3 ./TMvis/main.py```

`main.py` will generate a `results/db` folder which contains AF structure of TMbed predicted membrane alpha/beta proteins. Additionally, `db/**/pLDDT90F1` is a subset from `db` with alpha and beta proteins selected by following criteria: 
- Max. 2.700 base pairs long (length of one AlphaFold 2 PDB file) 
- pLDDT mean score per protein is higher than 90 (highly accurate structures). If needed you can change the threshold.

### Step 5: Predict membrane with PPM3  
1. Unpack the docker container. See `docker/` folder for details.
2. Run ```nextflow run ./PPM3/run_PPM3.nf -c custom.config```

### Step 6 (optional): Predict membrane with ANVIL
Note: to run ANVIL, you need an access key.
1. Unpack the docker container. See `docker/` folder for details on how to do that and where to get the access key.
2. Run ```nextflow run ./ANVIL/run_anvil.nf -c custom.config```

### Step 7: Visualize TMPs in 3D     
Run ```jupyter notebook ./TMvis/TMvis.ipynb```

## References 
1. AlphaFold - AlphaFold Jumper, John, Richard Evans, Alexander Pritzel, Tim Green, Michael Figurnov, Olaf Ronneberger, Kathryn Tunyasuvunakool, et al. 2021. “Highly Accurate Protein Structure Prediction with AlphaFold.” Nature 596 (7873): 583–89.
2. Alphafold DB - Varadi, Mihaly, Stephen Anyango, Mandar Deshpande, Sreenath Nair, Cindy Natassia, Galabina Yordanova, David Yuan, et al. 2022. “AlphaFold Protein Structure Database: Massively Expanding the Structural Coverage of Protein-Sequence Space with High-Accuracy Models.” Nucleic Acids Research 50 (D1): D439–44.
3. TMbed - TMbed Bernhofer, Michael, and Burkhard Rost. 2022. “TMbed – Transmembrane Proteins Predicted through Language Model Embeddings.” bioRxiv.
4. ProtT5 - A. Elnaggar et al., "ProtTrans: Towards Cracking the Language of Lifes Code Through Self-Supervised Deep Learning and High Performance Computing," in IEEE Transactions on Pattern Analysis and Machine Intelligence, doi: 10.1109/TPAMI.2021.3095381.
5. ANVIL - ANVIL Postic, Guillaume, Yassine Ghouzam, Vincent Guiraud, and Jean-Christophe Gelly. 2016. “Membrane Positioning for High- and Low-Resolution Protein Structures through a Binary Classification Approach.” Protein Engineering, Design & Selection: PEDS 29 (3): 87–91.
6. PPM3 - PPM3 Lomize, Mikhail A., Irina D. Pogozheva, Hyeon Joo, Henry I. Mosberg, and Andrei L. Lomize. 2012. “OPM Database and PPM Web Server: Resources for Positioning of Proteins in Membranes.” Nucleic Acids Research 40 (Database issue): D370–76.
7. UniProt - UniProt Consortium (2021). UniProt: the universal protein knowledgebase in 2021. Nucleic acids research, 49(D1), D480–D489.
