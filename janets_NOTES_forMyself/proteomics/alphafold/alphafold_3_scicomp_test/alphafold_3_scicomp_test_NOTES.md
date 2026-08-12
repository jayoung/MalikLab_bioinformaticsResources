# August 2026

https://sciwiki.fredhutch.org/compdemos/alphafold3/

https://github.com/google-deepmind/alphafold3?tab=readme-ov-file#installation-and-running-your-first-prediction


I created an example input json file as shown [here](https://github.com/google-deepmind/alphafold3?tab=readme-ov-file#installation-and-running-your-first-prediction).  

The sequence they use is 298aa long and is 100% identical to PDB [2PV7_A](https://www.rcsb.org/structure/2PV7) except that the PDB entry contains a bunch of X residues, which are selenomethionines ('MSE' "Selenomethionine is frequently used in X-ray crystallography for experimental phasing").   "Crystal structure of chorismate mutase / prephenate dehydrogenase (tyrA) (1574749) from Haemophilus influenzae RD at 2.00 A resolution". I think they substituted M for each selenomethionine

```
cd ~/FH_fast_storage/git_more_repos/MalikLab_bioinformaticsResources/janets_NOTES_forMyself/proteomics/alphafold/alphafold_3_scicomp_test

sbatch af3.sh
    # took ~15 mins
```

Annotated version of the input json:
```
{
  ## job name
  "name": "2PV7",
  ## you can specify >1 item to fold, 
  ## each item can be Protein | RNA | DNA | Ligand
  ## it's also possible to specify post-translational modifications, DNA/RNA modifictions, structural templates 
  ## can supply multiple sequence alignments but if you don't it'll generate them
  "sequences": [
    {
      "protein": {
        ## supplying two ids and one sequence implies a homomeric chain with two copies
        "id": ["A", "B"],
        "sequence": "GMRESYANENQFGFKTINSDIHKIVIVGGYGKLGGLFARYLRASGYPISILDREDWAVAESILANADVVIVSVPINLTLETIERLKPYLTENMLLADLTSVKREPLAKMLEVHTGAVLGLHPMFGADIASMAKQVVVRCDGRFPERYEWLLEQIQIWGAKIYQTNATEHDHNMTYIQALRHFSTFANGLHLSKQPINLANLLALSSPIYRLELAMIGRLFAQDAELYADIIMDKSENLAVIETLKQTYDEALTFFENNDRQGFIDAFHKVRDWFGDYSEQFLKESRQLLQQANDLKQG"
      }
    }
  ],
  ## you can provide >1 random seed and you'll get that number of predicted structures, each using the seed provided
  "modelSeeds": [1],
  "dialect": "alphafold3",
  ## JSON version - can be 1-4 (anything >1 allows specification of more things)
  "version": 1
}
```

Output - see [docs](https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md)

1 seed used, it outputs 5 samples from that seed (0-4)

The `ranking_scores.csv` file contains scores for each of those samples

The `2pv7_model.cif` is the top-ranking model file that can be opening in ChimeraX etc

this time, sample 4 has the highest ranking score, and the output files in the top-level folder are copies of the outupt files in the sample 4 directory:
```
diff 2pv7_model.cif seed-1_sample-4/model.cif 
diff 2pv7_confidences.json seed-1_sample-4/confidences.json
diff 2pv7_summary_confidences.json seed-1_sample-4/summary_confidences.json
```