# Lineage-Tracing-Simulation

Simulation of CRISPR-Cas9 Editing on Evolving Barcode and Accuracy of Lineage Tracing. The folder 'code' contains all the code for simulation. The folders `data/plots` contain the raw data and plots generated by our simulations.

This project is under collaboration with Dr. Yipeng Yang, Department of Mathematics and Statistics, University of Houston (yangy@uhcl.edu).


## References

   **Liu, F., et al. (2024).**  
   "Simulation of CRISPR-Cas9 editing on evolving barcode and accuracy of lineage tracing."  
   *Scientific Reports*, 14(1): 19213.  
   [Link to Article](https://doi.org/10.1038/s41598-024-19213-4)

   *Summary*: We designed a simulation program that mimics CRISPR-Cas9 editing on evolving barcodes and the double-strand break repair process during cell divisions. Barcode mutations build sequentially with each generation, resulting in unique mutation profiles for each cell. We sampled barcodes from leaf cells, reconstructed the lineage, and compared it to the original tree to test algorithm accuracy under different parameters. Our simulations highlight the importance of sampling size, barcode length, indel probabilities, and Cas9 activity for accurate lineage tracing. We found that sampling size and indel probabilities significantly impact accuracy, with large segment deletions in early generations potentially reducing lineage accuracy. These results offer recommendations for optimizing Cas9-mediated molecular barcodes in experiments.


## File Description for the Simulation Code

<details>
<summary>How to run the simulations?</summary>
<br>

1. Please put all code/main program/ files in the same folder to run the program.

2. Please run `main.m` to start this program.

3. The first portion in `main.m` is on parameter setting.

4. The second portion in `main.m` is on simulation of a single barcode:
    - `funbarnew` - RMP method
    - `funbarnewRMPNF` - RMPNF method
    - `funbarNBJ` - NBJ method
    - `funbarNBJNF` - NBJNF method

    Calculated lineage accuracy is stored in the corresponding vector `accuracy***`. When a certain method is used, please ensure to display the correct vector.

5. The third portion in `main.m` is on simulation of two barcodes:
    - `funbarNBJ2` - NBJ method
    - `funbarNBJ2trbk` - NBJNF method

    This portion requires one more parameter `propmi`, which stands for proportion of match on individual barcode.

6. One may modify the second or third portion to run a method multiple times at once. The following is an example to run the RMP method 10 times:

    ```matlab
    mcit = 10; % Monte Carlo simulation iteration number

    accuracyRMP = zeros(mcit, 2);  

    for i = 1:mcit
        [accuracyRMP(i,1), accuracyRMP(i,2)] = funbarnew(n, it, propm, ss, mupb, ins_sub, lgdelprob, divp, clive, pulse, trbk); % RMP method
    end

    disp('RMP - Accuracy of all internal nodes');
    accuracyRMP(:,1)

    disp('RMP - Accuracy of all dividing nodes');
    accuracyRMP(:,2)
    ```


</details>

<details>
<summary>Program/Code Function Descriptions </summary>
<br>
    
### btree.m
This file/function contains the implementation of a binary tree, a function to manipulate a binary tree structure.

### findm.m
This file/function finds and counts the matches between two matrices, specifically counting the total matched nodes and the total matched splitting nodes.

### funbarNBJ.m
This file/function using NBJ method for single barcode. It compares the similarity between cells in a binary tree structure, reconstructs the lineage tree based on similarity scores, and evaluates the alignment and matching of nodes within the tree.

### funbarNBJ2.m
**Functionality**: Similar to `funbarNBJ`, but using NBJ method for two indenpendent barcodes.
**Difference**: `funbarNBJ2` includes enhancements or modifications to the original `funbarNBJ` method, involving different similarity computation techniques and reconstruction steps.

### funbarNBJ2trbk.m
**Functionality**: An extension of `funbarNBJ2`, is using NBJNF method for 2 barcodes

### funbarNBJNF.m
**Functionality**: Similar to `funbarNBJ` but without certain features or enhancements present in `funbarNBJ2` or `funbarNBJ2trbk`.
**Difference**: A simplified or variant version of `funbarNBJ`, possibly without the features or modifications introduced in the subsequent versions.

### funbarnew.m
**Functionality**: This function performs similarity comparison between cells in a binary tree structure, reconstructs the lineage tree based on similarity scores, and evaluates the alignment and matching of nodes within the tree, including handling of collapsing columns and comparison among the leaves and nodes.
**Similarities**: Both `funbarnew` and `funbarNBJ` compare the similarity between cells in a binary tree structure, reconstruct the lineage tree, and evaluate node matching and alignment.
**Differences**: `funbarnew` includes additional steps for handling collapsing columns and provides more detailed evaluations among leaves and nodes, along with more complex reconstruction logic and detailed metrics for the lineage tree.

### funbarnewRMPNF.m
Similar to `funbarNBJNF` but uses the RMPNF method.

### genN.m
This file/function assigns generation numbers to each cell in a sample based on their similarity scores and relative positions within a tree structure.

### genN_unsorted.m
**Functionality**: This file/function assigns generation numbers without sorting cells based on similarity to the root, using random sampling and considering live leaves.
**Differences**: Unlike `genN`, which may sort cells by similarity to the root, `genN_unsorted` skips this sorting step and focuses on a random sample of live cells.

### genN_unsorted2.m
**Functionality**: This file/function assigns generation numbers to cells without sorting based on similarity to the root, focusing on a random sample of live leaves with a different structure for `checklive`.
**Differences**: Similar to `genN_unsorted`, it skips sorting by similarity to the root and uses a different structure for checking live cells, which influences how the generation numbers are assigned.

### genchild.m
This function generates a child barcode from a parent barcode by simulating mutations, cuts, insertions, deletions, and substitutions within the specified range.

### trim.m
This function trims the leading and trailing zeros from a vector and returns the resulting vector.

### valcheck.m
This function checks the validity of cell numbers in each generation of a lineage tree by ensuring the total number of cells does not exceed the maximum possible number at the deepest level.

### vcomp5.m and vcomp_cons.m
Both functions compare the similarity between two vectors `s` and `t`, and return a similarity score (`y`) and an optional matching string (`ost`). `vcomp5` performs similarity comparisons between vectors using a simpler 5x5 scoring matrix without handling consecutive matches or mismatches, while `vcomp_cons` uses a more complex 4x4 scoring matrix with additional handling for consecutive matches and mismatch penalties.

</details>
