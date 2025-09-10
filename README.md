# pyGCTA

`pyGCTA` is a Python package designed to create genetic relatedness matrices (GRMs) based on a factor model. It supports the generation of Factor GRMs (heterozygous, homozygous, and their interaction) from PLINK bfiles and performs REML analysis using GCTA. Users can also parse the results of REML analysis directly within this package.

---

## Key Features

1. **Factor GRM Creation**  
   - Generate Factor GRMs (heterozygous, homozygous, and interaction terms) from PLINK bfiles.

2. **REML Analysis with GCTA**  
   - Perform REML analysis using GCTA.  
   - Users can specify the path to the GCTA executable. Alternatively, users can run REML externally and use this package solely for parsing results.

3. **Preprocessing and Usage Examples**  
   - Basic preprocessing pipeline for input bfiles is provided in `notebooks/plink.ipynb`.  
   - Detailed usage examples, including GRM creation and REML analysis, are available in `notebooks/basic_usage.ipynb`.

---

## Installation

1. **Set Up Environment**  
   - Use the `environment.yml` file to configure the Python environment:
   ```bash
   conda env create -f environment.yml
   conda activate pygcta
   ```

2. **Dependencies**  
   - Ensure that PLINK and GCTA executables are available.  
   - Download PLINK and GCTA from their official sites:  
     - [PLINK](https://www.cog-genomics.org/plink/)  
     - [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#Download)  

---

## Usage

### 1. Preprocess Input Bfiles
Follow the preprocessing pipeline in `notebooks/plink.ipynb` to prepare your input bfiles.

### 2. Create Factor GRMs and Perform REML Analysis
```python
from gcta import GCTA

gcta = GCTA()

# Compute GRM
gcta.calculate_grm(
   bfile="<plink-bfile-header>",
   standardization="factor",
   afreq="<plink-afreq-file>",
   maf_threshold=0.01,
)

# Save GRM
gcta.save_grm(out_prefix="<output-file>")
```
Refer to `notebooks/basic_usage.ipynb` for detailed instructions on creating Factor GRMs.

### 3. Perform REML Analysis
You can perform REML analysis in two ways:  
1. **Externally**: Run REML analysis using GCTA outside this package and use `pyGCTA` to parse the results.  
2. **Internally**: Provide the path to the GCTA executable, and `pyGCTA` will handle the REML execution for you.

```python
from gcta import GCTA

gcta = GCTA()

# run REML
gcta.reml(
   gcta_path="<gcta-path>",
   mgrm_file="<mgrm-file>", # V1: hetero-grm, V2:hetero-homo-grm, V3: homo-grm
   phenotype_file="<phenotype-file>",
   out_prefix="<output-prefix>",
)

# compute h2 and other parameters based on delta methods
df_hsqF = gcta.load_factor_hsq(".hsq-file", save=True)
```

---

## License

This project is licensed under the [MIT License](LICENSE).

---

## References

On the Fisher's Additive Model: Revisiting Heritability through Genotype Encoding, __under review__