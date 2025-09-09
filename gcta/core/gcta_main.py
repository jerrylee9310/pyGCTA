"""
Main GCTA class - Python implementation
"""
import os
import time
import psutil
from ..utils.profiling import profile_step
import pandas as pd
import numpy as np

class GCTA:
    def __init__(self, debug=False):
        """Initialize GCTA instance"""
        print("GCTA Python v0.1.0 initialized")
        self.version = "0.1.0"
        self.debug = debug
        self.profile = {}
        
        # PLINK 파일 관련 attributes
        self.FAM = None
        self.BIM = None
        self.BED = None
        self.AFREQ = None
        self.keep_snp_indices = None
        self.keep_sample_indices = None
        
        # GRM 관련 attributes
        self.genotype_std = {}
        self.GRM = {}
        
    def print_profile(
        self, 
        indent=4
        ):
        import json
        print(json.dumps(self.profile, indent=indent))

    def calculate_grm(
        self,   
        bfile, # plink file prefix
        standardization, # 'add', 'dom', 'factor'
        afreq=None, # .afreq file
        keep=None, # sample id file to keep (optional)
        maf_threshold=0.01,
        ):
        
        try:
            # step 1. load plink files
            self._load_plink_files(bfile)
            
            if afreq is not None:
                # step 1.1 load afreq files
                self._load_afreq_file(afreq, maf_threshold)
            
            if keep is not None:
                # step 1.2 load keep files
                self._load_keep_file(keep)
            
            # step 2. compute grm
            self._calculate_grm(
                bfile=bfile,
                standardization=standardization,
            )
            
            print(f"|   Successfully computed GRM!")
            
        except Exception as e:
            raise Exception(e)
        
    def save_grm(
        self,
        out_prefix,
        ):
        from ..core.data_io import save_grm_files
        # input validation
        if self.GRM is None:
            raise Exception("No GRM to save")
        
        # sample_ids 준비 (FID, IID 튜플 리스트)
        if self.keep_sample_indices is not None:
            # keep subset 사용한 경우
            sample_ids = [(self.FAM.iloc[i]['FID'], self.FAM.iloc[i]['IID']) 
                         for i in self.keep_sample_indices]
        else:
            # 전체 샘플 사용한 경우
            sample_ids = [(row['FID'], row['IID']) for _, row in self.FAM.iterrows()]
            
        # SNP 개수 (실제 사용된 SNP 개수는 간단히 전체 BIM으로 처리)
        n_snps_used = len(self.keep_snp_indices)
        
        for std_types in self.GRM.keys():
            try:
                save_grm_files(
                    self.GRM[std_types], 
                    sample_ids, 
                    f"{out_prefix}.{std_types}", 
                    n_snps_used)
            except Exception as e:
                raise Exception(f"Error in {std_types} GRM saving: {str(e)}")
        
        print(f"|   Successfully saved GRM!")
    
    def remove_grm(
        self,
        grm_type,
        ):
        # input validation
        if grm_type not in self.GRM.keys():
            raise Exception(f"No such GRM type: {grm_type}")
        
        del self.GRM[grm_type]
        # del self.genotype_std[grm_type]
            
    def load_grm(
        self,
        grm_prefix: str, # grm prefix
        ):
        try:
            self._load_grm(grm_prefix)
            
        except Exception as e:
            raise Exception(f"Error in GRM loading: {str(e)}")

        
        print(f"|   GRM matrix loaded successfully")
       
    def load_factor_hsq(
        self,
        factor_hsq_file: str,
        save=False,
        ):
        try:
            from ..core.data_io import load_factor_hsq
            print(f"|   Must be (V1: het, V2: hethom, V3: hom)")
            self.hsqF = load_factor_hsq(factor_hsq_file)
            
            if save:
                from tabulate import tabulate
                save_path = factor_hsq_file + "F"
                metrics = ["h2(het)", "h2(hom)", "h2(F)", "rho", "hom/het"]

                df_out = pd.DataFrame({
                    "metric": metrics,
                    "mean": [self.hsqF[m].iloc[0] for m in metrics],
                    "SE": [self.hsqF[f"{m}_se"].iloc[0] for m in metrics]
                })

                # save as table format
                tbl_str = tabulate(df_out, headers="keys", tablefmt="plain", floatfmt=".6f", showindex=False)
                with open(save_path, "w") as f:
                    f.write(tbl_str)

                print(tbl_str)
            else:
                print(self.hsqF)
        
        except Exception as e:
            raise Exception(f"Error in factor hsq loading: {str(e)}")
            
         
    def reml(
        self,
        gcta_path: str,
        mgrm_file: str,
        phenotype_file: str,
        out_prefix: str,
        covar_file=None, # optional
        no_contrain=True, # True: --reml-no-constrain
        no_lrt=True,
        n_threads=5,
        ):
        
        try:
            gcta_cmd = [
                gcta_path,
                "--reml",
                "--mgrm", mgrm_file,
                "--pheno", phenotype_file,
                "--out", out_prefix,
            ]
            if covar_file is not None:
                gcta_cmd.extend(["--covar", covar_file])
            if no_contrain:
                gcta_cmd.append("--reml-no-constrain")
            if no_lrt:
                gcta_cmd.append("--reml-no-lrt")
            if n_threads > 1:
                gcta_cmd.extend(["--thread-num", str(n_threads)])
            
            cmd = " ".join(gcta_cmd)
            os.system(cmd)
            
        except Exception as e:
            raise Exception(f"Error in REML: {str(e)}")
        
    
    @profile_step("load_plink", "Load PLINK files")
    def _load_plink_files(self, bfile):
        """
        Load PLINK binary files (.bed/.bim/.fam)
        """
        from ..core.data_io import load_fam, load_bim
        fam_file = f"{bfile}.fam"
        bim_file = f"{bfile}.bim"
        bed_file = f"{bfile}.bed"
        
        # input validation
        for file_path in [fam_file, bim_file, bed_file]:
            if not os.path.exists(file_path):
                raise Exception(f"No such file: {file_path}")
                
        # read FAM, BIM files
        self.FAM = load_fam(fam_file)
        self.BIM = load_bim(bim_file)
        print(f"|   Loaded {len(self.FAM)} samples and {len(self.BIM)} SNPs")
        # return {"status": "success", 
                # "message": "PLINK files loaded successfully"}
    
    @profile_step("load_afreq", "Load AFREQ file")
    def _load_afreq_file(self, afreq, maf_threshold):
        """
        Load AFREQ file
        """
        from ..core.data_io import load_afreq
        # input validation
        if not os.path.exists(afreq):
            raise Exception(f"No such file: {afreq}")
        
        # load AFREQ file
        df_afreq = load_afreq(afreq)
        print(f"|   Loaded {len(df_afreq)} SNPs allele frequencies")
        
        # MAF thresholding
        df_afreq_filtered = df_afreq[(df_afreq["ALT_FREQS"] > maf_threshold) & (df_afreq["ALT_FREQS"] < 0.5)]
        print(f"|   {len(df_afreq_filtered)} SNPs after MAF thresholding")
        
        # snp indices to keep
        merged = self.BIM.reset_index().merge(
            df_afreq_filtered, 
            left_on=['SNP', 'A2', 'A1'],
            right_on=["ID", "REF", "ALT"],
            how="inner")
        self.keep_snp_indices = merged['index'].tolist()
        self.AFREQ = merged
        print(f"|   Using {len(self.keep_snp_indices)} common SNPs")
    
    @profile_step("load_keep", "Load keep file")
    def _load_keep_file(self, keep_file):
        # input validation
        if not os.path.exists(keep_file):
            raise Exception(f"No such file: {keep_file}")
        
        # load keep file
        keep_df = pd.read_csv(keep_file, sep='\s+', names=['FID', 'IID'])
        print(f"|   Loaded {len(keep_df)} samples from keep file")
        
        # sample indices to keep
        merged = self.FAM.reset_index().merge(keep_df, on=['FID', 'IID'])
        self.keep_sample_indices = merged['index'].tolist()
        print(f"|   Using {len(self.keep_sample_indices)} common samples")
        
    @profile_step("calculate_grm", "Calculate GRM")
    def _calculate_grm(
        self,
        bfile, 
        standardization, 
        ):
        from ..core.data_io import load_bed_chunk_vectorized
        from ..analysis.grm import calculate_grm_by_std
        #input validation
        
        # load bed file
        bed_file = f"{bfile}.bed"
        
        if self.BED is None:
            print(f"|   Loading BED file")
            self.BED = load_bed_chunk_vectorized(
                bed_file=bed_file, 
                n_samples=len(self.FAM), 
                n_snps=len(self.BIM),
                keep_sample_indices=self.keep_sample_indices,
                keep_snp_indices=self.keep_snp_indices,
            )
        
        # keep effective genotypes
        genotypes_effective = self.BED
        n_samples_effective = genotypes_effective.shape[0]
        n_snps_effective = genotypes_effective.shape[1]
        # n_samples_effective = len(self.keep_sample_indices)
        # n_snps_effective = len(self.keep_snp_indices)
        print(f"|   Using {n_samples_effective} samples and {n_snps_effective} SNPs")
        
        # compute grm
        # move this part to `../analysis/grm.py`
        grms = calculate_grm_by_std(
            genotype=genotypes_effective,
            afreq=self.AFREQ["ALT_FREQS"].to_numpy(),
            standardization=standardization,
        )
        if standardization == "factor":
            self.GRM["het"] = grms[0]
            self.GRM["hom"] = grms[1]
            self.GRM["hethom"] = grms[2]
        else:
            self.GRM[standardization] = grms[0]
        
        
    
    @profile_step("load_grm", "Load GRM")
    def _load_grm(self, grm_prefix):
        from ..core.data_io import load_grm_files
        grm_data = load_grm_files(grm_prefix)
        std_type = grm_prefix.split(".")[-1]
        self.GRM_LOADED = grm_data["grm_matrix"]
        self.FAM_LOADED = grm_data["sample_ids"]
        self.N_SNP_LOADED = grm_data["n_snps"]
        
        
        
if __name__ == "__main__":
    pass