"""
데이터 입출력 모듈
"""
import pandas as pd
import numpy as np
# from ..utils import profile_step

# LOAD PLINK BFILE
def load_fam(fam_file):
    """
    - FID: Family ID
    - IID: Individual ID
    - PID: Father ID
    - MID: Mother ID
    - SEX: Sex (1: Male, 2: Female, 0 or -9: Unknown)
    """
    
    # 컬럼 이름 정의
    columns = ['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENO']
    
    # 공백으로 구분된 파일 읽기
    fam_data = pd.read_csv(fam_file, sep=r'\s+', header=None, names=columns)
    
    return fam_data

def load_bim(bim_file):
    """
    - CHR: Chromosome
    - SNP: SNP ID
    - CM: Genetic map position
    - BP: Base-pair position
    - A1: First allele
    - A2: Second allele
    """
    # 컬럼 이름 정의
    columns = ['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2']
    
    # 공백으로 구분된 파일 읽기
    bim_data = pd.read_csv(bim_file, sep=r'\s+', header=None, names=columns)
    
    return bim_data

def load_afreq(afreq_file):
    """
    - #CHROM: Chromosome
    - ID: SNP ID
    - REF: Reference allele
    - ALT: Alternate allele
    - ALT_FREQS: Alternate allele frequency
    - OBS_CT: Number of observations
    """
    # columns = ["#CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT"]
    
    afreq_data = pd.read_csv(afreq_file, sep=r'\s+')
    # afreq_data.columns = columns
    
    return afreq_data

# simplest version of read BED
def load_bed(bed_file, n_samples, n_snps):
    """
    PLINK BED 파일을 읽어서 유전형 행렬을 반환하는 간단한 함수
    
    Parameters:
    -----------
    bed_file : str
        .bed 파일 경로
    fam_file : str  
        .fam 파일 경로 (샘플 수를 알기 위해)
    n_snps : int
        SNP 개수
        
    Returns:
    --------
    numpy.ndarray
        유전형 행렬 (samples x SNPs)
        0=aa, 1=Aa, 2=AA, -9=missing
    """
    # 결과 배열 초기화
    genotypes = np.full((n_samples, n_snps), -9, dtype=np.int8)
    
    # BED 파일 읽기
    with open(bed_file, 'rb') as f:
        # 처음 3바이트(매직넘버) 건너뛰기
        f.seek(3)
        
        bytes_per_snp = (n_samples + 3) // 4
        
        for snp_idx in range(n_snps):
            # SNP 데이터 읽기
            snp_bytes = f.read(bytes_per_snp)
            
            # 각 샘플의 유전형 디코딩
            for sample_idx in range(n_samples):
                byte_idx = sample_idx // 4
                bit_offset = (sample_idx % 4) * 2
                
                if byte_idx < len(snp_bytes):
                    # 2비트 추출
                    geno_code = (snp_bytes[byte_idx] >> bit_offset) & 3
                    
                    # PLINK 인코딩 변환
                    if geno_code == 0:
                        genotypes[sample_idx, snp_idx] = 2    # AA
                    elif geno_code == 1:
                        genotypes[sample_idx, snp_idx] = -9   # missing
                    elif geno_code == 2:
                        genotypes[sample_idx, snp_idx] = 1    # Aa
                    elif geno_code == 3:
                        genotypes[sample_idx, snp_idx] = 0    # aa
    
    return genotypes

# Load BED chunk (fastest)
def load_bed_chunk_vectorized(
    bed_file, 
    n_samples, 
    n_snps, 
    keep_sample_indices=None, 
    keep_snp_indices=None
    ):
    """
    NumPy 벡터화를 사용한 초고속 BED 읽기
    
    핵심 개선:
    1. 전체 SNP 블록을 한 번에 읽기
    2. NumPy 배열 연산으로 비트 디코딩
    3. 최소한의 Python 루프
    """
    
    if keep_sample_indices is None:
        keep_sample_indices = np.arange(n_samples)
    else:
        keep_sample_indices = np.array(keep_sample_indices)
        
    if keep_snp_indices is None:
        keep_snp_indices = np.arange(n_snps)
    else:
        keep_snp_indices = np.array(keep_snp_indices)
    
    n_selected_samples = len(keep_sample_indices)
    n_selected_snps = len(keep_snp_indices)
    bytes_per_snp = (n_samples + 3) // 4
    
    # 결과 배열
    genotypes = np.full((n_selected_samples, n_selected_snps), -9, dtype=np.int8)
    
    # 샘플 바이트 정보 미리 계산
    sample_byte_indices = keep_sample_indices // 4
    sample_bit_offsets = (keep_sample_indices % 4) * 2
    
    with open(bed_file, 'rb') as f:
        f.seek(3)
        
        # 큰 청크 단위로 처리 (메모리 사용량 조절)
        chunk_size = min(1000, n_selected_snps)  # 1000 SNP씩 처리
        
        for chunk_start in range(0, n_selected_snps, chunk_size):
            chunk_end = min(chunk_start + chunk_size, n_selected_snps)
            chunk_snp_indices = keep_snp_indices[chunk_start:chunk_end]
            
            # 연속된 SNP 범위 찾기
            min_snp = np.min(chunk_snp_indices)
            max_snp = np.max(chunk_snp_indices)
            
            # 범위 전체를 읽기
            range_start_byte = min_snp * bytes_per_snp
            range_size_snps = max_snp - min_snp + 1
            range_size_bytes = range_size_snps * bytes_per_snp
            
            f.seek(3 + range_start_byte)
            range_data = f.read(range_size_bytes)
            
            if len(range_data) != range_size_bytes:
                continue
            
            # NumPy 배열로 변환
            byte_array = np.frombuffer(range_data, dtype=np.uint8)
            
            # 각 SNP 처리
            for i, snp_idx in enumerate(chunk_snp_indices):
                local_snp_offset = (snp_idx - min_snp) * bytes_per_snp
                snp_bytes = byte_array[local_snp_offset:local_snp_offset + bytes_per_snp]
                
                # 벡터화된 디코딩
                result_col = chunk_start + i
                _decode_snp_vectorized(snp_bytes, sample_byte_indices, 
                                     sample_bit_offsets, genotypes[:, result_col])
    
    return genotypes

def _decode_snp_vectorized(snp_bytes, byte_indices, bit_offsets, output_col):
    """NumPy 벡터화로 SNP 디코딩"""
    
    # 유효한 바이트 인덱스만 처리
    valid_mask = byte_indices < len(snp_bytes)
    valid_byte_indices = byte_indices[valid_mask]
    valid_bit_offsets = bit_offsets[valid_mask]
    
    if len(valid_byte_indices) == 0:
        return
    
    # 해당 바이트들 추출
    selected_bytes = snp_bytes[valid_byte_indices]
    
    # 비트 시프트로 genotype 코드 추출
    geno_codes = (selected_bytes >> valid_bit_offsets) & 3
    
    # 벡터화된 변환
    output_col[valid_mask] = np.where(geno_codes == 0, 2,     # AA
                             np.where(geno_codes == 1, -9,    # missing  
                             np.where(geno_codes == 2, 1,     # Aa
                                      0)))                     # aa

     
# SAVE GRM FILE
def save_grm_files(grm_matrix, sample_ids, out_prefix, n_snps_used):
    """
    GRM을 GCTA 호환 파일 형식으로 저장
    
    Parameters:
    -----------
    grm_matrix : numpy.ndarray
        GRM 행렬 (n_samples x n_samples)
    sample_ids : list of tuples
        [(FID, IID), ...] 형태의 샘플 ID 목록
    out_prefix : str
        출력 파일 접두사
    n_snps_used : int
        GRM 계산에 사용된 SNP 개수
        
    Returns:
    --------
    bool
        저장 성공 여부
    """
    import struct
    
    try:
        n_samples = len(sample_ids)
        
        # print(f"| GRM 파일 저장: {n_samples}개체, {n_snps_used}개 SNP 사용")
        
        # 1. .grm.id 파일 저장 (개체 ID 목록)
        id_file = f"{out_prefix}.grm.id"
        print(f"|   Save: {id_file}")
        with open(id_file, 'w') as f:
            for fid, iid in sample_ids:
                f.write(f"{fid}\t{iid}\n")
        
        # 2. .grm.bin 파일 저장 (하삼각 행렬, float32)
        bin_file = f"{out_prefix}.grm.bin"
        print(f"|   Save: {bin_file}")
        with open(bin_file, 'wb') as f:
            for i in range(n_samples):
                for j in range(i + 1):  # 하삼각 행렬 (대각선 포함)
                    value = float(grm_matrix[i, j])
                    f.write(struct.pack('<f', value))  # little-endian float32
        
        # 3. .grm.N.bin 파일 저장 (각 쌍의 SNP 개수)
        n_file = f"{out_prefix}.grm.N.bin"
        print(f"|   Save: {n_file}")
        with open(n_file, 'wb') as f:
            for i in range(n_samples):
                for j in range(i + 1):  # 하삼각 행렬
                    # 모든 쌍에서 동일한 SNP 개수 사용 (단순화)
                    snp_count = float(n_snps_used)
                    f.write(struct.pack('<f', snp_count))  # little-endian float32
        
        return True
        
    except Exception as e:
        print(f"❌ GRM 파일 저장 오류: {e}")
        return False


# LOAD GRM FILES
def load_grm_files(grm_prefix):
    """
    GCTA 형식의 GRM 파일들을 불러오기
    
    Parameters:
    -----------
    grm_prefix : str
        GRM 파일 접두사 (예: "output.add" -> "output.add.grm.bin", "output.add.grm.id", "output.add.grm.N.bin")
        
    Returns:
    --------
    dict
        - 'grm_matrix': numpy.ndarray - GRM 행렬 (n_samples x n_samples)
        - 'sample_ids': list of tuples - [(FID, IID), ...] 형태의 샘플 ID
        - 'n_snps': numpy.ndarray - 각 쌍별 사용된 SNP 개수 행렬
    """
    import struct
    import os
    
    # 파일 경로 설정
    id_file = f"{grm_prefix}.grm.id"
    bin_file = f"{grm_prefix}.grm.bin"
    n_file = f"{grm_prefix}.grm.N.bin"
    
    # 파일 존재 여부 확인
    for file_path in [id_file, bin_file, n_file]:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
    
    # print(f"|   Loading GRM files with prefix: {grm_prefix}")
    
    # 1. .grm.id 파일 읽기 (샘플 ID)
    sample_ids = []
    with open(id_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                fid, iid = parts
                sample_ids.append((fid, iid))
    
    n_samples = len(sample_ids)
    print(f"|   Loaded {n_samples} sample IDs")
    
    # 2. .grm.bin 파일 읽기 (GRM 행렬)
    grm_matrix = np.zeros((n_samples, n_samples), dtype=np.float64)
    
    with open(bin_file, 'rb') as f:
        for i in range(n_samples):
            for j in range(i + 1):  # 하삼각 행렬
                data = f.read(4)  # float32 = 4 bytes
                if len(data) == 4:
                    value = struct.unpack('<f', data)[0]  # little-endian float32
                    grm_matrix[i, j] = value
                    grm_matrix[j, i] = value  # 대칭 행렬
                else:
                    raise ValueError(f"Unexpected end of file in {bin_file}")
    
    print(f"|   Loaded GRM matrix ({n_samples} x {n_samples})")
    
    # 3. .grm.N.bin 파일 읽기 (SNP 개수 행렬)
    n_snps_matrix = np.zeros((n_samples, n_samples), dtype=np.float64)
    
    with open(n_file, 'rb') as f:
        for i in range(n_samples):
            for j in range(i + 1):  # 하삼각 행렬
                data = f.read(4)  # float32 = 4 bytes
                if len(data) == 4:
                    value = struct.unpack('<f', data)[0]  # little-endian float32
                    n_snps_matrix[i, j] = value
                    n_snps_matrix[j, i] = value  # 대칭 행렬
                else:
                    raise ValueError(f"Unexpected end of file in {n_file}")
    
    print(f"|   Loaded SNP count matrix")
    
    return {
        'grm_matrix': grm_matrix,
        'sample_ids': sample_ids,
        'n_snps': n_snps_matrix
    }


# parse `.hsq` file
def load_factor_hsq(hsq_file):
    def parsing_factor_hsq(hsq_file):
        res = {}
        
        df = pd.read_csv(hsq_file, sep="\t")
        V_cols = [
            c for c in df["Source"].to_list() if c.startswith("V(") and "/Vp" not in c
            ]
        
        for V_col in V_cols:
            Var = df.loc[df["Source"] == V_col, "Variance"].values[0]
            Var_se = df.loc[df["Source"] == V_col, "SE"].values[0]
            res[f"{V_col}"] = Var
            res[f"{V_col}_se"] = Var_se

        return pd.DataFrame(res, index=[0])
    
    def make_hsqF(df_factor):
        res = {}
        
        # Variance components
        V_het = df_factor["V(G1)"].values  # heterozygote variance
        V_cross = df_factor["V(G2)"].values  # hetero X homo covariance
        V_hom = df_factor["V(G3)"].values  # homozygote variance
        V_e = df_factor["V(e)"].values  # error variance
        
        # Standard errors
        SE_het = df_factor["V(G1)_se"].values
        SE_cross = df_factor["V(G2)_se"].values
        SE_hom = df_factor["V(G3)_se"].values
        SE_e = df_factor["V(e)_se"].values
        
        # Total phenotypic variance
        V_p = V_het + V_hom + V_e
        V_G = V_het + V_hom
        
        # Calculate estimates
        h2_het = V_het / V_p  # Heterozygote heritability
        h2_hom = V_hom / V_p  # Homozygote heritability
        h2_FACTOR = h2_het + h2_hom  # Total heritability
        rho_hh = V_cross / np.sqrt(V_het * V_hom)  # Genetic correlation
        ratio_hh = h2_hom / h2_het  # Ratio
        
        # ===== Delta Method Calculations =====
        
        # 1. Standard Error for rho(het,hom) = V(G2) / sqrt(V(G1) * V(G3))
        sqrt_V1V3 = np.sqrt(V_het * V_hom)
        dCov_dV1 = -V_cross / (2 * V_het**(3/2) * V_hom**(1/2))
        dCov_dV2 = 1 / sqrt_V1V3
        dCov_dV3 = -V_cross / (2 * V_het**(1/2) * V_hom**(3/2))
        
        SE_rho_hh = np.sqrt(
            (dCov_dV1**2) * (SE_het**2) +
            (dCov_dV2**2) * (SE_cross**2) +
            (dCov_dV3**2) * (SE_hom**2)
        )
        
        # 2. Standard Error for h2(het) = V(G1) / V_p
        dh2het_dV1 = (V_hom + V_e) / (V_p**2)
        dh2het_dV3 = -V_het / (V_p**2)
        dh2het_dVe = -V_het / (V_p**2)
        
        SE_h2_het = np.sqrt(
            (dh2het_dV1**2) * (SE_het**2) +
            (dh2het_dV3**2) * (SE_hom**2) +
            (dh2het_dVe**2) * (SE_e**2)
        )
        
        # 3. Standard Error for h2(hom) = V(G3) / V_p
        dh2hom_dV1 = -V_hom / (V_p**2)
        dh2hom_dV3 = (V_het + V_e) / (V_p**2)
        dh2hom_dVe = -V_hom / (V_p**2)
        
        SE_h2_hom = np.sqrt(
            (dh2hom_dV1**2) * (SE_het**2) +
            (dh2hom_dV3**2) * (SE_hom**2) +
            (dh2hom_dVe**2) * (SE_e**2)
        )
        
        # 4. Standard Error for hom/het ratio = V(G3) / V(G1)
        dratio_dV1 = -V_hom / (V_het**2)
        dratio_dV3 = 1 / V_het
        
        SE_ratio_hh = np.sqrt(
            (dratio_dV1**2) * (SE_het**2) +
            (dratio_dV3**2) * (SE_hom**2)
        )
        
        # 5. Standard Error for h2(FACTOR) = h2(het) + h2(hom)
        dh2factor_dV1 = V_e / (V_p**2)
        dh2factor_dV3 = V_e / (V_p**2)
        dh2factor_dVe = -V_G / (V_p**2)
        
        SE_h2_FACTOR = np.sqrt(
            (dh2factor_dV1**2) * (SE_het**2) +
            (dh2factor_dV3**2) * (SE_hom**2) +
            (dh2factor_dVe**2) * (SE_e**2)
        )
        
        # Add results to dataframe
        res["rho"] = rho_hh
        res["rho_se"] = SE_rho_hh
        res["h2(het)"] = h2_het
        res["h2(het)_se"] = SE_h2_het
        res["h2(hom)"] = h2_hom
        res["h2(hom)_se"] = SE_h2_hom
        res["hom/het"] = ratio_hh
        res["hom/het_se"] = SE_ratio_hh
        res["h2(F)"] = h2_FACTOR
        res["h2(F)_se"] = SE_h2_FACTOR
        
        return pd.DataFrame(res, index=[0])
        # return df_factor[["pheno", 
        #                 "h2(het)", "h2(het)_se", 
        #                 "h2(hom)", "h2(hom)_se",
        #                 "Cov(hh)", "Cov(hh)_se",
        #                 "hom/het", "hom/het_se",
        #                 "h2(FACTOR)", "h2(FACTOR)_se"]]

    df_hsqF = make_hsqF(parsing_factor_hsq(hsq_file))
    return df_hsqF
   
if __name__ == "__main__":
    pass
    

