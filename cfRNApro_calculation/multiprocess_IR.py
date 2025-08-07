import numpy as np
import pandas as pd
import argparse
from concurrent.futures import ProcessPoolExecutor
import warnings
warnings.filterwarnings('ignore')

def compute_introns(exons):
    """计算内含子区域坐标"""
    exons_sorted = sorted(exons, key=lambda x: x[0])
    introns = []
    prev_end = exons_sorted[0][1]
    for exon in exons_sorted[1:]:
        if prev_end < exon[0]:
            introns.append((prev_end, exon[0]))
        prev_end = max(prev_end, exon[1])
    return introns

def get_region_coverage(full_coverage, region_mask):
    """根据掩码提取区域覆盖"""
    return full_coverage[region_mask]

def kl_divergence(P, Q):
    """KL散度计算 (增加数值稳定性)"""
    epsilon = 1e-10
    P = P / (P.sum() + epsilon)
    Q = Q / (Q.sum() + epsilon)
    return np.sum(P * np.log(P / Q + epsilon))

def process(_ref, _bed):
    results = {
        'name': [],
        'intron_ratio': []  # 内含子reads比例
    }
    
    # 按基因名分组处理
    for gene_name, group in _ref.groupby('name'):
        # 获取基因坐标信息
        exons = group[['start', 'end']].values
        gene_start = exons[:, 0].min()
        gene_end = exons[:, 1].max()
        gene_length = gene_end - gene_start
        strand = group['strand'].iloc[0]
        
        # 跳过无效基因
        if gene_length <= 0:
            continue
            
        # 生成全基因覆盖数组
        full_coverage = np.zeros(gene_length, dtype=int)
        valid_bed = _bed[(_bed['start'] < gene_end) & (_bed['end'] > gene_start)]
        for _, row in valid_bed.iterrows():
            start = max(row['start'], gene_start) - gene_start
            end = min(row['end'], gene_end) - gene_start
            full_coverage[start:end] += row['value']
        
        # 处理负链方向
        if strand == '-':
            full_coverage = full_coverage[::-1]
        
        # 生成exon掩码
        exon_mask = np.zeros(gene_length, dtype=bool)
        for start, end in exons:
            s = start - gene_start
            e = end - gene_start
            exon_mask[s:e] = True
        
        exon_cov = full_coverage[exon_mask]
        intron_cov = full_coverage[~exon_mask]
        
        # --- 计算指标 ---
        # 指标2: 内含子reads比例
        total_reads = exon_cov.sum() + intron_cov.sum()
        intron_ratio = intron_cov.sum() / (total_reads + 1e-10)  # 防止除零
        
        # 存储结果
        results['name'].append(gene_name)
        results['intron_ratio'].append(intron_ratio)
    
    return results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref', required=True, help='Gene annotation BED')
    parser.add_argument('--bed', required=True, help='Read coverage BED')
    parser.add_argument('--out', required=True, help='Output file')
    args = parser.parse_args()
    
    # 读取数据
    ref = pd.read_csv(args.ref, sep='\t', 
                    names=['chr','start','end','name','score','strand'])
    bed = pd.read_csv(args.bed, sep='\t',
                    names=['chr','start','end','value'])
    
    # 多进程处理
    final_results = {k: [] for k in ['name', 'intron_ratio']}
    
    with ProcessPoolExecutor() as executor:
        futures = []
        # 按染色体分组处理
        for chr_id in ref['chr'].unique():
            chr_ref = ref[ref['chr'] == chr_id]
            chr_bed = bed[bed['chr'] == chr_id]
            futures.append(executor.submit(process, chr_ref, chr_bed))
        
        # 收集结果
        for future in futures:
            chr_results = future.result()
            for k in final_results:
                final_results[k].extend(chr_results[k])
    
    # 生成结果文件
    df = pd.DataFrame(final_results)
    df = df[['name', 'intron_ratio']]  # 调整列顺序
    df.to_csv(args.out, sep='\t', index=False)

if __name__ == '__main__':
    main()