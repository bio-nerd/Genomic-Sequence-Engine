cat << 'EOF' > Fetch.py
import requests
import json
import re

def parse_natural_query(user_input):
    """Simple NLP/Intent parser for genomic sequence requests."""
    user_input_lower = user_input.lower().strip()
    
    # 1. Extract Gene Symbol (Capitalized word or standard gene pattern)
    genes = ["brca1", "brca2", "egfr", "tp53", "kras", "myc", "pten", "akt1", "alk", "braf"]
    found_gene = None
    for g in genes:
        if g in user_input_lower:
            found_gene = g.upper()
            break
            
    if not found_gene:
        # Fallback regex search for gene-like alphanumeric tokens
        match = re.search(r'\b[A-Za-z0-9]{3,8}\b', user_input)
        found_gene = match.group(0).upper() if match else "TP53"

    # 2. Extract Promoter / Upstream Window Intent
    upstream_bp = 0
    if "promoter" in user_input_lower or "upstream" in user_input_lower:
        upstream_match = re.search(r'(\d+)\s*(?:bp|base)?', user_input_lower)
        upstream_bp = int(upstream_match.group(1)) if upstream_match else 500

    return found_gene, upstream_bp

def calculate_gc_and_cpg(sequence, window_size=200):
    """Analyzes overall GC content and scans for CpG Islands."""
    gc_count = sequence.count('G') + sequence.count('C')
    gc_content = (gc_count / len(sequence)) * 100 if sequence else 0
    
    # CpG Island Detection: Window with GC > 50% & Obs/Exp ratio > 0.6
    cpg_islands = []
    for i in range(0, len(sequence) - window_size, 50):
        subseq = sequence[i:i+window_size]
        g = subseq.count('G')
        c = subseq.count('C')
        cg = subseq.count('CG')
        
        sub_gc = ((g + c) / window_size) * 100
        expected_cg = (c * g) / window_size
        obs_exp = (cg / expected_cg) if expected_cg > 0 else 0
        
        if sub_gc >= 50 and obs_exp >= 0.6:
            cpg_islands.append((i, i + window_size, sub_gc, obs_exp))
            
    return gc_content, cpg_islands

def find_orfs(sequence, min_len_aa=30):
    """Scans 3 forward reading frames for Open Reading Frames (ORFs)."""
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []
    
    for frame in range(3):
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon == start_codon:
                for j in range(i + 3, len(sequence) - 2, 3):
                    stop = sequence[j:j+3]
                    if stop in stop_codons:
                        length_bp = j + 3 - i
                        length_aa = length_bp // 3
                        if length_aa >= min_len_aa:
                            orfs.append({
                                'frame': frame + 1,
                                'start': i + 1,
                                'end': j + 3,
                                'length_aa': length_aa
                            })
                        break
    return orfs

def fetch_and_analyze(query):
    gene_symbol, upstream_bp = parse_natural_query(query)
    print(f"\n[*] Natural Query Parsed -> Gene: {gene_symbol} | Upstream: {upstream_bp} bp")
    print(f"[*] Querying Ensembl REST API...")
    
    lookup_url = f"https://rest.ensembl.org/lookup/symbol/human/{gene_symbol}?content-type=application/json"
    r = requests.get(lookup_url)
    
    if not r.ok:
        print(f"[!] Error: Could not resolve gene '{gene_symbol}'.")
        return
        
    gene_data = r.json()
    gene_id = gene_data['id']
    chrom = gene_data['seq_region_name']
    start, end, strand = gene_data['start'], gene_data['end'], gene_data['strand']
    
    if upstream_bp > 0:
        fetch_start = max(1, start - upstream_bp) if strand == 1 else start
        fetch_end = end if strand == 1 else end + upstream_bp
    else:
        fetch_start, fetch_end = start, end

    seq_url = f"https://rest.ensembl.org/sequence/region/human/{chrom}:{fetch_start}..{fetch_end}:{strand}?content-type=application/json"
    seq_r = requests.get(seq_url)
    
    if not seq_r.ok:
        print("[!] Error fetching sequence.")
        return
        
    sequence = seq_r.json()['seq'].upper()
    
    # Intelligence Layer Calculations
    gc_pct, cpg_islands = calculate_gc_and_cpg(sequence)
    orfs = find_orfs(sequence)
    
    # Save standard FASTA with metadata header
    filename = f"{gene_symbol}_annotated.fasta"
    fasta_header = f">{gene_symbol} | {gene_id} | chr{chrom}:{fetch_start}-{fetch_end} | GC={gc_pct:.2f}% | CpG_Islands={len(cpg_islands)} | Potential_ORFs={len(orfs)}\n"
    
    with open(filename, "w") as f:
        f.write(fasta_header)
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")
            
    # Print Intelligence Summary
    print("\n" + "="*50)
    print(f" SUCCESS: Generated '{filename}' ({len(sequence):,} bp)")
    print("="*50)
    print(f"▶ Overall GC Content: {gc_pct:.2f}%")
    print(f"▶ CpG Islands Detected: {len(cpg_islands)} regions")
    if cpg_islands:
        print(f"   ↳ First island at bp {cpg_islands[0][0]}-{cpg_islands[0][1]} (GC: {cpg_islands[0][2]:.1f}%, Obs/Exp: {cpg_islands[0][3]:.2f})")
    print(f"▶ Open Reading Frames (>=30 aa): {len(orfs)} found")
    if orfs:
        top_orf = max(orfs, key=lambda x: x['length_aa'])
        print(f"   ↳ Longest ORF: Frame +{top_orf['frame']} ({top_orf['start']}-{top_orf['end']} bp, {top_orf['length_aa']} aa)")
    print("="*50 + "\n")

if __name__ == "__main__":
    prompt = input("Ask the Genomic Assistant (e.g., 'Fetch BRCA1 promoter 500bp'): ")
    fetch_and_analyze(prompt)
EOF
)
