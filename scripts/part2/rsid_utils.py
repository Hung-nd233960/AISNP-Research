"""
rsid_utils.py - Utilities for converting rsID lists to BED files

This module provides functions to:
1. Load rsID lists from CSV files
2. Query genomic coordinates for rsIDs (via Ensembl REST API or local dbSNP)
3. Generate BED files with chr, start, end, rsID columns
"""

import pandas as pd
import numpy as np
import requests
import time
from pathlib import Path
from typing import List, Dict, Optional, Union
from tqdm import tqdm
import json


def load_rsid_list(file_path: str, column_name: Optional[str] = None) -> List[str]:
    """
    Load rsID list from a CSV or text file.
    
    Args:
        file_path: Path to file containing rsIDs (one per line or in a column)
        column_name: If CSV with header, specify the column name containing rsIDs
        
    Returns:
        List of rsID strings (e.g., ['rs123', 'rs456', ...])
    """
    file_path = Path(file_path)
    
    if file_path.suffix.lower() == '.csv':
        df = pd.read_csv(file_path)
        if column_name:
            rsids = df[column_name].dropna().astype(str).tolist()
        else:
            # Try common column names
            for col in ['rsid', 'rs_id', 'snp', 'snp_id', 'id', 'ID', 'rsID', 'RSID']:
                if col in df.columns:
                    rsids = df[col].dropna().astype(str).tolist()
                    break
            else:
                # Assume first column
                rsids = df.iloc[:, 0].dropna().astype(str).tolist()
    else:
        # Plain text file, one rsID per line
        with open(file_path, 'r') as f:
            rsids = [line.strip() for line in f if line.strip()]
    
    # Normalize rsIDs (ensure 'rs' prefix)
    normalized = []
    for rsid in rsids:
        rsid = rsid.strip()
        if rsid.isdigit():
            rsid = f"rs{rsid}"
        normalized.append(rsid)
    
    print(f"Loaded {len(normalized)} rsIDs from {file_path}")
    return normalized


def query_rsid_coordinates(rsid: str, assembly: str = 'GRCh38') -> Optional[Dict]:
    """
    Query Ensembl REST API for rsID genomic coordinates.
    
    Args:
        rsid: rsID string (e.g., 'rs123')
        assembly: Genome assembly ('GRCh38' or 'GRCh37')
        
    Returns:
        Dict with 'chr', 'start', 'end', 'rsid' or None if not found
    """
    # Use GRCh37 server for hg19, main server for GRCh38
    if assembly == 'GRCh37':
        server = "https://grch37.rest.ensembl.org"
    else:
        server = "https://rest.ensembl.org"
    
    endpoint = f"/variation/human/{rsid}?content-type=application/json"
    
    try:
        response = requests.get(server + endpoint, timeout=10)
        if response.status_code == 200:
            data = response.json()
            
            # Extract mappings
            if 'mappings' in data and len(data['mappings']) > 0:
                mapping = data['mappings'][0]
                chrom = mapping.get('seq_region_name', '')
                start = mapping.get('start', 0)
                end = mapping.get('end', 0)
                
                # Normalize chromosome name
                if not chrom.startswith('chr'):
                    chrom = f"chr{chrom}"
                
                return {
                    'chr': chrom,
                    'start': start,
                    'end': end,
                    'rsid': rsid
                }
        return None
    except Exception as e:
        print(f"Error querying {rsid}: {e}")
        return None


def batch_rsid_lookup(rsids: List[str], assembly: str = 'GRCh38', 
                      delay: float = 0.1, use_cache: bool = True,
                      cache_file: Optional[str] = None) -> pd.DataFrame:
    """
    Batch lookup of rsID coordinates with rate limiting and caching.
    
    Args:
        rsids: List of rsIDs
        assembly: Genome assembly
        delay: Delay between API calls (seconds)
        use_cache: Whether to use/update cache file
        cache_file: Path to cache file (CSV)
        
    Returns:
        DataFrame with chr, start, end, rsid columns
    """
    results = []
    cache = {}
    
    # Load cache if exists
    if use_cache and cache_file and Path(cache_file).exists():
        cache_df = pd.read_csv(cache_file)
        cache = {row['rsid']: row.to_dict() for _, row in cache_df.iterrows()}
        print(f"Loaded {len(cache)} cached entries")
    
    # Query missing rsIDs
    for rsid in tqdm(rsids, desc="Looking up rsIDs"):
        if rsid in cache:
            results.append(cache[rsid])
        else:
            coord = query_rsid_coordinates(rsid, assembly)
            if coord:
                results.append(coord)
                cache[rsid] = coord
            else:
                results.append({
                    'chr': None, 'start': None, 'end': None, 'rsid': rsid
                })
            time.sleep(delay)  # Rate limiting
    
    df = pd.DataFrame(results)
    
    # Save cache
    if use_cache and cache_file:
        cache_df = pd.DataFrame(list(cache.values()))
        cache_df.to_csv(cache_file, index=False)
        print(f"Cache saved to {cache_file}")
    
    return df


def rsid_to_bed(rsids: Union[str, List[str]], output_path: str,
                assembly: str = 'GRCh38', source_name: str = 'unknown') -> str:
    """
    Convert rsID list to BED file format.
    
    Args:
        rsids: Path to rsID file or list of rsIDs
        output_path: Output BED file path
        assembly: Genome assembly for coordinate lookup
        source_name: Name of the source (for cache naming)
        
    Returns:
        Path to generated BED file
    """
    # Load rsIDs if path provided
    if isinstance(rsids, str):
        rsid_list = load_rsid_list(rsids)
    else:
        rsid_list = rsids
    
    # Cache file
    cache_file = str(Path(output_path).parent / f"{source_name}_rsid_cache.csv")
    
    # Lookup coordinates
    print(f"Looking up {len(rsid_list)} rsIDs...")
    coords_df = batch_rsid_lookup(rsid_list, assembly, cache_file=cache_file)
    
    # Filter successful lookups
    valid_df = coords_df.dropna(subset=['chr', 'start', 'end'])
    print(f"Successfully resolved {len(valid_df)}/{len(rsid_list)} rsIDs")
    
    # Create BED format (0-based start)
    bed_df = pd.DataFrame({
        'chr': valid_df['chr'],
        'start': valid_df['start'].astype(int) - 1,  # BED is 0-based
        'end': valid_df['end'].astype(int),
        'rsid': valid_df['rsid']
    })
    
    # Sort by chromosome and position
    bed_df = bed_df.sort_values(['chr', 'start'])
    
    # Save BED file (no header, tab-separated)
    bed_df.to_csv(output_path, sep='\t', index=False, header=False)
    print(f"BED file saved: {output_path}")
    
    # Also save with header for reference
    bed_with_header = output_path.replace('.bed', '_annotated.csv')
    bed_df.to_csv(bed_with_header, index=False)
    print(f"Annotated CSV saved: {bed_with_header}")
    
    return output_path


def load_local_dbsnp(dbsnp_path: str) -> pd.DataFrame:
    """
    Load local dbSNP VCF for faster lookups (optional).
    
    Args:
        dbsnp_path: Path to dbSNP VCF file
        
    Returns:
        DataFrame indexed by rsID with coordinates
    """
    print(f"Loading local dbSNP from {dbsnp_path}...")
    # This would parse a local dbSNP VCF
    # For now, return empty - users can implement if they have local dbSNP
    raise NotImplementedError("Local dbSNP loading not implemented. Use API lookup.")


if __name__ == "__main__":
    # Example usage
    import sys
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        output_file = sys.argv[2] if len(sys.argv) > 2 else "output.bed"
        rsid_to_bed(input_file, output_file)
