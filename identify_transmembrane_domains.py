# Hydrophobic segment finder (Jupyter Notebook version)

def find_hydrophobic_segments(sequence, window_size=15, threshold=0.7):
    """
    Identify hydrophobic segments in an amino acid sequence.

    Parameters:
        sequence (str): The amino acid sequence (no spaces/newlines)
        window_size (int): Minimum length of hydrophobic window
        threshold (float): Fraction of hydrophobic residues required (0–1)

    Returns:
        list of dicts with keys: start, end, length, fraction, segment
    """

    hydrophobic_residues = set("AVILMFWYC")
    seq = sequence.upper().replace("\n", "").replace(" ", "")
    n = len(seq)
    hydrophobic_positions = [
        1 if aa in hydrophobic_residues else 0 for aa in seq
    ]

    segments = []
    i = 0
    while i <= n - window_size:
        window = hydrophobic_positions[i : i + window_size]
        if sum(window) / window_size >= threshold:
            # expand to include neighboring hydrophobic residues
            start = i
            end = i + window_size
            while end < n and seq[end] in hydrophobic_residues:
                end += 1
            segments.append((start, end))
            i = end  # skip ahead past this segment
        else:
            i += 1

    # Merge overlapping or adjacent segments
    merged = []
    for s, e in segments:
        if not merged or s > merged[-1][1]:
            merged.append([s, e])
        else:
            merged[-1][1] = max(merged[-1][1], e)

    # Summarize results
    results = []
    for s, e in merged:
        seg = seq[s:e]
        hydrophobic_fraction = sum(aa in hydrophobic_residues for aa in seg) / len(seg)
        results.append({
            "start": s + 1,  # 1-based
            "end": e,
            "length": len(seg),
            "fraction": round(hydrophobic_fraction, 2),
            "segment": seg
        })
    return results


def show_hydrophobic_segments(sequence, window_size=15, threshold=0.7):
    """
    Print formatted hydrophobic segments in a notebook-friendly way.
    """
    results = find_hydrophobic_segments(sequence, window_size, threshold)
    print(f"\nDetected {len(results)} hydrophobic segments:\n")
    for seg in results:
        print(
            f"{seg['start']:>5}-{seg['end']:<5} | "
            f"len={seg['length']:<3} | "
            f"hydro%={seg['fraction']*100:>5.1f}% | "
            f"{seg['segment']}"
        )

# Paste your amino acid sequence as a single string below
sequence = ""

# Run it
show_hydrophobic_segments(sequence)

