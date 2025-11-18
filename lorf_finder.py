
#!/usr/bin/env python3
import sys
import json
import argparse
from typing import Dict, List, Tuple, Iterable

GENETIC_CODES: Dict[int, Dict[str, str]] = {}
START_CODONS: Dict[int, set] = {}
STOP_AA = "*"


def _build_code(tt_id: int, aas: str, starts: str, base1: str, base2: str, base3: str) -> None:
    if not (len(aas) == len(base1) == len(base2) == len(base3) == 64):
        raise ValueError("Expected 64-length strings for aas and base rows")
    table = {}
    codons: List[str] = []
    for i in range(64):
        codon = base1[i] + base2[i] + base3[i]
        codons.append(codon)
        table[codon] = aas[i]
    GENETIC_CODES[tt_id] = table
    start_set = set()
    for i in range(64):
        if i < len(starts) and starts[i] == "M":
            start_set.add(codons[i])
    START_CODONS[tt_id] = start_set

# Built-in codes
_build_code(
    1,
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "---M------**--*----M---------------M---------------M------------",
    "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
    "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
    "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
)

_build_code(
    2,
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    "----------**--------------------MMMM----------**---M------------",
    "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
    "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
    "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
)

_build_code(
    11,
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "---M------**--*----M------------MMMM---------------M------------",
    "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
    "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
    "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
)


def load_custom_code_from_json(path: str, tt_id: int = 100) -> int:
    with open(path, "r") as fh:
        data = json.load(fh)
    # Expect either mapping of codon->aa, and optional starts list
    mapping = data.get("mapping") or data
    if not isinstance(mapping, dict):
        raise ValueError("JSON must contain a mapping object of codon to amino acid")
    code_table: Dict[str, str] = {}
    for codon, aa in mapping.items():
        c = codon.upper().replace("U", "T")
        if len(c) != 3:
            raise ValueError("Invalid codon: " + codon)
        code_table[c] = str(aa)[0]
    GENETIC_CODES[tt_id] = code_table
    starts = data.get("starts", ["ATG"]) if isinstance(data, dict) else ["ATG"]
    START_CODONS[tt_id] = set([s.upper().replace("U", "T") for s in starts])
    return tt_id


def load_custom_code_from_rows(aas: str, base1: str, base2: str, base3: str, starts: str = "") -> int:
    custom_id = 100
    if not (len(aas) == len(base1) == len(base2) == len(base3) == 64):
        raise ValueError("Inline aas/base rows must be 64 characters each")
    if not starts:
        starts = "".ljust(64, "-")
    _build_code(custom_id, aas, starts, base1, base2, base3)
    return custom_id


def read_fasta(path: str) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    header = None
    parts: List[str] = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = ("".join(parts)).upper().replace(" ", "").replace("	", "")
                header = line[1:].strip()
                parts = []
            else:
                parts.append(line)
        if header is not None:
            seqs[header] = ("".join(parts)).upper().replace(" ", "").replace("	", "")
    return seqs


def write_fasta(records: Iterable[Tuple[str, str]], out_path: str, width: int = 60) -> None:
    with open(out_path, "w") as fh:
        for header, seq in records:
            fh.write(">" + header + "
")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i+width] + "
")


def revcomp(seq: str) -> str:
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join([comp.get(b, "N") for b in seq.upper()[::-1]])


def translate(seq: str, table_id: int = 1) -> str:
    code = GENETIC_CODES.get(table_id)
    if code is None:
        raise ValueError("Unknown translation table id: " + str(table_id))
    prot = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = code.get(codon, "X")
        prot.append(aa)
    return "".join(prot)


def find_lorfs(seq: str, table_id: int, min_aa: int = 100, require_start: bool = True) -> List[Tuple[int, int, str, str, str]]:
    stops = set([c for c, a in GENETIC_CODES[table_id].items() if a == "*"])
    starts = START_CODONS.get(table_id, set(["ATG"]))

    def scan_frames(dna: str, strand_label: str) -> List[Tuple[int, int, str, str, str]]:
        results: List[Tuple[int, int, str, str, str]] = []
        n = len(dna)
        for frame in range(3):
            i = frame
            while i + 2 < n:
                codon = dna[i:i+3]
                if (not require_start) or (codon in starts):
                    j = i
                    while j + 2 < n:
                        c = dna[j:j+3]
                        if c in stops and j > i:
                            break
                        j += 3
                    end_index = j
                    nt_len = end_index - i
                    if nt_len >= 3:
                        aa = translate(dna[i:end_index], table_id)
                        if len(aa) >= min_aa:
                            if strand_label == "+":
                                start_nt = i + 1
                                end_nt = end_index
                            else:
                                L = len(seq)
                                start_nt = L - end_index + 1
                                end_nt = L - i
                            results.append((start_nt, end_nt, strand_label, str(frame), aa))
                    i = j + 3
                else:
                    i += 3
        return results

    forward = scan_frames(seq, "+")
    reverse = scan_frames(revcomp(seq), "-")
    return forward + reverse


def lorfs_from_fasta(in_fasta: str, out_fasta: str, table_id: int = 1, min_aa: int = 100, require_start: bool = True) -> List[Tuple[str, str]]:
    seqs = read_fasta(in_fasta)
    outputs: List[Tuple[str, str]] = []
    for header, dna in seqs.items():
        lorfs = find_lorfs(dna, table_id=table_id, min_aa=min_aa, require_start=require_start)
        lorfs.sort(key=lambda x: len(x[4]), reverse=True)
        for idx, (s, e, strand, frame, aa) in enumerate(lorfs):
            new_header = header + "|lorf_" + str(idx + 1) + "|" + "strand=" + strand + "|frame=" + frame + "|start=" + str(s) + "|end=" + str(e) + "|lenAA=" + str(len(aa)) + "|table=" + str(table_id)
            outputs.append((new_header, aa))
    write_fasta(outputs, out_fasta)
    return outputs


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Find Long ORFs and output translated protein FASTA")
    p.add_argument("input_fasta", help="Input nucleotide FASTA file")
    p.add_argument("output_fasta", help="Output amino acid FASTA file")
    p.add_argument("--table", type=int, default=1, help="Built-in NCBI translation table id (default 1)")
    p.add_argument("--table_file", type=str, default=None, help="Path to JSON file defining custom translation table")
    p.add_argument("--aas", type=str, default=None, help="64-char amino acid string for inline custom table")
    p.add_argument("--base1", type=str, default=None, help="64-char base1 row for inline custom table")
    p.add_argument("--base2", type=str, default=None, help="64-char base2 row for inline custom table")
    p.add_argument("--base3", type=str, default=None, help="64-char base3 row for inline custom table")
    p.add_argument("--starts", type=str, default=None, help="Optional 64-char starts mask with M at initiator positions")
    p.add_argument("--min_aa", type=int, default=100, help="Minimum amino acid length for LORFs (default 100)")
    p.add_argument("--no_start", action="store_true", help="Do not require a start codon to begin ORFs")
    return p


def main(argv: List[str]) -> None:
    ap = build_argparser()
    args = ap.parse_args(argv)
    table_id = args.table

    # Custom table precedence: table_file > inline aas/base rows > built-in
    if args.table_file:
        table_id = load_custom_code_from_json(args.table_file, tt_id=100)
    elif args.aas and args.base1 and args.base2 and args.base3:
        table_id = load_custom_code_from_rows(
            args.aas, args.base1, args.base2, args.base3, starts=(args.starts or "")
        )

    req_start = not args.no_start
    lorfs_from_fasta(args.input_fasta, args.output_fasta, table_id=table_id, min_aa=args.min_aa, require_start=req_start)

if __name__ == "__main__":
    main(sys.argv[1:])
