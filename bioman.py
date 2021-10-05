import re
import os, os.path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import pandas
import io

def fasta2json(fasta_seq, fasta_file=""):
    with io.BytesIO(fasta_seq.encode()) as f:
        for rec in SeqIO.parse(io.TextIOWrapper(f), "fasta"):
            yield {"fasta_file": os.path.basename(fasta_file), 
                    "seq_id": rec.id, 
                    "description": rec.description, 
                    "seq": str(rec.seq)}

def blast2json(s, outfmt=""):
    if outfmt == "":
        outfmt = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    cols = outfmt.split()
    for line in s.split("\n"):
        if line.strip():
            sp = line.split("\t")
            yield dict(zip(cols, sp))

def json2fasta(json_obj):
    return SeqRecord(
        id = json_obj["seq_id"],
        description = json_obj["description"],
        name = json_obj["seq_id"],
        seq = Seq(json_obj["seq"]))

def json2blast(json_obj):
    return pandas.DataFrame(json_obj)

def parse_orf_faa(fasta_json):
    parser = re.compile(r"lcl\|(ORF\d+)_([^:|\s]+):(\d+):(\d+)\s([\w\s,]+)")
    ORF_name, original_seq_id, start, end, ORF_description = parser.findall(fasta_json["description"]).pop()
    start, end = int(start), int(end)
    ORF_from, ORF_to = min(start, end) + 1, max(start, end) + 1
    fasta_json["ORF_name"] = ORF_name
    fasta_json["original_seq_id"] = original_seq_id
    fasta_json["ORF_from"] = ORF_from
    fasta_json["ORF_to"] = ORF_to
    fasta_json["ORF_isforward"] = end > start
    fasta_json["ORF_ispartial"] = ORF_description.endswith("partial")    
    return fasta_json

def call_orf_finder(query):
    with open(".tmp_in", "w") as f:
        f.write(query)
    ret = ""
    p = subprocess.Popen(["ORFfinder", "-s", "2", "-in", ".tmp_in"], stdout = subprocess.PIPE)
    outs, errs = p.communicate()
    ret = outs.decode()
    return ret

def call_diamond_blast(query, program, database, outfmt="6", timeout=600):
    cmd = ["diamond", program, "-d", database, "-f"]  + outfmt.split()
    ret = ""
    proc = subprocess.Popen(cmd, stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    try:
        outs, errs = proc.communicate(query.encode(), timeout)
    except:
        pass # handle error
    finally:
        proc.kill()
    ret = outs.decode()
    return ret

def find(db, query_dict):
    for doc in db:
        matched = True
        for key in query_dict:
            if key in doc:
                if doc[key] != query_dict[key]:
                    matched = False
                    break
            else:
                matched = False
                break
        if matched:
            yield doc