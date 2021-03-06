import re
import os, os.path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import pandas
import io
import pymongo

from scipy.cluster.hierarchy import linkage, cut_tree
from scipy.spatial.distance import squareform

def fasta2json(fasta_seq, fasta_file=""):
    with io.BytesIO(fasta_seq.encode()) as f:
        for rec in SeqIO.parse(io.TextIOWrapper(f), "fasta"):
            yield {"fasta_file": os.path.basename(fasta_file), 
                    "seq_id": rec.id, 
                    "description": rec.description, 
                    "seq": str(rec.seq)}

def blast2json(s, fmt=""):
    lines = s.split("\n")
    if fmt == "":
        fmt = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    cols = fmt.split()
    for line in lines:
        if line.strip():
            sp = line.split("\t")
            yield dict(zip(cols, sp))

def json2fasta(json_obj):
    return SeqRecord(
            id = json_obj["seq_id"],
            description = json_obj["description"],
            name = json_obj["seq_id"],
            seq = Seq(json_obj["seq"]))

def json2table(json_obj):
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

def parse_msa_dist(file_name, cut_off=90, criteria="single"):
    """ Parse multi-sequences alignment percent identity matrix,
        and assign species (groups) according to similarity.

        Arguments:
            file_name -- the input file should be in csv format, 
                the geneious % identity matirx of MSA is preferred.
            cut_off -- the % identity cut off to group sequences.
            criteria -- see 'scipy.cluster.hierarchy.linkage'

        Return:
            a dict map virus_name to group_label
    """
    df = pandas.read_csv(file_name, index_col=0)
    dist_matrix = squareform(100-df.fillna(100))
    Z = linkage(dist_matrix, criteria)
    group_labels = cut_tree(Z, height=100-cut_off).transpose().flatten()
    
    return dict(zip(df.index, group_labels))

def call_orf_finder(query):
    with open(".tmp_in", "w") as f:
        f.write(query)
    ret = ""
    try:
        proc = subprocess.Popen(["ORFfinder", "-s", "2", "-in", ".tmp_in"], stdout = subprocess.PIPE)
        outs, errs = proc.communicate(300)
        ret = outs.decode()
    except:
        print("Error") # handle error
    finally:
        proc.kill()
    return ret

def call_getorf(query):
    with open(".tmp_in", "w") as f:
        f.write(query)
    ret = ""
    try:
        proc = subprocess.Popen(["getorf", ".tmp_in", ".tmp_out"])
        with open(".tmp_out") as f:
            ret = f.read()
    except:
        print("Error") # handle error
    finally:
        proc.kill()
    return ret

def call_diamond_blast(query, program, database, outfmt="", options="", timeout=600):
    query = query.replace("-", "")
    cmd = ["diamond", program, "-d", database, "-f", "6"] + outfmt.split() + options.split()
    ret = ""
    try:
        proc = subprocess.Popen(cmd, stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        outs, errs = proc.communicate(query.encode(), timeout)
        ret = outs.decode()
    except:
        print("Error") # handle error
    finally:
        proc.kill()
    return ret

def call_cd_hit_est(fasta, threshold):
    with open(".tmp_in", "w") as f:
        f.write(fasta)
    ret = ""
    try:
        proc = subprocess.Popen(["cd-hit-est", "-i", ".tmp_in", "-o", ".tmp_out", "-c", str(threshold)], stdout=subprocess.PIPE)
        outs, errs = proc.communicate()
        with open(".tmp_out") as f:
            ret = f.read()
    except:
        print("Error") # handle error
    finally:
        proc.kill()
    return ret, outs.decode()

def call_cd_hit(fasta, threshold):
    with open(".tmp_in", "w") as f:
        f.write(fasta)
    ret = ""
    try:
        proc = subprocess.Popen(["cd-hit", "-i", ".tmp_in", "-o", ".tmp_out", "-c", str(threshold)], stdout=subprocess.PIPE)
        outs, errs = proc.communicate()
        with open(".tmp_out") as f:
            ret = f.read()
    except:
        print("Error") # handle error
    finally:
        proc.kill()
    return ret, outs.decode()

def call_cap3(fasta, options=""):
    with open(".tmp_in", "w") as f:
        f.write(fasta)
    contigs = ""
    singlets = ""
    contigs_info = ""
    try:
        proc = subprocess.Popen(["cap3", ".tmp_in"] + options.split(), stdout=subprocess.PIPE)
        outs, errs = proc.communicate()
        contigs_info = outs.decode()
        with open(".tmp_in.cap.singlets") as f:
            singlets = f.read()
        with open(".tmp_in.cap.contigs") as f:
            contigs = f.read()
    except:
        print("Error") # handle error
    finally:
        proc.kill()
    return contigs, singlets, contigs_info

def find(db, query_dict):
    for doc in db:
        matched = True
        for key in query_dict:
            if key in doc:
                if type(query_dict[key]).__name__ == "str":
                    criteria = lambda x: x == query_dict[key]
                elif type(query_dict[key]).__name__ == "function":
                    criteria = lambda x: query_dict[key](x)
                elif type(query_dict[key]).__name__ == "Pattern":
                    criteria = lambda x: query_dict[key].match(x)
                
                if not criteria(doc[key]):
                    matched = False
                    break
            else:
                matched = False
                break
        if matched:
            yield doc

def inner_join(left, right, on):
    ret_list = []
    for key in on:
        for ldoc in left:
            for rdoc in find(right, {on[key]: lambda x: x == ldoc[key]}):
                ndoc = {}
                ndoc.update(ldoc)
                ndoc.update(rdoc)
                ret_list.append(ndoc)
    return ret_list

class BioMongoDBClient:
    
    def __init__(self, host="localhost", port=27017, db_name="temp", overwrite=True):
        self.client = pymongo.MongoClient("mongodb://{}:{}/".format(host, port))
        self.db = self.client[db_name]
        dblist = self.client.list_database_names()
        if db_name in dblist:
            if overwrite: # WARNING: We will clear and overwrite it!
                collist = self.db.list_collection_names()
                for col in collist:
                    if col != "system.indexes":
                        self.db[col].drop()

    def insert(self, doc, col_name="temp"):
        col = self.db[col_name]
        ret = None
        if type(doc).__name__ == "list":
            ret = col.insert_many(doc)
        elif type(doc).__name__ == "dict":
            ret = col.insert_one(doc)
        else:
            pass # invalid type
        return ret

    def find(self, query=None, selector=None, col_name="temp"):
        col = self.db[col_name]
        return col.find(query, selector)
    
    def delete(self, query, smode=False, col_name="temp"):
        col = self.db[col_name]
        if smode:
            ret = col.delete_one(query)
        else:
            ret = col.delete_many(query)
        return ret

    def update(self, query, new_doc, smode=False, col_name="temp"):
        col = self.db[col_name]
        if smode:
            ret = col.update_one(query, new_doc)
        else:
            ret = col.update_many(query, new_doc)
        return ret