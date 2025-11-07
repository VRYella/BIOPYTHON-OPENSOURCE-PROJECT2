# == Imports and Streamlit setup ==
import streamlit as st; from Bio import SeqIO; from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio import pairwise2; from Bio.pairwise2 import format_alignment
import io; import pandas as pd; import matplotlib.pyplot as plt; import base64; import textwrap; import traceback; import re
from Bio import Entrez  # For NCBI sequence fetch

# == Page configuration and main title ==
st.set_page_config(page_title="Biopython Toolkit", layout="wide")
st.title("Biopython Toolkit — Streamlit App")
st.markdown("""Elegant toolbox for common sequence tasks using **Biopython**.<br>
Upload/paste FASTA, fetch sequences from NCBI, inspect, compute GC, translate, ORF, align, search motifs, visualize codon/GC plots.
""", unsafe_allow_html=True)
st.sidebar.header("Quick actions")
st.sidebar.markdown("Upload, paste FASTA or fetch NCBI sequence. Select tools below.", unsafe_allow_html=True)
st.sidebar.markdown("---")

# == Helper: Parse FASTA/pasted text/NCBI fetch ==
def parse_fasta_bytes(uploaded_file): # Parse uploaded file
    try:
        if uploaded_file is None: return []
        content = uploaded_file.read().decode('utf-8'); handle = io.StringIO(content)
        return list(SeqIO.parse(handle, "fasta"))
    except Exception as e: st.error("Failed to parse FASTA. " + str(e)); return []

def parse_fasta_text(text): # Parse paste box
    try:
        if not text.strip(): return []
        handle = io.StringIO(text)
        return list(SeqIO.parse(handle, "fasta"))
    except Exception as e: st.error(f"Could not parse pasted sequence: {e}"); return []

def fetch_ncbi(uid, db="nucleotide"): # Fetch by NCBI accession
    try:
        Entrez.email = "user@example.com"
        handle = Entrez.efetch(db=db, id=uid, rettype="fasta", retmode="text")
        return list(SeqIO.parse(handle, "fasta"))
    except Exception as e: st.error(f"NCBI fetch failed: {e}"); return []

def download_link(data: bytes, fname: str, text: str): # File download link
    b64 = base64.b64encode(data).decode(); return f'<a href="data:file/txt;base64,{b64}" download="{fname}">{text}</a>'

# == Data Input Panel ==
with st.expander("Input Sequences (FASTA, NCBI, Paste)", expanded=True):
    tab_uploaded, tab_paste, tab_ncbi = st.tabs(["Upload FASTA", "Paste FASTA", "Fetch from NCBI"])
    with tab_uploaded: uploaded = st.file_uploader("Upload a FASTA file", type=["fa","fasta","txt"], label_visibility="collapsed")
    with tab_paste: pasted = st.text_area("Paste FASTA data here", height=120, label_visibility="collapsed")
    with tab_ncbi:
        ncbi_id = st.text_input("Enter NCBI accession(s) (comma separated)", key="ncbi_input")
        ncbi_btn = st.button("Fetch from NCBI", key="ncbi_btn")
        ncbi_seqs = []
        if ncbi_btn and ncbi_id.strip():
            ids = [x.strip() for x in re.split(r",|;|\s", ncbi_id) if x.strip()]
            for uid in ids: ncbi_seqs += fetch_ncbi(uid)
# == Example FASTA display ==
with st.expander("Example FASTA"):
    st.code(">seq1_example\nATGCGTACGTTAGCGTAGCTAGCTAGCGTAGCTAGCTGACTGATCGATCGATCGTAGCTAG\n>seq2_example\nATGAAATTTGGGCCCTTTAAACCCGGGATGCTAGCTAGCTAA")

# == Gather all sequences ==
records = []
if uploaded: records += parse_fasta_bytes(uploaded)
if pasted and not uploaded: records += parse_fasta_text(pasted)
if ncbi_btn: records += ncbi_seqs

# == Toolbox control panel ==
st.header("Analysis Tools")
col1, col2 = st.columns(2)
with col1:
    show_upload = st.checkbox("Show sequences", True)
    show_summary = st.checkbox("Summary (length, GC%)", True)
    show_revtrans = st.checkbox("Reverse compl./translate", True)
with col2:
    show_orf = st.checkbox("ORF/six-frame finder", True)
    show_align = st.checkbox("Pairwise alignment", True)
    show_motif = st.checkbox("Motif (regex) search", True)
    show_codon = st.checkbox("Codon/GC plots", True)

# == Uploaded / Pasted / NCBI sequence display block ==
if show_upload:
    st.subheader("Sequences loaded")
    if not records: st.info("No sequences loaded. Add via upload, paste, or NCBI.")
    else:
        for rec in records:
            with st.expander(f"{rec.id} | {len(rec.seq)} bp"):
                st.write(f"<b>Description:</b> {rec.description}", unsafe_allow_html=True)
                st.code(str(rec.seq[:1000]) + ("..." if len(rec.seq)>1000 else ""))
                st.markdown(download_link(f">{rec.id}\n{str(rec.seq)}\n".encode(), f"{rec.id}.fasta", "Download FASTA"), unsafe_allow_html=True)

# == Sequence summary block: length, GC ==
if show_summary:
    st.subheader("Sequence summary")
    if not records: st.info("Add sequence above for summary.")
    else:
        rows = [{"ID": r.id, "Length": len(r.seq), "GC%": round(gc_fraction(r.seq)*100,3)} for r in records]
        df = pd.DataFrame(rows); st.dataframe(df, use_container_width=True); st.markdown(download_link(df.to_csv(index=False).encode(), "summary.csv", "Download CSV"), unsafe_allow_html=True)

# == Reverse Complement & Translation UI ==
if show_revtrans:
    st.subheader("Reverse Complement & Translation")
    seq_choice = st.selectbox("Select sequence", [r.id for r in records], key="rev_choice")
    table = st.checkbox("Show translation details", False)
    if seq_choice:
        rec = next(r for r in records if r.id==seq_choice); seq = rec.seq
        try:
            st.markdown("**Reverse Complement**"); st.code(str(seq.reverse_complement()))
            st.markdown("**Translation (std Table)**"); prot = seq.translate(to_stop=False)
            st.code(str(prot)); st.markdown("**Protein MW (Da)**");
            try: st.write(round(molecular_weight(prot),3))
            except Exception: st.write("Error computing MW.")
            if table: st.write("Translation (first 200 aa):"); st.code(str(prot[:200]))
        except Exception as e: st.error("Translation/rev.comp error: " + str(e))

# == Six-frame ORF finder ==
if show_orf:
    st.subheader("Six-frame ORF/Translation")
    seq_choice = st.selectbox("Sequence for ORF", [r.id for r in records], key="orf_choice2")
    min_orf_len = st.number_input("Min ORF length (aa)", 10, 10000, 30)
    search_start_meth = st.selectbox("ORF method", ["Start=M, end=Stop", "Any open w/o internal stop"], key="orf_method")
    if seq_choice:
        rec = next(r for r in records if r.id==seq_choice); seq = rec.seq; frames = []
        def find_orfs_simple(s, strand, frame):
            prot = s.translate(to_stop=False); orfs = []
            if search_start_meth=="Start=M, end=Stop":
                for i, aa in enumerate(prot):
                    if aa == "M":
                        for j in range(i+1, len(prot)):
                            if prot[j] == "*":
                                length = j - i + 1
                                if length >= min_orf_len: 
                                    orfs.append({"strand": strand, "frame": frame, "aa_start": i, "aa_end": j, "length_aa": length,
                                                 "protein": str(prot[i:j+1]), "nuc_start": (i*3)+frame, "nuc_end": (j*3)+frame+3})
                                break
            else:
                start = None
                for i, aa in enumerate(prot):
                    if aa != "*" and start is None: start = i
                    if aa == "*" and start is not None:
                        length = i - start
                        if length >= min_orf_len:
                            orfs.append({"strand": strand, "frame": frame, "aa_start": start, "aa_end": i-1, "length_aa": length,
                                "protein": str(prot[start:i]), "nuc_start": (start*3)+frame, "nuc_end": (i*3)+frame})
                        start = None
            return orfs
        seq_str = seq
        for strand, s in [(+1, seq_str), (-1, seq_str.reverse_complement())]:
            for frame in range(3): s_frame = s[frame:]; frames += find_orfs_simple(s_frame, strand, frame)
        if frames:
            df_orfs = pd.DataFrame(frames).sort_values("length_aa", ascending=False)
            st.dataframe(df_orfs, use_container_width=True)
            st.markdown(download_link(df_orfs.to_csv(index=False).encode(), "orfs.csv", "Download ORFs CSV"), unsafe_allow_html=True)
        else: st.write("No ORFs found. Lower min length.")

# == Pairwise alignment tool ==
if show_align:
    st.subheader("Pairwise alignment (global/local)")
    seqs = [r.id for r in records]
    if len(seqs) < 2: st.info("Load >=2 sequences for alignment.")
    else:
        s1 = st.selectbox("Sequence 1", seqs, key="s1"); s2 = st.selectbox("Sequence 2", seqs, key="s2")
        method = st.radio("Alignment type", ["global", "local"], key="align_method")
        rec1 = next(r for r in records if r.id==s1); rec2 = next(r for r in records if r.id==s2)
        if st.button("Align!", key="align_btn"):
            try:
                alns = pairwise2.align.globalxx(str(rec1.seq), str(rec2.seq), one_alignment_only=True) if method=="global" else pairwise2.align.localxx(str(rec1.seq), str(rec2.seq), one_alignment_only=True)
                if alns:
                    aln = alns[0]
                    st.text(format_alignment(*aln))
                    aln_text = format_alignment(*aln)
                    st.markdown(download_link(aln_text.encode(), f"alignment_{s1}_vs_{s2}.txt", "Download result"), unsafe_allow_html=True)
                else: st.write("No alignment found.")
            except Exception as e: st.error("Alignment failed: " + str(e)); st.text(traceback.format_exc())

# == Motif (restriction/regex) search ==
if show_motif:
    st.subheader("Motif / Restriction regex")
    motif = st.text_input("Motif (DNA regex, e.g. GAATTC)", key="motif_txt")
    seq_choice = st.selectbox("Seq for motif", [r.id for r in records], key="motif_seq")
    if motif and seq_choice:
        rec = next(r for r in records if r.id==seq_choice); seq = str(rec.seq)
        try:
            matches = [(m.start()+1, m.group(0)) for m in re.finditer(motif, seq)]
            st.write(f"Occurrences: {len(matches)}")
            if matches:
                dfm = pd.DataFrame(matches, columns=["start","motif"])
                st.dataframe(dfm, use_container_width=True)
                st.markdown(download_link(dfm.to_csv(index=False).encode(), "motif_hits.csv", "Download hits"), unsafe_allow_html=True)
        except re.error as e: st.error("Invalid regex: " + str(e))

# == Codon usage and GC distribution plotting ==
if show_codon:
    st.subheader("Codon usage & GC% plot")
    seq_choice = st.selectbox("Seq for analysis", [r.id for r in records], key="codon_seq")
    if seq_choice:
        rec = next(r for r in records if r.id==seq_choice); seq = str(rec.seq)
        codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]; codon_counts = {}
        for c in codons: codon_counts[c] = codon_counts.get(c,0)+1 if len(c)==3 else codon_counts.get(c,0)
        df_codon = pd.DataFrame(list(codon_counts.items()), columns=["Codon","Count"]).sort_values("Count", ascending=False)
        st.dataframe(df_codon, use_container_width=True)
        st.markdown(download_link(df_codon.to_csv(index=False).encode(), "codon_usage.csv", "Download codons"), unsafe_allow_html=True)
        fig, ax = plt.subplots(figsize=(7,3)); ax.bar(df_codon['Codon'], df_codon['Count'])
        ax.set_xlabel("Codon"); ax.set_ylabel("Count"); ax.set_title("Codon usage (frame 0)")
        plt.xticks(rotation=90); st.pyplot(fig)
        window = st.slider("GC sliding window (bp)", 10, 200, 50); seq_upper = seq.upper(); gc_values, positions = [], []
        for i in range(0, max(1, len(seq_upper)-window+1), window//2):
            win = seq_upper[i:i+window]; positions.append(i+1); gc_values.append(gc_fraction(win)*100)
        fig2, ax2 = plt.subplots(figsize=(7,2.5)); ax2.plot(positions, gc_values)
        ax2.set_xlabel("Pos (bp)"); ax2.set_ylabel("GC%"); ax2.set_title("Sliding-window GC%"); st.pyplot(fig2)

# == Sidebar information and developer export block ==
st.sidebar.markdown("---"); st.sidebar.markdown("Made with Biopython & Streamlit. No network except NCBI fetch.")
st.sidebar.markdown("Download full app package below.")

# == Developer export tool ==
import zipfile, os; from io import BytesIO
def make_zip_bytes(files):
    buf = BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for fpath, arcname in files:
            zf.writestr(arcname, open(fpath, "rb").read())
    return buf.getvalue()

if st.sidebar.button("Create app zip (download)"):
    try:
        files = [
            ("/mnt/data/streamlit_biopython_app.py", "streamlit_biopython_app.py"),
            ("/mnt/data/requirements.txt", "requirements.txt"),
            ("/mnt/data/README_streamlit_app.txt", "README_streamlit_app.txt")
        ]
        zip_bytes = make_zip_bytes(files)
        st.markdown(download_link(zip_bytes, "biopython_streamlit_app.zip", "Download ZIP"), unsafe_allow_html=True)
    except Exception as e: st.error("Zip creation failed: " + str(e))
# End of application
