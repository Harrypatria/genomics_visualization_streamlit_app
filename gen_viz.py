import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm
import matplotlib.colors as mcolors
import io
import base64
from Bio import SeqIO, Seq, SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pygenomeviz as pgv
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import tempfile
import os
import random
import gffutils
from dna_features_viewer import BiopythonTranslator, GraphicFeature, GraphicRecord
import networkx as nx
import requests

def load_genbank_file(gb_file):
    """
    Load GenBank file using BioPython
    
    Parameters:
    -----------
    gb_file : file
        GenBank file object
        
    Returns:
    --------
    list
        List of SeqRecord objects
    """
    return list(SeqIO.parse(gb_file, "genbank"))

# Set page configuration
st.set_page_config(
    page_title="Genomic Visualization Hub",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
    .main {
        background-color: #f8f9fa;
    }
    .stApp {
        max-width: 1200px;
        margin: 0 auto;
    }
    h1, h2, h3 {
        color: #2c3e50;
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 24px;
    }
    .stTabs [data-baseweb="tab"] {
        height: 50px;
        white-space: pre-wrap;
        background-color: #f1f3f6;
        border-radius: 4px 4px 0px 0px;
        gap: 1px;
        padding-top: 10px;
        padding-bottom: 10px;
    }
    .stTabs [aria-selected="true"] {
        background-color: #4e89ae;
        color: white;
    }
    .highlight {
        background-color: #f0f7fb;
        border-left: 5px solid #3498db;
        padding: 10px;
        margin: 10px 0;
    }
    .info-box {
        background-color: #e8f4f8;
        border-radius: 5px;
        padding: 15px;
        margin: 10px 0;
        border-left: 5px solid #3498db;
    }
    .warning-box {
        background-color: #fff3e0;
        border-radius: 5px;
        padding: 15px;
        margin: 10px 0;
        border-left: 5px solid #ff9800;
    }
    .feature-card {
        padding: 20px;
        border-radius: 5px;
        background-color: white;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        margin-bottom: 20px;
    }
</style>
""", unsafe_allow_html=True)

# Utility functions
def create_download_link(file_content, file_name, file_label="Download File"):
    """Generate a download link for any file content"""
    b64 = base64.b64encode(file_content).decode()
    href = f'<a href="data:file/txt;base64,{b64}" download="{file_name}">{file_label}</a>'
    return href

def parse_genbank(uploaded_file):
    """Parse a GenBank file and return a list of SeqRecord objects"""
    with tempfile.NamedTemporaryFile(delete=False, suffix=".gb") as tmp:
        tmp.write(uploaded_file.getvalue())
        tmp_path = tmp.name

    try:
        # Use the new function instead of load_genbank
        records = load_genbank_file(tmp_path)
        return records
    finally:
        os.unlink(tmp_path)

def parse_fasta(uploaded_file):
    """Parse a FASTA file and return a list of SeqRecord objects"""
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmp:
        tmp.write(uploaded_file.getvalue())
        tmp_path = tmp.name

    try:
        records = list(SeqIO.parse(tmp_path, "fasta"))
        return records
    finally:
        os.unlink(tmp_path)

def parse_gff(uploaded_file, fasta_file=None):
    """Parse a GFF file and create a database"""
    with tempfile.NamedTemporaryFile(delete=False, suffix=".gff") as tmp:
        tmp.write(uploaded_file.getvalue())
        tmp_path = tmp.name

    try:
        # Create a database
        db_path = f"{tmp_path}.db"
        db = gffutils.create_db(tmp_path, dbfn=db_path, force=True, merge_strategy='create_unique')
        
        # If FASTA file is provided, create SeqRecord objects
        records = []
        if fasta_file is not None:
            fasta_records = parse_fasta(fasta_file)
            
            # Iterate through sequences and add features from GFF
            for fasta_record in fasta_records:
                record = SeqRecord.SeqRecord(
                    seq=fasta_record.seq,
                    id=fasta_record.id,
                    name=fasta_record.name,
                    description=fasta_record.description,
                    features=[]
                )
                
                # Add features from GFF
                try:
                    for feature in db.region(seqid=fasta_record.id):
                        qualifiers = {k: v for k, v in feature.attributes.items()}
                        seq_feature = SeqFeature(
                            location=FeatureLocation(feature.start, feature.end, strand=feature.strand),
                            type=feature.featuretype,
                            qualifiers=qualifiers
                        )
                        record.features.append(seq_feature)
                except ValueError:
                    pass  # No features for this sequence
                
                records.append(record)
        
        return db, records
    finally:
        os.unlink(tmp_path)
        if os.path.exists(f"{tmp_path}.db"):
            os.unlink(f"{tmp_path}.db")

def get_feature_color(feature_type):
    """Return a consistent color for each feature type"""
    color_map = {
        'CDS': '#2C3E50',
        'gene': '#3498DB',
        'mRNA': '#1ABC9C',
        'exon': '#F1C40F',
        'tRNA': '#E74C3C',
        'rRNA': '#9B59B6',
        'repeat_region': '#E67E22',
        'mobile_element': '#34495E',
        'regulatory': '#16A085',
        'misc_feature': '#95A5A6'
    }
    
    return color_map.get(feature_type, '#7F8C8D')

def get_random_color():
    """Generate a random color in hex format"""
    r = lambda: random.randint(0, 200)  # Not too light
    return '#%02X%02X%02X' % (r(), r(), r())

def find_homologous_regions(records, min_identity=70, min_length=100):
    """
    Mock function to find homologous regions between sequences
    In a real implementation, you would use tools like BLAST or similar
    
    Returns a list of tuples (record1_id, start1, end1, record2_id, start2, end2, identity)
    """
    # For demonstration purposes, just return some mock homology data
    homology_data = []
    
    if len(records) < 2:
        return homology_data
    
    for i in range(len(records)):
        for j in range(i+1, len(records)):
            rec1 = records[i]
            rec2 = records[j]
            
            # Create some mock homology regions
            # In a real implementation, you would use sequence alignment tools
            seq1_len = len(rec1.seq)
            seq2_len = len(rec2.seq)
            
            # Generate 1-3 random homologous regions
            for _ in range(random.randint(1, 3)):
                # Random start position in first sequence
                start1 = random.randint(0, max(0, seq1_len - min_length))
                # Random length of homologous region
                region_length = random.randint(min_length, min(seq1_len - start1, 5000))
                end1 = start1 + region_length
                
                # Create a corresponding region in the second sequence
                # For demo, shift it a bit
                shift = random.randint(-200, 200)
                start2 = max(0, min(seq2_len - region_length, start1 + shift))
                end2 = start2 + region_length
                
                # Random identity score (for demonstration)
                identity = random.randint(min_identity, 100)
                
                homology_data.append((rec1.id, start1, end1, rec2.id, start2, end2, identity))
    
    return homology_data

def get_feature_types_multiselect(record, default_types=None):
    """
    Create a multiselect for feature types with validated defaults
    
    Parameters:
    -----------
    record : SeqRecord
        The sequence record to get features from
    default_types : list, optional
        List of feature types to use as defaults
        
    Returns:
    --------
    list
        Selected feature types
    """
    if default_types is None:
        default_types = ["CDS", "gene"]
        
    available_types = list(set(f.type for f in record.features)) if hasattr(record, 'features') else []
    default_selection = [ft for ft in default_types if ft in available_types]
    
    return st.multiselect(
        "Feature types to display:",
        options=available_types,
        default=default_selection
    )

# Main Application
def main():
    st.title("🧬 Advanced Genomic Visualization Hub")
    
    st.markdown("""
    <div class="info-box">
    This application provides advanced tools for visualizing genomic data. Upload your sequences
    in GenBank, FASTA, or GFF format to explore different visualization options.
    </div>
    """, unsafe_allow_html=True)
    
    # Sidebar
    st.sidebar.header("Settings")
    
    # Upload files with new data source options
    with st.sidebar.expander("📤 DATA SOURCE", expanded=True):
        data_source = st.radio("Select data source:", ["Upload File", "URL"])
        
        uploaded_file = None
        fasta_file = None
        
        if data_source == "Upload File":
            file_type = st.radio("Select file type:", ["GenBank", "FASTA", "GFF"])
            uploaded_file = st.file_uploader(f"Upload {file_type} file", type=["gb", "gbk", "fasta", "fa", "fna", "gff", "gff3"])
            
            if file_type == "GFF" and uploaded_file is not None:
                fasta_file = st.file_uploader("Upload corresponding FASTA file", type=["fasta", "fa", "fna"])
        
        else:  # URL
            file_type = st.radio("Select file type:", ["GenBank", "FASTA", "GFF"])
            url = st.text_input(f"Enter {file_type} URL:", "https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/pBAD30.gb")
            
            if url:
                try:
                    import requests
                    response = requests.get(url)
                    if response.status_code == 200:
                        uploaded_file = io.BytesIO(response.content)
                        uploaded_file.name = url.split('/')[-1]
                    else:
                        st.error(f"Failed to fetch URL: {response.status_code}")
                except Exception as e:
                    st.error(f"Error accessing URL: {str(e)}")
            
            if file_type == "GFF" and uploaded_file is not None:
                fasta_url = st.text_input("Enter corresponding FASTA URL:")
                if fasta_url:
                    try:
                        response = requests.get(fasta_url)
                        if response.status_code == 200:
                            fasta_file = io.BytesIO(response.content)
                            fasta_file.name = fasta_url.split('/')[-1]
                        else:
                            st.error(f"Failed to fetch FASTA URL: {response.status_code}")
                    except Exception as e:
                        st.error(f"Error accessing FASTA URL: {str(e)}")
            
            # Add sample URLs for convenience
            st.markdown("### Sample URLs:")
            if file_type == "GenBank":
                st.markdown("""
                * [pBAD30 Plasmid](https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/pBAD30.gb)
                * [E. coli K12 segment](https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_000913.gb)
                * [Mycoplasma](https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_000908.gb)
                """)
            elif file_type == "FASTA":
                st.markdown("""
                * [pBAD30 Plasmid](https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/pBAD30.fasta)
                * [Mycoplasma](https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_000908.fasta)
                """)
            elif file_type == "GFF":
                st.markdown("""
                * GFF Example URL (need to provide both GFF and FASTA URLs)
                """)

    # Parse uploaded file
    records = []
    gff_db = None
    if uploaded_file is not None:
        if file_type == "GenBank":
            records = parse_genbank(uploaded_file)
        elif file_type == "FASTA":
            records = parse_fasta(uploaded_file)
        elif file_type == "GFF":
            if fasta_file is not None:
                gff_db, records = parse_gff(uploaded_file, fasta_file)
            else:
                gff_db, _ = parse_gff(uploaded_file)
                st.sidebar.warning("To visualize features from GFF, please upload a corresponding FASTA file.")
    
    # Display information about loaded data
    if records:
        with st.sidebar.expander("📊 LOADED DATA", expanded=True):
            st.write(f"Number of sequences: {len(records)}")
            if len(records) > 0:
                seq_info = pd.DataFrame({
                    "ID": [rec.id for rec in records],
                    "Length": [len(rec.seq) for rec in records],
                    "Features": [len(rec.features) if hasattr(rec, 'features') else 0 for rec in records]
                })
                st.dataframe(seq_info)
    
    # Visualization options
    tabs = st.tabs([
        "📌 Genome Map", 
        "🔄 Multi-Genome Comparison", 
        "🧩 Feature Analysis", 
        "📊 GC Content & Coverage", 
        "🌐 Synteny Network"
    ])
    
    # 1. Genome Map Tab
    with tabs[0]:
        st.header("Genomic Feature Map")
        
        if not records:
            st.info("Please upload a sequence file to visualize genomic features.")
        else:
            col1, col2 = st.columns([3, 1])
            
            with col2:
                st.subheader("Visualization Settings")
                
                # Select sequence to visualize
                sequence_options = [rec.id for rec in records]
                selected_sequence = st.selectbox("Select sequence:", sequence_options)
                
                # Get the selected record
                selected_record = next((rec for rec in records if rec.id == selected_sequence), None)
                
                if selected_record:
                    # Get all available feature types from the record
                    available_feature_types = list(set(f.type for f in selected_record.features)) if hasattr(selected_record, 'features') else []

                    # Set default values only if they exist in available options
                    default_types = ["CDS", "gene"]
                    default_selection = [ft for ft in default_types if ft in available_feature_types]

                    # Create the multiselect with validated defaults
                    feature_types = st.multiselect(
                        "Feature types to display:",
                        options=available_feature_types,
                        default=default_selection
                    )
                    
                    plot_height = st.slider("Plot height:", 300, 1000, 500)
                    show_labels = st.checkbox("Show feature labels", value=True)
                    label_angle = st.slider("Label angle:", 0, 90, 45)
                    plot_style = st.selectbox("Plot style:", ["Linear", "Circular"])
                
            with col1:
                if selected_record and hasattr(selected_record, 'features') and len(selected_record.features) > 0:
                    st.subheader(f"Genome Map for {selected_record.id}")
                    
                    # Filter features by type
                    filtered_features = [
                        f for f in selected_record.features 
                        if f.type in feature_types
                    ] if feature_types else selected_record.features
                    
                    if not filtered_features:
                        st.warning("No features of the selected types found in this sequence.")
                    else:
                        # Create visualization
                        if plot_style == "Linear":
                            # Use DNA Features Viewer for linear visualization
                            class CustomTranslator(BiopythonTranslator):
                                def compute_feature_color(self, feature):
                                    return get_feature_color(feature.type)
                                
                                def compute_feature_label(self, feature):
                                    if not show_labels:
                                        return None
                                    for key in ['gene', 'locus_tag', 'product', 'label', 'name']:
                                        if key in feature.qualifiers:
                                            return feature.qualifiers[key][0]
                                    return f"{feature.type}_{feature.location.start}..{feature.location.end}"
                            
                            translator = CustomTranslator()
                            graphic_record = translator.translate_record(selected_record)
                            
                            # Linear plot
                            fig, ax = plt.subplots(1, figsize=(12, plot_height/100))
                            graphic_record.plot(ax=ax, with_ruler=True)
                            
                            if show_labels:
                                for feature in graphic_record.features:
                                    if feature.label is not None:
                                        ax.text(
                                            feature.start + (feature.end - feature.start) / 2,
                                            0,
                                            feature.label,
                                            horizontalalignment='center',
                                            verticalalignment='center' if feature.strand == 0 else ('top' if feature.strand > 0 else 'bottom'),
                                            fontsize=8,
                                            rotation=label_angle,
                                            transform=ax.get_xaxis_transform()
                                        )
                            
                            # Add legend
                            legend_elements = [
                                mpatches.Patch(color=get_feature_color(f_type), label=f_type)
                                for f_type in set(f.type for f in filtered_features)
                            ]
                            ax.legend(handles=legend_elements, loc='upper center', 
                                     bbox_to_anchor=(0.5, -0.15), ncol=min(5, len(legend_elements)))
                            
                            plt.tight_layout()
                            st.pyplot(fig)
                            
                            # Save plot option
                            buf = io.BytesIO()
                            fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
                            buf.seek(0)
                            st.download_button(
                                label="Download PNG",
                                data=buf,
                                file_name=f"{selected_record.id}_genome_map.png",
                                mime="image/png"
                            )
                            
                            # Additionally add SVG option
                            buf_svg = io.BytesIO()
                            fig.savefig(buf_svg, format='svg', bbox_inches='tight')
                            buf_svg.seek(0)
                            st.download_button(
                                label="Download SVG",
                                data=buf_svg,
                                file_name=f"{selected_record.id}_genome_map.svg",
                                mime="image/svg+xml",
                                key="svg_download"
                            )
                            
                        else:  # Circular plot
                            # Use matplotlib for circular visualization
                            fig, ax = plt.subplots(
                                1, 1, figsize=(10, 10), 
                                subplot_kw={'projection': 'polar'}
                            )
                            
                            # Total sequence length
                            seq_length = len(selected_record.seq)
                            
                            # Calculate positions in radians
                            feature_types_set = set(f.type for f in filtered_features)
                            type_to_track = {t: i for i, t in enumerate(feature_types_set)}
                            num_tracks = len(type_to_track)
                            
                            # Tracks for different feature types (from inner to outer)
                            track_width = 1.0 / (num_tracks + 1)  # +1 for some spacing
                            
                            # Plot each feature
                            for feature in filtered_features:
                                start = int(feature.location.start)
                                end = int(feature.location.end)
                                track = type_to_track[feature.type]
                                
                                # Convert to radians
                                start_rad = 2 * np.pi * start / seq_length
                                end_rad = 2 * np.pi * end / seq_length
                                
                                # Inner and outer radius for this track
                                inner_r = 0.3 + track * track_width
                                outer_r = inner_r + track_width * 0.8
                                
                                # Draw feature arc
                                theta = np.linspace(start_rad, end_rad, 100)
                                radii = np.ones_like(theta) * outer_r
                                ax.fill_between(theta, inner_r, radii, color=get_feature_color(feature.type), alpha=0.8)
                                
                                # Add label if enabled
                                if show_labels and (end - start) / seq_length > 0.01:  # Only label features that are large enough
                                    label = None
                                    for key in ['gene', 'locus_tag', 'product', 'label', 'name']:
                                        if key in feature.qualifiers:
                                            label = feature.qualifiers[key][0]
                                            break
                                    
                                    if label:
                                        mid_rad = (start_rad + end_rad) / 2
                                        mid_r = (inner_r + outer_r) / 2
                                        
                                        ax.text(
                                            mid_rad, mid_r, label,
                                            horizontalalignment='center',
                                            verticalalignment='center',
                                            rotation=np.rad2deg(mid_rad) - 90,
                                            fontsize=8,
                                            rotation_mode='anchor'
                                        )
                            
                            # Set up the axes
                            ax.set_theta_zero_location('N')
                            ax.set_theta_direction(-1)  # Clockwise
                            ax.set_rticks([])  # No radial ticks
                            
                            # Add a grid for the major divisions
                            ax.set_xticks(np.linspace(0, 2*np.pi, 13)[:-1])  # 12 divisions (0-indexed)
                            ax.grid(True)
                            
                            # Add sequence position labels at the edge
                            positions = np.linspace(0, seq_length, 13)[:-1]  # 12 positions (0-indexed)
                            for i, pos in enumerate(positions):
                                angle = 2 * np.pi * pos / seq_length
                                ax.text(
                                    angle, 1.05, f"{int(pos):,}",
                                    horizontalalignment='center',
                                    verticalalignment='center',
                                    fontsize=8
                                )
                            
                            # Add legend
                            legend_elements = [
                                mpatches.Patch(color=get_feature_color(f_type), label=f_type)
                                for f_type in feature_types_set
                            ]
                            ax.legend(handles=legend_elements, loc='center')
                            
                            plt.title(f"Circular Genome Map - {selected_record.id}")
                            st.pyplot(fig)
                            
                            # Save plot option
                            buf = io.BytesIO()
                            fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
                            buf.seek(0)
                            st.download_button(
                                label="Download PNG",
                                data=buf,
                                file_name=f"{selected_record.id}_circular_map.png",
                                mime="image/png"
                            )
                else:
                    st.info("This sequence has no annotated features.")
    
    # 2. Multi-Genome Comparison Tab
    with tabs[1]:
        st.header("Multi-Genome Comparison")
        
        if len(records) < 2:
            st.info("Please upload a file with multiple sequences to compare genomes.")
        else:
            col1, col2 = st.columns([3, 1])
            
            with col2:
                st.subheader("Comparison Settings")
                
                # Select sequences to compare
                sequence_options = [rec.id for rec in records]
                selected_sequences = st.multiselect(
                    "Select sequences to compare:",
                    options=sequence_options,
                    default=sequence_options[:min(5, len(sequence_options))]
                )
                
                # Homology settings
                st.write("Homology Settings:")
                min_identity = st.slider("Minimum identity %:", 50, 100, 70)
                min_length = st.slider("Minimum alignment length:", 50, 1000, 100)
                
                # Visualization settings
                link_color = st.selectbox("Link color scheme:", ["Viridis", "Plasma", "Blues", "Reds"])
                show_feature_labels = st.checkbox("Show feature labels", value=False)
                
                # Features to highlight
                feature_types_to_show = st.multiselect(
                    "Feature types to highlight:",
                    options=list(set(f.type for rec in records for f in rec.features if hasattr(rec, 'features'))),
                    default=["CDS"] if any(hasattr(rec, 'features') and any(f.type == "CDS" for f in rec.features) for rec in records) else []
                )
            
            with col1:
                if selected_sequences and len(selected_sequences) >= 2:
                    st.subheader("Genome Comparison Visualization")
                    
                    # Get selected records
                    selected_records = [rec for rec in records if rec.id in selected_sequences]
                    
                    # Find homologous regions between selected sequences
                    homology_data = find_homologous_regions(
                        selected_records, 
                        min_identity=min_identity,
                        min_length=min_length
                    )
                    
                    # Use pygenomeviz to create the comparison visualization
                    fig = pgv.Visualization(name="Genome Comparison")
                    
                    # Add tracks for each selected genome
                    tracks = {}
                    cmap = plt.get_cmap('tab10')
                    
                    for i, record in enumerate(selected_records):
                        # Create track
                        tracks[record.id] = fig.add_feature_track(
                            name=record.id,
                            data=[
                                (record.id, 0, len(record.seq))
                            ],
                            track_color=mcolors.rgb2hex(cmap(i % 10))
                        )
                        
                        # Add features if available
                        if hasattr(record, 'features'):
                            for feature in record.features:
                                if feature.type in feature_types_to_show:
                                    start = int(feature.location.start)
                                    end = int(feature.location.end)
                                    strand = feature.location.strand
                                    
                                    # Get label if available
                                    label = None
                                    if show_feature_labels:
                                        for key in ['gene', 'locus_tag', 'product', 'name']:
                                            if key in feature.qualifiers:
                                                label = feature.qualifiers[key][0]
                                                break
                                    
                                    tracks[record.id].add_feature(
                                        start=start,
                                        end=end,
                                        strand=1 if strand == 1 else -1 if strand == -1 else 0,
                                        label=label,
                                        color=get_feature_color(feature.type)
                                    )
                    
                    # Add links between homologous regions
                    color_map = {
                        'Viridis': cm.viridis,
                        'Plasma': cm.plasma,
                        'Blues': cm.Blues,
                        'Reds': cm.Reds
                    }
                    
                    colormap = color_map.get(link_color, cm.viridis)
                    
                    for homology in homology_data:
                        rec1_id, start1, end1, rec2_id, start2, end2, identity = homology
                        
                        # Scale color by identity
                        color = mcolors.rgb2hex(colormap(identity / 100))
                        
                        # Add link
                        fig.add_link(
                            track_name_A=rec1_id,
                            start_A=start1,
                            end_A=end1,
                            track_name_B=rec2_id,
                            start_B=start2,
                            end_B=end2,
                            color=color,
                            alpha=0.7
                        )
                    
                    # Create the plot
                    fig_width = min(15, 8 + len(selected_sequences) * 0.5)
                    fig_height = min(10, 3 + len(selected_sequences) * 0.5)
                    
                    fig.plot(
                        feature_track_ratio=0.25,
                        link_track_ratio=1.0,
                        tick_style="axis",
                        tick_labelsize=8,
                        feature_labelsize=6,
                        track_labelsize=10,
                        fig_width=fig_width,
                        fig_height=fig_height
                    )
                    
                    # Display the plot
                    st.pyplot(fig.get_figure())
                    
                    # Add color legend for identity
                    fig_legend, ax_legend = plt.subplots(figsize=(6, 0.8))
                    cbar = plt.colorbar(
                        plt.cm.ScalarMappable(cmap=colormap), 
                        cax=ax_legend, 
                        orientation='horizontal'
                    )
                    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
                    cbar.set_ticklabels([f"{min_identity}%", "", "", "", "100%"])
                    cbar.set_label('Sequence Identity')
                    st.pyplot(fig_legend)
                    
                    # Feature legend if needed
                    if feature_types_to_show:
                        fig_feat_legend, ax_feat_legend = plt.subplots(figsize=(6, 0.5))
                        ax_feat_legend.axis('off')
                        legend_elements = [
                            mpatches.Patch(color=get_feature_color(f_type), label=f_type)
                            for f_type in feature_types_to_show
                        ]
                        ax_feat_legend.legend(
                            handles=legend_elements, 
                            loc='center', 
                            ncol=min(5, len(feature_types_to_show)),
                            frameon=False
                        )
                        st.pyplot(fig_feat_legend)
                    
                    # Save options
                    buffer = io.BytesIO()
                    fig.get_figure().savefig(buffer, format='png', dpi=300, bbox_inches='tight')
                    buffer.seek(0)
                    
                    st.download_button(
                        label="Download Comparison Image (PNG)",
                        data=buffer,
                        file_name="genome_comparison.png",
                        mime="image/png"
                    )
                    
                    # Add SVG option
                    buffer_svg = io.BytesIO()
                    fig.get_figure().savefig(buffer_svg, format='svg', bbox_inches='tight')
                    buffer_svg.seek(0)
                    
                    st.download_button(
                        label="Download Comparison Image (SVG)",
                        data=buffer_svg,
                        file_name="genome_comparison.svg",
                        mime="image/svg+xml",
                        key="comparison_svg"
                    )
    
    # 3. Feature Analysis Tab
    with tabs[2]:
        st.header("Feature Analysis")
        
        if not records or not any(hasattr(rec, 'features') and len(rec.features) > 0 for rec in records):
            st.info("Please upload a file with annotated features to analyze.")
        else:
            # Filter records that have features
            records_with_features = [rec for rec in records if hasattr(rec, 'features') and len(rec.features) > 0]
            
            col1, col2 = st.columns([3, 1])
            
            with col2:
                st.subheader("Analysis Settings")
                
                # Select sequence to analyze
                sequence_options = [rec.id for rec in records_with_features]
                selected_sequence = st.selectbox(
                    "Select sequence to analyze:",
                    options=sequence_options,
                    key="feature_analysis_seq"
                )
                
                # Get the selected record
                selected_record = next((rec for rec in records if rec.id == selected_sequence), None)
                
                if selected_record:
                    # Get all feature types
                    feature_types = list(set(f.type for f in selected_record.features))
                    
                    # Select feature types to analyze
                    selected_feature_types = st.multiselect(
                        "Feature types to analyze:",
                        options=feature_types,
                        default=["CDS"] if "CDS" in feature_types else [feature_types[0]] if feature_types else []
                    )
                    
                    # Analysis options
                    analysis_type = st.radio(
                        "Analysis type:",
                        ["Feature Distribution", "Feature Length Statistics", "Feature Density", "Strand Bias"]
                    )
            
            with col1:
                if selected_record and selected_feature_types:
                    # Filter features by selected types
                    filtered_features = [
                        f for f in selected_record.features 
                        if f.type in selected_feature_types
                    ]
                    
                    if not filtered_features:
                        st.warning("No features of the selected types found in this sequence.")
                    else:
                        # Analyze features based on selected analysis type
                        if analysis_type == "Feature Distribution":
                            st.subheader(f"Feature Distribution on {selected_record.id}")
                            
                            # Create a dataframe with feature positions
                            feature_data = []
                            for feature in filtered_features:
                                feature_type = feature.type
                                start = int(feature.location.start)
                                end = int(feature.location.end)
                                strand = feature.location.strand
                                strand_str = "+" if strand == 1 else "-" if strand == -1 else "none"
                                
                                # Get feature name if available
                                name = None
                                for key in ['gene', 'locus_tag', 'product', 'name']:
                                    if key in feature.qualifiers:
                                        name = feature.qualifiers[key][0]
                                        break
                                if not name:
                                    name = f"{feature_type}_{start}..{end}"
                                
                                feature_data.append({
                                    'Type': feature_type,
                                    'Start': start,
                                    'End': end,
                                    'Length': end - start,
                                    'Strand': strand_str,
                                    'Name': name
                                })
                            
                            df = pd.DataFrame(feature_data)
                            
                            # Create interactive plotly figure
                            fig = px.scatter(
                                df,
                                x='Start',
                                y='Type',
                                color='Strand',
                                size='Length',
                                hover_name='Name',
                                hover_data={
                                    'Start': True,
                                    'End': True,
                                    'Length': True,
                                    'Type': False,  # Already in y-axis
                                    'Strand': True
                                },
                                color_discrete_map={
                                    '+': '#1f77b4',  # Forward strand - blue
                                    '-': '#d62728',  # Reverse strand - red
                                    'none': '#7f7f7f'  # No strand - gray
                                },
                                title=f"Feature Distribution on {selected_record.id}",
                                labels={
                                    'Start': 'Position (bp)',
                                    'Type': 'Feature Type',
                                    'Length': 'Feature Length (bp)'
                                },
                                height=500
                            )
                            
                            # Add sequence length line
                            fig.add_shape(
                                type="line",
                                x0=0,
                                y0=-0.5,
                                x1=len(selected_record.seq),
                                y1=-0.5,
                                line=dict(color="black", width=2)
                            )
                            
                            # Add ticks every 10% of the sequence
                            seq_len = len(selected_record.seq)
                            for i in range(0, 11):
                                pos = i * seq_len / 10
                                fig.add_shape(
                                    type="line",
                                    x0=pos,
                                    y0=-0.5,
                                    x1=pos,
                                    y1=-0.7,
                                    line=dict(color="black", width=1)
                                )
                                
                                # Add label for every 25%
                                if i % 2.5 == 0:
                                    fig.add_annotation(
                                        x=pos,
                                        y=-0.9,
                                        text=f"{int(pos):,} bp",
                                        showarrow=False,
                                        font=dict(size=10)
                                    )
                            
                            # Update layout
                            fig.update_layout(
                                xaxis_title="Position (bp)",
                                yaxis_title="Feature Type",
                                legend_title="Strand",
                                height=400 + len(selected_feature_types) * 50
                            )
                            
                            st.plotly_chart(fig, use_container_width=True)
                            
                        elif analysis_type == "Feature Length Statistics":
                            st.subheader(f"Feature Length Statistics for {selected_record.id}")
                            
                            # Create a dataframe with feature lengths
                            feature_data = []
                            for feature in filtered_features:
                                feature_type = feature.type
                                start = int(feature.location.start)
                                end = int(feature.location.end)
                                length = end - start
                                strand = feature.location.strand
                                strand_str = "+" if strand == 1 else "-" if strand == -1 else "none"
                                
                                # Get feature name if available
                                name = None
                                for key in ['gene', 'locus_tag', 'product', 'name']:
                                    if key in feature.qualifiers:
                                        name = feature.qualifiers[key][0]
                                        break
                                if not name:
                                    name = f"{feature_type}_{start}..{end}"
                                
                                feature_data.append({
                                    'Type': feature_type,
                                    'Length': length,
                                    'Strand': strand_str,
                                    'Name': name
                                })
                            
                            df = pd.DataFrame(feature_data)
                            
                            # Create a box plot of feature lengths by type
                            fig1 = px.box(
                                df,
                                x='Type',
                                y='Length',
                                color='Type',
                                points="all",
                                hover_name='Name',
                                labels={
                                    'Type': 'Feature Type',
                                    'Length': 'Feature Length (bp)'
                                },
                                title=f"Feature Length Distribution by Type"
                            )
                            
                            # Update layout
                            fig1.update_layout(
                                xaxis_title="Feature Type",
                                yaxis_title="Length (bp)",
                                showlegend=False
                            )
                            
                            st.plotly_chart(fig1, use_container_width=True)
                            
                            # Create a histogram of feature lengths
                            fig2 = px.histogram(
                                df,
                                x='Length',
                                color='Type',
                                marginal="rug",
                                hover_name='Name',
                                barmode='overlay',
                                opacity=0.7,
                                labels={
                                    'Length': 'Feature Length (bp)',
                                    'count': 'Count'
                                },
                                title=f"Feature Length Histogram"
                            )
                            
                            # Update layout
                            fig2.update_layout(
                                xaxis_title="Length (bp)",
                                yaxis_title="Count",
                                legend_title="Feature Type"
                            )
                            
                            st.plotly_chart(fig2, use_container_width=True)
                            
                            # Display summary statistics
                            st.subheader("Summary Statistics")
                            
                            stats = df.groupby('Type')['Length'].agg([
                                ('Count', 'count'),
                                ('Mean', 'mean'),
                                ('Median', 'median'),
                                ('Min', 'min'),
                                ('Max', 'max'),
                                ('Std', 'std')
                            ]).reset_index()
                            
                            # Format numeric columns
                            for col in ['Mean', 'Median', 'Std']:
                                stats[col] = stats[col].round(2)
                            
                            st.dataframe(stats, use_container_width=True)
                            
                        elif analysis_type == "Feature Density":
                            st.subheader(f"Feature Density Analysis for {selected_record.id}")
                            
                            # Sequence length
                            seq_len = len(selected_record.seq)
                            
                            # Window size options
                            window_size = st.slider(
                                "Window size (bp):",
                                min_value=1000,
                                max_value=seq_len // 10,
                                value=min(5000, seq_len // 10),
                                step=1000
                            )
                            
                            # Step size options
                            step_size = st.slider(
                                "Step size (bp):",
                                min_value=window_size // 10,
                                max_value=window_size,
                                value=window_size // 2,
                                step=window_size // 10
                            )
                            
                            # Calculate feature density across the genome
                            positions = list(range(0, seq_len - window_size + 1, step_size))
                            density_data = {feature_type: [] for feature_type in selected_feature_types}
                            
                            for pos in positions:
                                # Count features in this window for each type
                                for feature_type in selected_feature_types:
                                    # Filter features of this type
                                    type_features = [f for f in filtered_features if f.type == feature_type]
                                    
                                    # Count features that overlap with this window
                                    count = sum(
                                        1 for f in type_features
                                        if (int(f.location.start) <= pos + window_size and 
                                            int(f.location.end) >= pos)
                                    )
                                    
                                    density_data[feature_type].append(count)
                            
                            # Create a dataframe for plotting
                            df_density = pd.DataFrame({
                                'Position': [p + window_size/2 for p in positions]  # Center of the window
                            })
                            
                            for feature_type in selected_feature_types:
                                df_density[feature_type] = density_data[feature_type]
                            
                            # Create an interactive plotly figure
                            fig = go.Figure()
                            
                            # Add a trace for each feature type
                            for i, feature_type in enumerate(selected_feature_types):
                                color = get_feature_color(feature_type)
                                
                                fig.add_trace(go.Scatter(
                                    x=df_density['Position'],
                                    y=df_density[feature_type],
                                    mode='lines',
                                    name=feature_type,
                                    line=dict(color=color, width=2),
                                    fill='tozeroy',
                                    fillcolor=f'rgba{tuple(list(mcolors.to_rgb(color)) + [0.2])}'
                                ))
                            
                            # Update layout
                            fig.update_layout(
                                title=f"Feature Density (Window Size: {window_size:,} bp, Step Size: {step_size:,} bp)",
                                xaxis_title="Position (bp)",
                                yaxis_title="Feature Count",
                                hovermode="x unified",
                                legend_title="Feature Type"
                            )
                            
                            st.plotly_chart(fig, use_container_width=True)
                            
                            # Add explanation
                            st.markdown("""
                            <div class="info-box">
                            <strong>About Feature Density Analysis:</strong><br>
                            This analysis shows the distribution of features along the genome using a sliding window approach.
                            Each point represents the number of features that overlap with a window of the specified size.
                            The x-axis represents the position along the genome, and the y-axis shows the feature count.
                            </div>
                            """, unsafe_allow_html=True)
                            
                        elif analysis_type == "Strand Bias":
                            st.subheader(f"Strand Bias Analysis for {selected_record.id}")
                            
                            # Create a dataframe with feature strand information
                            feature_data = []
                            for feature in filtered_features:
                                feature_type = feature.type
                                start = int(feature.location.start)
                                end = int(feature.location.end)
                                length = end - start
                                strand = feature.location.strand
                                strand_str = "Forward" if strand == 1 else "Reverse" if strand == -1 else "None"
                                
                                feature_data.append({
                                    'Type': feature_type,
                                    'Start': start,
                                    'End': end,
                                    'Length': length,
                                    'Strand': strand_str
                                })
                            
                            df = pd.DataFrame(feature_data)
                            
                            # Overall strand bias
                            st.subheader("Overall Strand Distribution")
                            
                            strand_counts = df['Strand'].value_counts().reset_index()
                            strand_counts.columns = ['Strand', 'Count']
                            
                            # Calculate percentages
                            total = strand_counts['Count'].sum()
                            strand_counts['Percentage'] = (strand_counts['Count'] / total * 100).round(2)
                            
                            # Create pie chart
                            fig1 = px.pie(
                                strand_counts,
                                values='Count',
                                names='Strand',
                                title="Overall Strand Distribution",
                                color='Strand',
                                color_discrete_map={
                                    'Forward': '#1f77b4',  # Blue
                                    'Reverse': '#d62728',  # Red
                                    'None': '#7f7f7f'      # Gray
                                },
                                hover_data=['Percentage']
                            )
                            
                            st.plotly_chart(fig1, use_container_width=True)
                            
                            # Strand bias by feature type
                            st.subheader("Strand Distribution by Feature Type")
                            
                            # Create grouped bar chart
                            strand_type_counts = df.groupby(['Type', 'Strand']).size().reset_index(name='Count')
                            
                            fig2 = px.bar(
                                strand_type_counts,
                                x='Type',
                                y='Count',
                                color='Strand',
                                barmode='group',
                                title="Strand Distribution by Feature Type",
                                color_discrete_map={
                                    'Forward': '#1f77b4',  # Blue
                                    'Reverse': '#d62728',  # Red
                                    'None': '#7f7f7f'      # Gray
                                },
                            )
                            
                            # Update layout
                            fig2.update_layout(
                                xaxis_title="Feature Type",
                                yaxis_title="Count",
                                legend_title="Strand"
                            )
                            
                            st.plotly_chart(fig2, use_container_width=True)
                            
                            # Calculate replichore skew if we have both forward and reverse features
                            if 'Forward' in df['Strand'].values and 'Reverse' in df['Strand'].values:
                                st.subheader("Replichore Skew Analysis")
                                
                                st.markdown("""
                                <div class="info-box">
                                <strong>About Replichore Skew:</strong><br>
                                Replichore skew analysis examines the distribution of genes on the leading and lagging strands 
                                of DNA replication. In many bacteria, essential and highly expressed genes tend to be coded on 
                                the leading strand to avoid head-on collisions between DNA and RNA polymerases.
                                </div>
                                """, unsafe_allow_html=True)
                                
                                # Cumulative strand bias along the genome
                                # Window size options
                                seq_len = len(selected_record.seq)
                                window_size = st.slider(
                                    "Window size (bp):",
                                    min_value=1000,
                                    max_value=seq_len // 5,
                                    value=min(10000, seq_len // 5),
                                    step=1000,
                                    key="skew_window_size"
                                )
                                
                                # Step size options
                                step_size = st.slider(
                                    "Step size (bp):",
                                    min_value=window_size // 10,
                                    max_value=window_size,
                                    value=window_size // 2,
                                    step=window_size // 10,
                                    key="skew_step_size"
                                )
                                
                                # Calculate strand bias across the genome
                                positions = list(range(0, seq_len - window_size + 1, step_size))
                                forward_counts = []
                                reverse_counts = []
                                skew_values = []
                                
                                for pos in positions:
                                    # Count features in this window for each strand
                                    forward_features = [
                                        f for f in filtered_features
                                        if f.location.strand == 1 and
                                        int(f.location.start) <= pos + window_size and 
                                        int(f.location.end) >= pos
                                    ]
                                    
                                    reverse_features = [
                                        f for f in filtered_features
                                        if f.location.strand == -1 and
                                        int(f.location.start) <= pos + window_size and 
                                        int(f.location.end) >= pos
                                    ]
                                    
                                    forward_count = len(forward_features)
                                    reverse_count = len(reverse_features)
                                    
                                    forward_counts.append(forward_count)
                                    reverse_counts.append(reverse_count)
                                    
                                    # Calculate GC skew
                                    if forward_count + reverse_count > 0:
                                        skew = (forward_count - reverse_count) / (forward_count + reverse_count)
                                    else:
                                        skew = 0
                                    
                                    skew_values.append(skew)
                                
                                # Create a dataframe for plotting
                                df_skew = pd.DataFrame({
                                    'Position': [p + window_size/2 for p in positions],  # Center of the window
                                    'Forward': forward_counts,
                                    'Reverse': reverse_counts,
                                    'Skew': skew_values
                                })
                                
                                # Create an interactive plotly figure for the skew
                                fig3 = go.Figure()
                                
                                fig3.add_trace(go.Scatter(
                                    x=df_skew['Position'],
                                    y=df_skew['Skew'],
                                    mode='lines',
                                    name='Strand Skew',
                                    line=dict(color='#2ca02c', width=2),  # Green
                                    fill='tozeroy',
                                    fillcolor='rgba(44, 160, 44, 0.2)'
                                ))
                                
                                # Add a horizontal line at y=0
                                fig3.add_shape(
                                    type="line",
                                    x0=0,
                                    y0=0,
                                    x1=seq_len,
                                    y1=0,
                                    line=dict(color="black", width=1, dash="dash")
                                )
                                
                                # Calculate cumulative skew
                                cumulative_skew = np.cumsum(skew_values)
                                
                                # Add cumulative skew trace
                                fig3.add_trace(go.Scatter(
                                    x=df_skew['Position'],
                                    y=cumulative_skew,
                                    mode='lines',
                                    name='Cumulative Skew',
                                    line=dict(color='#ff7f0e', width=2),
                                    visible='legendonly'  # Hidden by default
                                ))
                                
                                # Update layout
                                fig3.update_layout(
                                    title=f"Strand Skew Analysis (Window Size: {window_size:,} bp, Step Size: {step_size:,} bp)",
                                    xaxis_title="Position (bp)",
                                    yaxis_title="Strand Skew (F-R)/(F+R)",
                                    hovermode="x unified"
                                )
                                
                                # Add annotations for potential origin and terminus
                                if len(cumulative_skew) > 0:
                                    origin_idx = np.argmin(cumulative_skew)
                                    terminus_idx = np.argmax(cumulative_skew)
                                    
                                    origin_pos = df_skew.iloc[origin_idx]['Position']
                                    terminus_pos = df_skew.iloc[terminus_idx]['Position']
                                    
                                    # Add vertical lines and annotations
                                    fig3.add_vline(
                                        x=origin_pos,
                                        line_width=1,
                                        line_dash="dash",
                                        line_color="blue",
                                        annotation_text="Potential Origin",
                                        annotation_position="top"
                                    )
                                    
                                    fig3.add_vline(
                                        x=terminus_pos,
                                        line_width=1,
                                        line_dash="dash",
                                        line_color="red",
                                        annotation_text="Potential Terminus",
                                        annotation_position="top"
                                    )
                                
                                st.plotly_chart(fig3, use_container_width=True)
                                
                                # Add explanation
                                st.markdown("""
                                <div class="info-box">
                                <strong>Interpretation:</strong><br>
                                - <b>Strand Skew</b>: (Forward - Reverse)/(Forward + Reverse) for each window.
                                - <b>Cumulative Skew</b>: The running sum of skew values, which often helps identify the origin and terminus of replication.
                                - <b>Potential Origin</b>: Often located at the minimum of the cumulative skew curve.
                                - <b>Potential Terminus</b>: Often located at the maximum of the cumulative skew curve.
                                
                                Note that this is a simplistic approach and actual origin/terminus identification
                                would require additional analysis, such as GC skew or DnaA box identification.
                                </div>
                                """, unsafe_allow_html=True)
                                
                                # Now add counts by strand
                                fig4 = go.Figure()
                                
                                fig4.add_trace(go.Scatter(
                                    x=df_skew['Position'],
                                    y=df_skew['Forward'],
                                    mode='lines',
                                    name='Forward',
                                    line=dict(color='#1f77b4', width=2),  # Blue
                                    fill='tozeroy',
                                    fillcolor='rgba(31, 119, 180, 0.2)'
                                ))
                                
                                fig4.add_trace(go.Scatter(
                                    x=df_skew['Position'],
                                    y=df_skew['Reverse'],
                                    mode='lines',
                                    name='Reverse',
                                    line=dict(color='#d62728', width=2),  # Red
                                    fill='tozeroy',
                                    fillcolor='rgba(214, 39, 40, 0.2)'
                                ))
                                
                                # Update layout
                                fig4.update_layout(
                                    title=f"Feature Counts by Strand (Window Size: {window_size:,} bp)",
                                    xaxis_title="Position (bp)",
                                    yaxis_title="Feature Count",
                                    hovermode="x unified"
                                )
                                
                                st.plotly_chart(fig4, use_container_width=True)
    
    # 4. GC Content & Coverage Tab
    with tabs[3]:
        st.header("GC Content & Coverage Analysis")
        
        if not records:
            st.info("Please upload a sequence file to analyze GC content.")
        else:
            col1, col2 = st.columns([3, 1])
            
            with col2:
                st.subheader("Analysis Settings")
                
                # Select sequence to analyze
                sequence_options = [rec.id for rec in records]
                selected_sequence = st.selectbox(
                    "Select sequence:",
                    options=sequence_options,
                    key="gc_analysis_seq"
                )
                
                # Get the selected record
                selected_record = next((rec for rec in records if rec.id == selected_sequence), None)
                
                if selected_record:
                    # GC Content settings
                    st.subheader("GC Content Settings")
                    
                    window_size = st.slider(
                        "Window size (bp):",
                        min_value=100,
                        max_value=min(10000, len(selected_record.seq) // 10),
                        value=1000,
                        step=100
                    )
                    
                    step_size = st.slider(
                        "Step size (bp):",
                        min_value=window_size // 10,
                        max_value=window_size,
                        value=window_size // 2,
                        step=window_size // 10
                    )
                    
                    # Coverage analysis
                    st.subheader("Coverage Analysis")
                    st.info("For coverage analysis, please upload a read alignment file (coming in future update).")
                    
                    # Additional analysis options
                    show_skew = st.checkbox("Show GC Skew", value=True)
                    show_cumulative = st.checkbox("Show Cumulative GC Skew", value=True)
                    highlight_regions = st.checkbox("Highlight Feature Regions", value=False)
            
            with col1:
                if selected_record:
                    st.subheader(f"GC Content Analysis for {selected_record.id}")
                    
                    # Calculate GC content across the genome
                    seq_str = str(selected_record.seq).upper()
                    seq_len = len(seq_str)
                    positions = list(range(0, seq_len - window_size + 1, step_size))
                    gc_content = []
                    gc_skew = []
                    
                    for pos in positions:
                        window = seq_str[pos:pos + window_size]
                        g_count = window.count('G')
                        c_count = window.count('C')
                        a_count = window.count('A')
                        t_count = window.count('T')
                        
                        # Calculate GC content
                        gc = (g_count + c_count) / len(window) * 100 if len(window) > 0 else 0
                        gc_content.append(gc)
                        
                        # Calculate GC skew
                        skew = (g_count - c_count) / (g_count + c_count) if (g_count + c_count) > 0 else 0
                        gc_skew.append(skew)
                    
                    # Create a dataframe for plotting
                    df_gc = pd.DataFrame({
                        'Position': [p + window_size/2 for p in positions],  # Center of the window
                        'GC': gc_content,
                        'GC_Skew': gc_skew
                    })
                    
                    # Create combined plot with GC content and skew
                    fig = make_subplots(
                        rows=2 if show_skew else 1,
                        cols=1,
                        shared_xaxes=True,
                        vertical_spacing=0.1,
                        subplot_titles=(
                            "GC Content",
                            "GC Skew" if show_skew else None
                        )
                    )
                    
                    # Add GC content trace
                    fig.add_trace(
                        go.Scatter(
                            x=df_gc['Position'],
                            y=df_gc['GC'],
                            mode='lines',
                            name='GC Content',
                            line=dict(color='#1f77b4', width=2),
                            fill='tozeroy',
                            fillcolor='rgba(31, 119, 180, 0.2)'
                        ),
                        row=1, col=1
                    )
                    
                    # Calculate average GC content
                    avg_gc = np.mean(gc_content)
                    
                    # Add a horizontal line for average GC content
                    fig.add_shape(
                        type="line",
                        x0=0,
                        y0=avg_gc,
                        x1=seq_len,
                        y1=avg_gc,
                        line=dict(color="red", width=1, dash="dash"),
                        row=1, col=1
                    )
                    
                    # Add annotation for average GC content
                    fig.add_annotation(
                        x=seq_len * 0.02,
                        y=avg_gc,
                        text=f"Avg: {avg_gc:.2f}%",
                        showarrow=False,
                        font=dict(color="red"),
                        row=1, col=1
                    )
                    
                    # Highlight feature regions if requested
                    if highlight_regions and hasattr(selected_record, 'features'):
                        # Filter for interesting features (genes, CDS, etc.)
                        important_features = [
                            f for f in selected_record.features 
                            if f.type in ['gene', 'CDS', 'rRNA', 'tRNA']
                        ]
                        
                        # Add rectangular highlights for features
                        for feature in important_features[:100]:  # Limit to 100 features to avoid overcrowding
                            start = int(feature.location.start)
                            end = int(feature.location.end)
                            
                            color = get_feature_color(feature.type)
                            alpha = 0.1
                            
                            # Add shape to highlight feature
                            fig.add_shape(
                                type="rect",
                                x0=start,
                                y0=0,
                                x1=end,
                                y1=100,  # GC content is in percentage
                                fillcolor=f"rgba{tuple(list(mcolors.to_rgb(color)) + [alpha])}",
                                line=dict(width=0),
                                layer="below",
                                row=1, col=1
                            )
                    
                    # Add GC skew if requested
                    if show_skew:
                        # Add GC skew trace
                        fig.add_trace(
                            go.Scatter(
                                x=df_gc['Position'],
                                y=df_gc['GC_Skew'],
                                mode='lines',
                                name='GC Skew',
                                line=dict(color='#2ca02c', width=2),
                                fill='tozeroy',
                                fillcolor='rgba(44, 160, 44, 0.2)'
                            ),
                            row=2, col=1
                        )
                        
                        # Add a horizontal line at y=0 for skew
                        fig.add_shape(
                            type="line",
                            x0=0,
                            y0=0,
                            x1=seq_len,
                            y1=0,
                            line=dict(color="black", width=1, dash="dash"),
                            row=2, col=1
                        )
                        
                        # Add cumulative GC skew if requested
                        if show_cumulative:
                            cumulative_skew = np.cumsum(gc_skew)
                            
                            # Normalize to fit on the same scale
                            max_abs_skew = max(abs(np.min(cumulative_skew)), abs(np.max(cumulative_skew)))
                            if max_abs_skew > 0:
                                normalized_skew = cumulative_skew / max_abs_skew
                            else:
                                normalized_skew = cumulative_skew
                            
                            fig.add_trace(
                                go.Scatter(
                                    x=df_gc['Position'],
                                    y=normalized_skew,
                                    mode='lines',
                                    name='Cumulative GC Skew',
                                    line=dict(color='#ff7f0e', width=2)
                                ),
                                row=2, col=1
                            )
                            
                            # Add potential origin and terminus
                            if len(cumulative_skew) > 0:
                                origin_idx = np.argmin(cumulative_skew)
                                terminus_idx = np.argmax(cumulative_skew)
                                
                                origin_pos = df_gc.iloc[origin_idx]['Position']
                                terminus_pos = df_gc.iloc[terminus_idx]['Position']
                                
                                # Add vertical lines for origin and terminus
                                fig.add_vline(
                                    x=origin_pos,
                                    line_width=1,
                                    line_dash="dash",
                                    line_color="blue",
                                    annotation_text="Potential Origin",
                                    annotation_position="top",
                                    row=2, col=1
                                )
                                
                                fig.add_vline(
                                    x=terminus_pos,
                                    line_width=1,
                                    line_dash="dash",
                                    line_color="red",
                                    annotation_text="Potential Terminus",
                                    annotation_position="top",
                                    row=2, col=1
                                )
                    
                    # Update layout
                    height = 400 if not show_skew else 700
                    
                    fig.update_layout(
                        title=f"GC Content Analysis (Window Size: {window_size:,} bp, Step Size: {step_size:,} bp)",
                        height=height,
                        showlegend=True,
                        legend=dict(
                            orientation="h",
                            yanchor="bottom",
                            y=1.02,
                            xanchor="right",
                            x=1
                        )
                    )
                    
                    # Update y-axis titles
                    fig.update_yaxes(title_text="GC Content (%)", row=1, col=1)
                    if show_skew:
                        fig.update_yaxes(title_text="GC Skew (G-C)/(G+C)", row=2, col=1)
                    
                    # Update x-axis title (only for the bottom plot)
                    fig.update_xaxes(title_text="Position (bp)", row=2 if show_skew else 1, col=1)
                    
                    st.plotly_chart(fig, use_container_width=True)
                    
                    # Add explanation
                    st.markdown("""
                    <div class="info-box">
                    <strong>About GC Content Analysis:</strong><br>
                    - <b>GC Content</b>: The percentage of G and C nucleotides in a DNA sequence.
                    - <b>GC Skew</b>: (G-C)/(G+C), which measures the relative abundance of G vs. C on a single strand.
                    - <b>Cumulative GC Skew</b>: The running sum of GC skew values. In bacterial genomes, the minimum and
                      maximum of this curve often correspond to the origin and terminus of replication.
                    
                    GC content and skew analyses are useful for:
                    - Identifying potential origins of replication
                    - Detecting horizontally transferred genomic islands
                    - Studying evolutionary relationships
                    - Analyzing coding vs. non-coding regions
                    </div>
                    """, unsafe_allow_html=True)
                    
                    # Summary statistics
                    st.subheader("Summary Statistics")
                    
                    col_stats1, col_stats2 = st.columns(2)
                    
                    with col_stats1:
                        st.metric("Average GC Content", f"{avg_gc:.2f}%")
                        st.metric("Minimum GC Content", f"{min(gc_content):.2f}%")
                        st.metric("Maximum GC Content", f"{max(gc_content):.2f}%")
                    
                    with col_stats2:
                        if show_skew:
                            st.metric("Average GC Skew", f"{np.mean(gc_skew):.4f}")
                            st.metric("Minimum GC Skew", f"{min(gc_skew):.4f}")
                            st.metric("Maximum GC Skew", f"{max(gc_skew):.4f}")
    
    # 5. Synteny Network Tab
    with tabs[4]:
        st.header("Synteny Network Analysis")
        
        if len(records) < 2:
            st.info("Please upload a file with multiple sequences to analyze synteny.")
        else:
            col1, col2 = st.columns([3, 1])
            
            with col2:
                st.subheader("Network Settings")
                
                # Homology settings for network
                st.write("Homology Settings:")
                min_identity = st.slider(
                    "Minimum identity %:",
                    50, 100, 70,
                    key="network_min_identity"
                )
                min_length = st.slider(
                    "Minimum alignment length:",
                    50, 1000, 100,
                    key="network_min_length"
                )
                
                # Select sequences for network
                sequence_options = [rec.id for rec in records]
                selected_seqs = st.multiselect(
                    "Select sequences for network:",
                    options=sequence_options,
                    default=sequence_options[:min(5, len(sequence_options))]
                )
                
                # Network visualization options
                network_type = st.radio(
                    "Network type:",
                    ["Genome-Genome Network", "Gene-Gene Network"]
                )
                
                if network_type == "Gene-Gene Network":
                    # Feature selection for gene network
                    feature_types = st.multiselect(
                        "Feature types to include:",
                        options=list(set(f.type for rec in records if hasattr(rec, 'features') for f in rec.features)),
                        default=["CDS"] if any(hasattr(rec, 'features') and any(f.type == "CDS" for f in rec.features) for rec in records) else []
                    )
            
            with col1:
                if selected_seqs and len(selected_seqs) >= 2:
                    st.subheader("Synteny Network Visualization")
                    
                    # Get selected records
                    selected_records = [rec for rec in records if rec.id in selected_seqs]
                    
                    if network_type == "Genome-Genome Network":
                        # Find homologous regions between selected sequences
                        homology_data = find_homologous_regions(
                            selected_records, 
                            min_identity=min_identity,
                            min_length=min_length
                        )
                        
                        # Create network using networkx
                        G = nx.Graph()
                        
                        # Add nodes (genomes)
                        for record in selected_records:
                            G.add_node(record.id, size=len(record.seq), type='genome')
                        
                        # Add edges based on homology
                        for hom in homology_data:
                            rec1_id, _, _, rec2_id, _, _, identity = hom
                            
                            # Add or update edge
                            if G.has_edge(rec1_id, rec2_id):
                                G[rec1_id][rec2_id]['weight'] += 1
                                G[rec1_id][rec2_id]['identity'] = max(G[rec1_id][rec2_id]['identity'], identity)
                            else:
                                G.add_edge(rec1_id, rec2_id, weight=1, identity=identity)
                        
                        # Create a plot of the network
                        plt.figure(figsize=(10, 10))
                        
                        # Calculate node sizes based on sequence length
                        sizes = [len(records[sequence_options.index(node)].seq) / 50000 for node in G.nodes()]
                        sizes = [max(10, min(30, s)) for s in sizes]  # Limit size range
                        
                        # Calculate edge widths based on weight
                        widths = [G[u][v]['weight'] * 0.5 for u, v in G.edges()]
                        
                        # Calculate edge colors based on identity
                        edge_colors = [G[u][v]['identity'] / 100 for u, v in G.edges()]
                        
                        # Create positions for nodes using spring layout
                        pos = nx.spring_layout(G, seed=42)
                        
                        # Draw the network
                        nodes = nx.draw_networkx_nodes(
                            G, pos,
                            node_size=sizes,
                            node_color='skyblue',
                            alpha=0.8
                        )
                        
                        edges = nx.draw_networkx_edges(
                            G, pos,
                            width=widths,
                            alpha=0.7,
                            edge_color=edge_colors,
                            edge_cmap=plt.cm.viridis
                        )
                        
                        # Add labels
                        nx.draw_networkx_labels(
                            G, pos,
                            font_size=10,
                            font_weight='bold'
                        )
                        
                        # Add colorbar for edge colors
                        sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(vmin=min_identity/100, vmax=1))
                        sm.set_array([])
                        plt.colorbar(sm, label='Sequence Identity', shrink=0.7)
                        
                        plt.axis('off')
                        plt.tight_layout()
                        
                        st.pyplot(plt)
                        
                        # Add network stats
                        st.subheader("Network Statistics")
                        
                        col_net1, col_net2 = st.columns(2)
                        
                        with col_net1:
                            st.metric("Number of Genomes", len(G.nodes()))
                            st.metric("Number of Connections", len(G.edges()))
                            
                            # Most connected genome
                            if len(G.nodes()) > 0:
                                most_connected = max(G.degree, key=lambda x: x[1])[0]
                                st.metric("Most Connected Genome", most_connected)
                        
                        with col_net2:
                            # Average number of connections
                            avg_degree = sum(dict(G.degree()).values()) / len(G.nodes()) if len(G.nodes()) > 0 else 0
                            st.metric("Avg. Connections per Genome", f"{avg_degree:.2f}")
                            
                            # Average sequence identity
                            avg_identity = sum(G[u][v]['identity'] for u, v in G.edges()) / len(G.edges()) if len(G.edges()) > 0 else 0
                            st.metric("Avg. Sequence Identity", f"{avg_identity:.2f}%")
                        
                        # Display edge table
                        st.subheader("Connections Table")
                        
                        edge_data = []
                        for u, v in G.edges():
                            edge_data.append({
                                'Genome 1': u,
                                'Genome 2': v,
                                'Homologous Regions': G[u][v]['weight'],
                                'Max Identity (%)': G[u][v]['identity']
                            })
                        
                        edge_df = pd.DataFrame(edge_data)
                        st.dataframe(edge_df)
                    
                    else:  # Gene-Gene Network
                        if not feature_types:
                            st.warning("Please select at least one feature type for the gene network.")
                        else:
                            st.info("Creating gene-gene network based on sequence similarity...")
                            
                            # Collect all features of selected types
                            all_features = []
                            for record in selected_records:
                                if hasattr(record, 'features'):
                                    for feature in record.features:
                                        if feature.type in feature_types:
                                            # Extract sequence
                                            feature_seq = feature.extract(record.seq)
                                            
                                            # Get feature name
                                            name = None
                                            for key in ['gene', 'locus_tag', 'product', 'name']:
                                                if key in feature.qualifiers:
                                                    name = feature.qualifiers[key][0]
                                                    break
                                            
                                            if not name:
                                                name = f"{feature.type}_{record.id}_{feature.location.start}..{feature.location.end}"
                                            
                                            all_features.append({
                                                'name': name,
                                                'type': feature.type,
                                                'seq': str(feature_seq),
                                                'genome': record.id,
                                                'start': int(feature.location.start),
                                                'end': int(feature.location.end),
                                                'length': len(feature_seq)
                                            })
                            
                            # For demo purposes, create a mock similarity network
                            # In a real implementation, you would use sequence alignment tools
                            
                            # Create graph
                            G = nx.Graph()
                            
                            # Add nodes for each feature
                            for i, feature in enumerate(all_features):
                                G.add_node(
                                    feature['name'],
                                    type=feature['type'],
                                    genome=feature['genome'],
                                    length=feature['length']
                                )
                            
                            # Add edges based on simulated similarity
                            # In reality, you would use sequence alignment to determine similarity
                            for i in range(len(all_features)):
                                for j in range(i+1, len(all_features)):
                                    # Skip if from same genome
                                    if all_features[i]['genome'] == all_features[j]['genome']:
                                        continue
                                    
                                    # Mock similarity score
                                    # In reality, this would be based on sequence alignment
                                    len_diff = abs(all_features[i]['length'] - all_features[j]['length'])
                                    if len_diff > all_features[i]['length'] * 0.3:  # Skip if length differs by >30%
                                        continue
                                    
                                    # Random similarity score between 60% and 100%
                                    # In a real implementation, this would be from sequence alignment
                                    similarity = random.uniform(60, 100)
                                    
                                    if similarity >= min_identity:
                                        G.add_edge(
                                            all_features[i]['name'],
                                            all_features[j]['name'],
                                            weight=similarity
                                        )
                            
                            # Remove nodes without any connections
                            isolated_nodes = list(nx.isolates(G))
                            G.remove_nodes_from(isolated_nodes)
                            
                            if len(G.nodes()) == 0:
                                st.warning("No significant gene similarities found with the current settings.")
                            else:
                                # Create a plot of the network
                                plt.figure(figsize=(12, 12))
                                
                                # Calculate node sizes based on sequence length
                                sizes = [G.nodes[node]['length'] / 100 for node in G.nodes()]
                                sizes = [max(10, min(50, s)) for s in sizes]  # Limit size range
                                
                                # Calculate node colors based on type
                                types = [G.nodes[node]['type'] for node in G.nodes()]
                                unique_types = list(set(types))
                                type_colors = {t: get_feature_color(t) for t in unique_types}
                                node_colors = [type_colors[G.nodes[node]['type']] for node in G.nodes()]
                                
                                # Calculate edge widths based on similarity
                                widths = [(G[u][v]['weight'] - min_identity) / (100 - min_identity) * 2 + 0.5 for u, v in G.edges()]
                                
                                # Create positions for nodes using spring layout
                                pos = nx.spring_layout(G, seed=42, k=0.2)
                                
                                # Draw the network
                                nodes = nx.draw_networkx_nodes(
                                    G, pos,
                                    node_size=sizes,
                                    node_color=node_colors,
                                    alpha=0.8
                                )
                                
                                edges = nx.draw_networkx_edges(
                                    G, pos,
                                    width=widths,
                                    alpha=0.6,
                                    edge_color='gray'
                                )
                                
                                # Add labels for the largest nodes only (to avoid overcrowding)
                                node_sizes = {node: size for node, size in zip(G.nodes(), sizes)}
                                largest_nodes = sorted(node_sizes.items(), key=lambda x: x[1], reverse=True)[:min(30, len(G.nodes()))]
                                largest_nodes = [node for node, _ in largest_nodes]
                                
                                labeldict = {node: node for node in largest_nodes}
                                
                                nx.draw_networkx_labels(
                                    G, pos,
                                    labels=labeldict,
                                    font_size=8,
                                    font_weight='bold'
                                )
                                
                                # Add legend for node types
                                legend_elements = [
                                    mpatches.Patch(color=type_colors[t], label=t)
                                    for t in unique_types
                                ]
                                plt.legend(handles=legend_elements, loc='upper right')
                                
                                plt.axis('off')
                                plt.tight_layout()
                                
                                st.pyplot(plt)
                                
                                # Network statistics
                                st.subheader("Network Statistics")
                                
                                col_net1, col_net2 = st.columns(2)
                                
                                with col_net1:
                                    st.metric("Number of Genes", len(G.nodes()))
                                    st.metric("Number of Connections", len(G.edges()))
                                    
                                    # Connected components (gene families)
                                    components = list(nx.connected_components(G))
                                    st.metric("Number of Gene Families", len(components))
                                
                                with col_net2:
                                    # Average degree
                                    avg_degree = sum(dict(G.degree()).values()) / len(G.nodes()) if len(G.nodes()) > 0 else 0
                                    st.metric("Avg. Connections per Gene", f"{avg_degree:.2f}")
                                    
                                    # Average similarity
                                    avg_similarity = sum(G[u][v]['weight'] for u, v in G.edges()) / len(G.edges()) if len(G.edges()) > 0 else 0
                                    st.metric("Avg. Sequence Similarity", f"{avg_similarity:.2f}%")
                                    
                                    # Largest gene family
                                    largest_component = max(components, key=len) if components else []
                                    st.metric("Largest Gene Family Size", len(largest_component))
                                
                                # Display the largest gene families
                                st.subheader("Largest Gene Families")
                                
                                # Sort components by size
                                sorted_components = sorted(components, key=len, reverse=True)
                                
                                # Display the top 5 components
                                for i, component in enumerate(sorted_components[:5]):
                                    with st.expander(f"Gene Family {i+1} ({len(component)} genes)"):
                                        # Get the genes in this family
                                        genes = list(component)
                                        
                                        # Create a table with gene information
                                        gene_data = []
                                        for gene in genes:
                                            gene_data.append({
                                                'Gene': gene,
                                                'Type': G.nodes[gene]['type'],
                                                'Genome': G.nodes[gene]['genome'],
                                                'Length': G.nodes[gene]['length']
                                            })
                                        
                                        gene_df = pd.DataFrame(gene_data)
                                        st.dataframe(gene_df)
                                        
                                        # Create a subgraph for this component
                                        subgraph = G.subgraph(genes)
                                        
                                        # Create a plot of the subgraph
                                        plt.figure(figsize=(8, 6))
                                        
                                        # Calculate node sizes and colors
                                        sub_sizes = [G.nodes[node]['length'] / 100 for node in subgraph.nodes()]
                                        sub_sizes = [max(20, min(100, s)) for s in sub_sizes]
                                        
                                        sub_colors = [type_colors[G.nodes[node]['type']] for node in subgraph.nodes()]
                                        
                                        # Calculate edge widths
                                        sub_widths = [(subgraph[u][v]['weight'] - min_identity) / (100 - min_identity) * 2 + 1 for u, v in subgraph.edges()]
                                        
                                        # Create positions
                                        sub_pos = nx.spring_layout(subgraph, seed=42)
                                        
                                        # Draw the network
                                        nx.draw_networkx_nodes(
                                            subgraph, sub_pos,
                                            node_size=sub_sizes,
                                            node_color=sub_colors,
                                            alpha=0.8
                                        )
                                        
                                        nx.draw_networkx_edges(
                                            subgraph, sub_pos,
                                            width=sub_widths,
                                            alpha=0.7,
                                            edge_color='gray'
                                        )
                                        
                                        nx.draw_networkx_labels(
                                            subgraph, sub_pos,
                                            font_size=8,
                                            font_weight='bold'
                                        )
                                        
                                        plt.axis('off')
                                        plt.tight_layout()
                                        
                                        st.pyplot(plt)

# Run the app
if __name__ == "__main__":
    main()
    
