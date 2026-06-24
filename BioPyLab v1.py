#!/usr/bin/env python3
# BioPyLab.py — Single-file Professional Bioinformatics Toolkit v3.0
# Author: Amir Mahdi Mahboubi
# License: MIT
# Python 3.10+  |  Dependencies: matplotlib, reportlab  (openpyxl optional)
# --------------------------------------------------------------------------

from __future__ import annotations
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import threading
import queue
import re
import csv
import json
import logging
import io
import os
import sys
import time
from collections import Counter, defaultdict
from typing import Dict, List, Tuple, Optional, Any, Callable

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib import colors
from xml.sax.saxutils import escape as xml_escape

# Optional Excel support
try:
    import openpyxl
    HAS_OPENPYXL = True
except ImportError:
    HAS_OPENPYXL = False

# ====================== CONSTANTS & DATA ======================
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}
STOP_CODONS = {'TAA', 'TAG', 'TGA'}
IUPAC_CODES = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT',
    'M': 'AC', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG',
    'N': 'ACGT'
}

# ---------- Kyte-Doolittle hydropathy indices ----------
KD_HYDROPATHY = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

# ---------- Instability index dipeptide weights (Guruprasad et al. 1990) ----------
DIWV = {
    'AA': 1.0, 'AC': 44.94, 'AD': 1.0, 'AE': 1.0, 'AF': 1.0, 'AG': -0.15, 'AH': -8.8, 'AI': 1.0, 'AK': -1.05, 'AL': 1.0, 'AM': 1.0, 'AN': 1.0, 'AP': 20.26, 'AQ': 1.0, 'AR': 1.0, 'AS': 1.0, 'AT': 1.0, 'AV': 1.0, 'AW': 1.0, 'AY': 1.0,
    'CA': 1.0, 'CC': 1.0, 'CD': 1.0, 'CE': 1.0, 'CF': 1.0, 'CG': 1.0, 'CH': 1.0, 'CI': 1.0, 'CK': 1.0, 'CL': 1.0, 'CM': -19.49, 'CN': 1.0, 'CP': 1.0, 'CQ': 1.0, 'CR': 1.0, 'CS': 1.0, 'CT': 1.0, 'CV': 1.0, 'CW': 1.0, 'CY': 1.0,
    'DA': 1.0, 'DC': 1.0, 'DD': 1.0, 'DE': 1.0, 'DF': 1.0, 'DG': 1.0, 'DH': 1.0, 'DI': 1.0, 'DK': 1.0, 'DL': 1.0, 'DM': 1.0, 'DN': 1.0, 'DP': 1.0, 'DQ': 1.0, 'DR': 1.0, 'DS': 1.0, 'DT': 1.0, 'DV': 1.0, 'DW': -11.48, 'DY': 1.0,
    'EA': 1.0, 'EC': 1.0, 'ED': 1.0, 'EE': 1.0, 'EF': 1.0, 'EG': 1.0, 'EH': 1.0, 'EI': 1.0, 'EK': 1.0, 'EL': 1.0, 'EM': 1.0, 'EN': 1.0, 'EP': 1.0, 'EQ': 1.0, 'ER': 1.0, 'ES': 1.0, 'ET': 1.0, 'EV': 1.0, 'EW': 1.0, 'EY': 1.0,
    'FA': 1.0, 'FC': 1.0, 'FD': 13.34, 'FE': 1.0, 'FF': 1.0, 'FG': 1.0, 'FH': 1.0, 'FI': 1.0, 'FK': -0.02, 'FL': 1.0, 'FM': 1.0, 'FN': 1.0, 'FP': 1.0, 'FQ': 1.0, 'FR': 1.0, 'FS': 1.0, 'FT': 1.0, 'FV': 1.0, 'FW': 1.0, 'FY': 33.48,
    'GA': -0.89, 'GC': 1.0, 'GD': 1.0, 'GE': 1.0, 'GF': 1.0, 'GG': -7.8, 'GH': 1.0, 'GI': 1.0, 'GK': 1.0, 'GL': 1.0, 'GM': 1.0, 'GN': -7.55, 'GP': 1.0, 'GQ': 1.0, 'GR': 1.0, 'GS': 1.0, 'GT': 1.0, 'GV': 1.0, 'GW': 1.0, 'GY': 1.0,
    'HA': 1.0, 'HC': 1.0, 'HD': 1.0, 'HE': 1.0, 'HF': 1.0, 'HG': -0.41, 'HH': 1.0, 'HI': 1.0, 'HK': 1.0, 'HL': 1.0, 'HM': 1.0, 'HN': 1.0, 'HP': 1.0, 'HQ': 1.0, 'HR': 1.0, 'HS': 1.0, 'HT': 1.0, 'HV': 1.0, 'HW': 1.0, 'HY': 1.0,
    'IA': 1.0, 'IC': 1.0, 'ID': 1.0, 'IE': 1.0, 'IF': 1.0, 'IG': 1.0, 'IH': 1.0, 'II': 1.0, 'IK': 1.0, 'IL': 1.0, 'IM': 1.0, 'IN': 1.0, 'IP': 1.0, 'IQ': 1.0, 'IR': 1.0, 'IS': 1.0, 'IT': 1.0, 'IV': 1.0, 'IW': 1.0, 'IY': 1.0,
    'KA': 1.0, 'KC': 1.0, 'KD': 1.0, 'KE': 1.0, 'KF': 1.0, 'KG': -4.14, 'KH': 1.0, 'KI': 1.0, 'KK': 1.0, 'KL': 1.0, 'KM': 1.0, 'KN': 1.0, 'KP': 1.0, 'KQ': 1.0, 'KR': 1.0, 'KS': 1.0, 'KT': 1.0, 'KV': 1.0, 'KW': 1.0, 'KY': 1.0,
    'LA': 1.0, 'LC': 1.0, 'LD': 1.0, 'LE': 1.0, 'LF': 1.0, 'LG': 1.0, 'LH': 1.0, 'LI': 1.0, 'LK': 1.0, 'LL': 1.0, 'LM': 1.0, 'LN': 1.0, 'LP': 1.0, 'LQ': 1.0, 'LR': 1.0, 'LS': 1.0, 'LT': 1.0, 'LV': 1.0, 'LW': 1.0, 'LY': 1.0,
    'MA': 1.0, 'MC': 1.0, 'MD': 1.0, 'ME': 1.0, 'MF': 1.0, 'MG': 1.0, 'MH': 1.0, 'MI': 1.0, 'MK': 1.0, 'ML': 1.0, 'MM': 1.0, 'MN': 1.0, 'MP': 1.0, 'MQ': 1.0, 'MR': 1.0, 'MS': 1.0, 'MT': 1.0, 'MV': 1.0, 'MW': 1.0, 'MY': 1.0,
    'NA': 1.0, 'NC': 1.0, 'ND': 1.0, 'NE': 1.0, 'NF': 1.0, 'NG': -1.04, 'NH': 1.0, 'NI': 1.0, 'NK': 1.0, 'NL': 1.0, 'NM': 1.0, 'NN': 1.0, 'NP': 1.0, 'NQ': 1.0, 'NR': 1.0, 'NS': 1.0, 'NT': 1.0, 'NV': 1.0, 'NW': 1.0, 'NY': 1.0,
    'PA': 1.0, 'PC': 1.0, 'PD': 1.0, 'PE': 1.0, 'PF': 1.0, 'PG': 1.0, 'PH': 1.0, 'PI': 1.0, 'PK': 1.0, 'PL': 1.0, 'PM': 1.0, 'PN': 1.0, 'PP': 1.0, 'PQ': 1.0, 'PR': 1.0, 'PS': 1.0, 'PT': 1.0, 'PV': 1.0, 'PW': 1.0, 'PY': 1.0,
    'QA': 1.0, 'QC': 1.0, 'QD': 1.0, 'QE': 1.0, 'QF': 1.0, 'QG': 1.0, 'QH': 1.0, 'QI': 1.0, 'QK': 1.0, 'QL': 1.0, 'QM': 1.0, 'QN': 1.0, 'QP': 1.0, 'QQ': 1.0, 'QR': 1.0, 'QS': 1.0, 'QT': 1.0, 'QV': 1.0, 'QW': 1.0, 'QY': 1.0,
    'RA': 1.0, 'RC': 1.0, 'RD': 1.0, 'RE': 1.0, 'RF': 1.0, 'RG': 1.0, 'RH': 1.0, 'RI': 1.0, 'RK': 1.0, 'RL': 1.0, 'RM': 1.0, 'RN': 1.0, 'RP': 1.0, 'RQ': 1.0, 'RR': 1.0, 'RS': 1.0, 'RT': 1.0, 'RV': 1.0, 'RW': 1.0, 'RY': 1.0,
    'SA': 1.0, 'SC': 1.0, 'SD': 1.0, 'SE': 1.0, 'SF': 1.0, 'SG': 1.0, 'SH': 1.0, 'SI': 1.0, 'SK': 1.0, 'SL': 1.0, 'SM': 1.0, 'SN': 1.0, 'SP': 1.0, 'SQ': 1.0, 'SR': 1.0, 'SS': 1.0, 'ST': 1.0, 'SV': 1.0, 'SW': 1.0, 'SY': 1.0,
    'TA': 1.0, 'TC': 1.0, 'TD': 1.0, 'TE': 1.0, 'TF': 1.0, 'TG': 1.0, 'TH': 1.0, 'TI': 1.0, 'TK': 1.0, 'TL': 1.0, 'TM': 1.0, 'TN': 1.0, 'TP': 1.0, 'TQ': 1.0, 'TR': 1.0, 'TS': 1.0, 'TT': 1.0, 'TV': 1.0, 'TW': 1.0, 'TY': 1.0,
    'VA': 1.0, 'VC': 1.0, 'VD': 1.0, 'VE': 1.0, 'VF': 1.0, 'VG': 1.0, 'VH': 1.0, 'VI': 1.0, 'VK': 1.0, 'VL': 1.0, 'VM': 1.0, 'VN': 1.0, 'VP': 1.0, 'VQ': 1.0, 'VR': 1.0, 'VS': 1.0, 'VT': 1.0, 'VV': 1.0, 'VW': 1.0, 'VY': 1.0,
    'WA': 1.0, 'WC': 1.0, 'WD': 1.0, 'WE': 1.0, 'WF': 1.0, 'WG': 1.0, 'WH': 1.0, 'WI': 1.0, 'WK': 1.0, 'WL': 1.0, 'WM': 1.0, 'WN': 1.0, 'WP': 1.0, 'WQ': 1.0, 'WR': 1.0, 'WS': 1.0, 'WT': 1.0, 'WV': 1.0, 'WW': 1.0, 'WY': 1.0,
    'YA': 1.0, 'YC': 1.0, 'YD': 1.0, 'YE': 1.0, 'YF': 1.0, 'YG': 1.0, 'YH': 1.0, 'YI': 1.0, 'YK': 1.0, 'YL': 1.0, 'YM': 1.0, 'YN': 1.0, 'YP': 1.0, 'YQ': 1.0, 'YR': 1.0, 'YS': 1.0, 'YT': 1.0, 'YV': 1.0, 'YW': 1.0, 'YY': 1.0,
}
def _dipeptide_weight(dipep: str) -> float:
    return DIWV.get(dipep, 1.0)

# ---------- Chou-Fasman propensities ----------
CHOU_FASMAN = {
    'A': (1.42, 0.83), 'R': (0.98, 0.93), 'N': (0.67, 0.89), 'D': (0.73, 0.54),
    'C': (0.70, 1.19), 'E': (1.39, 0.50), 'Q': (1.17, 0.75), 'G': (0.43, 0.75),
    'H': (1.05, 0.87), 'I': (1.22, 1.60), 'L': (1.21, 1.30), 'K': (1.00, 0.74),
    'M': (1.45, 1.05), 'F': (1.33, 1.38), 'P': (0.57, 0.55), 'S': (0.79, 0.75),
    'T': (0.82, 1.19), 'W': (1.14, 1.37), 'Y': (0.69, 1.47), 'V': (1.14, 1.70)
}

# ---------- Amino acid properties for pI ----------
AA_PKA = {
    'A': 0.0, 'R': 12.48, 'N': 0.0, 'D': 3.65, 'C': 8.18, 'E': 4.25,
    'Q': 0.0, 'G': 0.0, 'H': 6.00, 'I': 0.0, 'L': 0.0, 'K': 10.53,
    'M': 0.0, 'F': 0.0, 'P': 0.0, 'S': 0.0, 'T': 0.0, 'W': 0.0,
    'Y': 10.07, 'V': 0.0
}
N_TERM_PKA = 8.0
C_TERM_PKA = 3.1

AA_MW = {
    'A': 71.08, 'R': 156.19, 'N': 114.11, 'D': 115.09, 'C': 103.15, 'E': 129.12,
    'Q': 128.13, 'G': 57.05, 'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.20, 'F': 147.18, 'P': 97.12, 'S': 87.08, 'T': 101.11, 'W': 186.21,
    'Y': 163.18, 'V': 99.13
}

# ---------- Restriction Enzymes (100+) ----------
ENZYMES = [
    ("AatII",   "GACGTC", 5), ("Acc65I",  "GGTACC", 1), ("AflII",   "CTTAAG", 1),
    ("AgeI",    "ACCGGT", 1), ("ApaI",    "GGGCCC", 5), ("ApoI",    "RAATTY", 1),
    ("AscI",    "GGCGCGCC", 2), ("AseI",    "ATTAAT", 2), ("AsiSI",   "GCGATCGC", 4),
    ("BamHI",   "GGATCC", 1), ("BclI",    "TGATCA", 2), ("BglII",   "AGATCT", 1),
    ("BmtI",    "GCTAGC", 5), ("BsaI",    "GGTCTC", 1), ("BsaHI",   "GRCGYC", 1),
    ("BsaWI",   "WCCGGW", 1), ("BsiWI",   "CGTACG", 1), ("BsmBI",   "CGTCTC", 1),
    ("BspDI",   "ATCGAT", 2), ("BspHI",   "TCATGA", 1), ("BsrGI",   "TGTACA", 1),
    ("BstBI",   "TTCGAA", 2), ("BstEII",  "GGTNACC", 1), ("BstXI",   "CCANNNNNNTGG", 6),
    ("Bsu36I",  "CCTNAGG", 4), ("ClaI",    "ATCGAT", 2), ("DraI",    "TTTAAA", 3),
    ("EagI",    "CGGCCG", 1), ("EcoRI",   "GAATTC", 1), ("EcoRV",   "GATATC", 3),
    ("FseI",    "GGCCGGCC", 6), ("HhaI",    "GCGC", 2), ("HindIII", "AAGCTT", 1),
    ("HpaI",    "GTTAAC", 3), ("KasI",    "GGCGCC", 1), ("KpnI",    "GGTACC", 5),
    ("MfeI",    "CAATTG", 1), ("MluI",    "ACGCGT", 2), ("MscI",    "TGGCCA", 3),
    ("NcoI",    "CCATGG", 2), ("NdeI",    "CATATG", 2), ("NheI",    "GCTAGC", 1),
    ("NotI",    "GCGGCCGC", 2), ("NruI",    "TCGCGA", 3), ("NsiI",    "ATGCAT", 5),
    ("PacI",    "TTAATTAA", 5), ("PciI",    "ACATGT", 2), ("PmeI",    "GTTTAAAC", 4),
    ("PmlI",    "CACGTG", 3), ("PstI",    "CTGCAG", 5), ("PvuI",    "CGATCG", 4),
    ("PvuII",   "CAGCTG", 3), ("SacI",    "GAGCTC", 5), ("SacII",   "CCGCGG", 4),
    ("SalI",    "GTCGAC", 1), ("SbfI",    "CCTGCAGG", 6), ("ScaI",    "AGTACT", 3),
    ("SfoI",    "GGCGCC", 1), ("SmaI",    "CCCGGG", 3), ("SnaBI",   "TACGTA", 3),
    ("SpeI",    "ACTAGT", 1), ("SphI",    "GCATGC", 5), ("SspI",    "AATATT", 3),
    ("StuI",    "AGGCCT", 3), ("StyI",    "CCWWGG", 1), ("SwaI",    "ATTTAAAT", 4),
    ("TspMI",   "CCCGGG", 3), ("XbaI",    "TCTAGA", 1), ("XhoI",    "CTCGAG", 1),
    ("XmaI",    "CCCGGG", 1), ("ZraI",    "GACGTC", 5),
    ("BsaAI",   "YACGTR", 3), ("BseRI",   "GAGGAG", 2), ("BsgI",    "GTGCAG", 2),
    ("BsmI",    "GAATGC", 1), ("BspMI",   "ACCTGC", 4), ("EarI",    "CTCTTC", 1),
    ("FokI",    "GGATG", 9), ("HgaI",    "GACGC", 5), ("MboII",   "GAAGA", 8),
    ("PflMI",   "CCANNNNNTGG", 6), ("SapI",    "GCTCTTC", 1), ("Tth111I", "GACNNNGTC", 2),
    ("AarI",    "CACCTGC", 4), ("AcuI",    "CTGAAG", 4), ("AlwI",    "GGATC", 4),
    ("BbsI",    "GAAGAC", 2), ("BciVI",   "GTATCC", 5), ("BfuAI",   "ACCTGC", 4),
    ("BmrI",    "ACTGGG", 2), ("BpmI",    "CTGGAG", 4), ("BpuEI",   "CTTGAG", 4),
    ("BsaXI",   "ACNNNNNCTCC", 2), ("BseYI",   "CCCAGC", 4), ("BsmFI",   "GGGAC", 5),
    ("BspQI",   "GCTCTTC", 1), ("BsrDI",   "GCAATG", 2), ("EciI",    "GGCGGA", 2),
    ("Esp3I",   "CGTCTC", 1), ("FauI",    "CCCGC", 4), ("Hin4I",   "GAYNNNNNVTC", 5),
    ("HphI",    "GGTGA", 8), ("MlyI",    "GAGTC", 5), ("MmeI",    "TCCRAC", 6),
    ("MnlI",    "CCTC", 6), ("NmeAIII", "GCCGAG", 5), ("PleI",    "GAGTC", 4),
    ("PsrI",    "GAACNNNNNNTAC", 2), ("SfaNI",   "GCATC", 5), ("TspRI",   "CAGTG", 4),
]

# ====================== HELPER CLASSES ======================
class LogCapture(logging.Handler):
    def __init__(self, text_widget: scrolledtext.ScrolledText):
        super().__init__()
        self.text_widget = text_widget

    def emit(self, record: logging.LogRecord):
        msg = self.format(record)
        def append():
            self.text_widget.insert(tk.END, msg + '\n')
            self.text_widget.see(tk.END)
        if self.text_widget.winfo_exists():
            self.text_widget.after(0, append)

class ThreadedTask:
    def __init__(self, app, target: Callable, args=(), kwargs=None,
                 progress_callback: Optional[Callable] = None):
        self.app = app
        self.target = target
        self.args = args
        self.kwargs = kwargs or {}
        self.progress_callback = progress_callback
        self.thread = threading.Thread(target=self._run, daemon=True)
        self.queue = queue.Queue()

    def _run(self):
        try:
            result = self.target(*self.args, **self.kwargs)
            self.queue.put(("success", result))
        except Exception as e:
            logging.exception("Task error")
            self.queue.put(("error", str(e)))

    def start(self, on_complete: Callable = None):
        def check():
            try:
                msg, val = self.queue.get(block=False)
            except queue.Empty:
                self.app.root.after(100, check)
                return
            if msg == "success":
                if on_complete:
                    on_complete(val)
            else:
                messagebox.showerror("Error", f"Task failed: {val}")
            self.app.set_status("Ready")
        self.app.set_status("Running...")
        self.thread.start()
        self.app.root.after(100, check)

# ====================== MAIN APPLICATION ======================
class BioPyLabApp:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("BioPyLab v3.0 – Professional Bioinformatics Toolkit")
        self.root.geometry("1400x900")
        self.root.minsize(1100, 700)

        self.dark_mode = tk.BooleanVar(value=False)
        self.current_theme = "light"

        self.sequences: List[Tuple[str, str]] = []
        self.active_seq_index = 0
        self.analysis_results: Dict[str, Any] = {}
        self.log_stream = io.StringIO()

        self._setup_logging()
        self._setup_styles()
        self._build_toolbar()
        self._build_main_notebook()
        self._build_statusbar()
        self._apply_theme()

    def _setup_logging(self):
        self.logger = logging.getLogger("BioPyLab")
        self.logger.setLevel(logging.DEBUG)
        ch = logging.StreamHandler(sys.stdout)
        ch.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
        self.logger.addHandler(ch)

    def _setup_styles(self):
        style = ttk.Style()
        style.theme_use('clam')
        self._configure_theme(style, dark=False)

    def _configure_theme(self, style, dark: bool):
        if dark:
            bg = '#2e2e2e'
            fg = '#ffffff'
            entry_bg = '#3c3c3c'
            select_bg = '#4a6984'
        else:
            bg = '#f0f0f0'
            fg = '#000000'
            entry_bg = '#ffffff'
            select_bg = '#0078d7'
        style.configure('.', background=bg, foreground=fg, fieldbackground=entry_bg, selectbackground=select_bg)
        style.configure('TLabel', background=bg, foreground=fg)
        style.configure('TButton', background=bg, foreground=fg)
        style.configure('TNotebook', background=bg, borderwidth=0)
        style.configure('TNotebook.Tab', background=bg, foreground=fg, padding=[12, 4], font=('Segoe UI', 9, 'bold'))
        style.map('TNotebook.Tab', background=[('selected', select_bg)], foreground=[('selected', 'white')])
        style.configure('TFrame', background=bg)
        style.configure('TProgressbar', background=select_bg)
        self.root.tk_setPalette(background=bg, foreground=fg, selectBackground=select_bg,
                                selectForeground='white', activeBackground=select_bg)
        self.root.option_add('*Text.Background', entry_bg)
        self.root.option_add('*Text.Foreground', fg)

    def _apply_theme(self):
        dark = self.dark_mode.get()
        style = ttk.Style()
        self._configure_theme(style, dark)
        self.current_theme = "dark" if dark else "light"
        if dark:
            plt.style.use('dark_background')
        else:
            plt.style.use('default')
        self.set_status(f"Theme switched to {self.current_theme} mode")

    def _build_toolbar(self):
        toolbar = ttk.Frame(self.root)
        toolbar.pack(side=tk.TOP, fill=tk.X, padx=5, pady=4)
        ttk.Button(toolbar, text="Load Sequence(s)", command=self.load_sequences).pack(side=tk.LEFT, padx=2)
        ttk.Button(toolbar, text="Clear All", command=self.clear_all).pack(side=tk.LEFT, padx=2)
        ttk.Button(toolbar, text="Export Results", command=self.export_current_results).pack(side=tk.LEFT, padx=2)
        ttk.Separator(toolbar, orient=tk.VERTICAL).pack(side=tk.LEFT, padx=8, fill=tk.Y)
        ttk.Checkbutton(toolbar, text="Dark Mode", variable=self.dark_mode, command=self._apply_theme).pack(side=tk.LEFT, padx=5)
        ttk.Label(toolbar, text="   BioPyLab v3.0", font=('Segoe UI', 10, 'italic')).pack(side=tk.RIGHT, padx=10)

    def _build_main_notebook(self):
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.tab_input = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_input, text="Input & Sequences")
        self._build_input_tab()

        self.tab_basic = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_basic, text="Basic Analysis")
        self._build_basic_tab()

        self.tab_orf = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_orf, text="ORF Finder")
        self._build_orf_tab()

        self.tab_mut = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_mut, text="Mutations")
        self._build_mutation_tab()

        self.tab_kmer = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_kmer, text="K-mer Analysis")
        self._build_kmer_tab()

        self.tab_codon = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_codon, text="Codon Usage")
        self._build_codon_tab()

        self.tab_enz = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_enz, text="Restriction Map")
        self._build_enzyme_tab()

        self.tab_aln = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_aln, text="Alignment")
        self._build_alignment_tab()

        self.tab_prot = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_prot, text="Protein Analysis")
        self._build_protein_tab()

        self.tab_charts = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_charts, text="Charts")
        self._build_charts_tab()

        self.tab_log = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_log, text="Log & Report")
        self._build_log_tab()

    def _build_statusbar(self):
        self.status_var = tk.StringVar(value="Ready")
        status_bar = ttk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W)
        status_bar.pack(side=tk.BOTTOM, fill=tk.X)

    def set_status(self, text: str):
        self.status_var.set(text)
        self.root.update_idletasks()

    # -------------------- Tab Builders --------------------
    def _build_input_tab(self):
        frm = ttk.Frame(self.tab_input)
        frm.pack(fill=tk.BOTH, expand=True)
        btnf = ttk.Frame(frm)
        btnf.pack(fill=tk.X, padx=5, pady=5)
        ttk.Button(btnf, text="Open FASTA / Multi-FASTA", command=self.load_fasta).pack(side=tk.LEFT, padx=2)
        ttk.Button(btnf, text="Open FASTQ", command=self.load_fastq).pack(side=tk.LEFT, padx=2)
        ttk.Button(btnf, text="Open GenBank", command=self.load_genbank).pack(side=tk.LEFT, padx=2)
        ttk.Button(btnf, text="Validate Sequences", command=self.validate_sequences).pack(side=tk.LEFT, padx=2)

        paned = ttk.PanedWindow(frm, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        left = ttk.Frame(paned)
        right = ttk.Frame(paned)
        paned.add(left, weight=1)
        paned.add(right, weight=3)

        ttk.Label(left, text="Loaded Sequences").pack(anchor=tk.W)
        self.seq_listbox = tk.Listbox(left, selectmode=tk.SINGLE, exportselection=False)
        self.seq_listbox.pack(fill=tk.BOTH, expand=True, pady=5)
        self.seq_listbox.bind('<<ListboxSelect>>', self._on_seq_select)

        ttk.Label(right, text="Selected Sequence").pack(anchor=tk.W)
        self.seq_display = scrolledtext.ScrolledText(right, wrap=tk.NONE, font=('Courier', 10))
        self.seq_display.pack(fill=tk.BOTH, expand=True)

    def _build_basic_tab(self):
        frm = ttk.Frame(self.tab_basic)
        frm.pack(fill=tk.BOTH, expand=True)
        ttk.Button(frm, text="Run Basic Analysis", command=lambda: self._run_threaded(self.basic_analysis)).pack(pady=5)
        self.result_basic = scrolledtext.ScrolledText(frm, wrap=tk.WORD, font=('Consolas', 10))
        self.result_basic.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    def _build_orf_tab(self):
        frm = ttk.Frame(self.tab_orf)
        frm.pack(fill=tk.BOTH, expand=True)
        top = ttk.Frame(frm)
        top.pack(fill=tk.X, padx=5, pady=5)
        ttk.Label(top, text="Min ORF length (bp):").pack(side=tk.LEFT)
        self.min_orf_var = tk.IntVar(value=30)
        ttk.Entry(top, textvariable=self.min_orf_var, width=8).pack(side=tk.LEFT, padx=5)
        ttk.Button(top, text="Find ORFs (6-frame)", command=lambda: self._run_threaded(self.orf_analysis)).pack(side=tk.LEFT, padx=5)
        ttk.Button(top, text="Export ORFs to CSV", command=self.export_orfs_csv).pack(side=tk.LEFT, padx=5)
        self.result_orf = scrolledtext.ScrolledText(frm, wrap=tk.WORD, font=('Consolas', 10))
        self.result_orf.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    def _build_mutation_tab(self):
        frm = ttk.Frame(self.tab_mut)
        frm.pack(fill=tk.BOTH, expand=True)
        ttk.Button(frm, text="Compare Selected Sequences", command=lambda: self._run_threaded(self.mutation_analysis)).pack(pady=5)
        self.result_mut = scrolledtext.ScrolledText(frm, wrap=tk.WORD, font=('Consolas', 10))
        self.result_mut.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    def _build_kmer_tab(self):
        frm = ttk.Frame(self.tab_kmer)
        frm.pack(fill=tk.BOTH, expand=True)
        top = ttk.Frame(frm)
        top.pack(fill=tk.X, padx=5, pady=5)
        ttk.Label(top, text="K value:").pack(side=tk.LEFT)
        self.k_val_var = tk.IntVar(value=4)
        ttk.Entry(top, textvariable=self.k_val_var, width=5).pack(side=tk.LEFT, padx=5)
        ttk.Button(top, text="Run K-mer Analysis", command=lambda: self._run_threaded(self.kmer_analysis)).pack(side=tk.LEFT, padx=5)
        ttk.Button(top, text="Show K-mer Frequency Plot", command=self.plot_kmer_freq).pack(side=tk.LEFT, padx=5)
        self.result_kmer = scrolledtext.ScrolledText(frm, wrap=tk.WORD, font=('Consolas', 10))
        self.result_kmer.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    def _build_codon_tab(self):
        frm = ttk.Frame(self.tab_codon)
        frm.pack(fill=tk.BOTH, expand=True)
        ttk.Button(frm, text="Show Codon Usage Table", command=lambda: self._run_threaded(self.codon_usage)).pack(pady=5)
        ttk.Button(frm, text="Plot Codon Frequency", command=self.plot_codon_freq).pack(pady=5)
        self.result_codon = scrolledtext.ScrolledText(frm, wrap=tk.WORD, font=('Consolas', 10))
        self.result_codon.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    def _build_enzyme_tab(self):
        frm = ttk.Frame(self.tab_enz)
        frm.pack(fill=tk.BOTH, expand=True)
        ttk.Button(frm, text="Restriction Map", command=lambda: self._run_threaded(self.restriction_map)).pack(pady=5)
        self.result_enz = scrolledtext.ScrolledText(frm, wrap=tk.WORD, font=('Consolas', 10))
        self.result_enz.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    def _build_alignment_tab(self):
        frm = ttk.Frame(self.tab_aln)
        frm.pack(fill=tk.BOTH, expand=True)
        top = ttk.Frame(frm)
        top.pack(fill=tk.X, padx=5, pady=5)
        ttk.Label(top, text="Match:").pack(side=tk.LEFT)
        self.match_var = tk.IntVar(value=1)
        ttk.Entry(top, textvariable=self.match_var, width=5).pack(side=tk.LEFT, padx=5)
        ttk.Label(top, text="Mismatch:").pack(side=tk.LEFT)
        self.mismatch_var = tk.IntVar(value=-1)
        ttk.Entry(top, textvariable=self.mismatch_var, width=5).pack(side=tk.LEFT, padx=5)
        ttk.Label(top, text="Gap:").pack(side=tk.LEFT)
        self.gap_var = tk.IntVar(value=-2)
        ttk.Entry(top, textvariable=self.gap_var, width=5).pack(side=tk.LEFT, padx=5)
        ttk.Button(top, text="Global (NW)", command=lambda: self._run_threaded(self.global_alignment)).pack(side=tk.LEFT, padx=5)
        ttk.Button(top, text="Local (SW)", command=lambda: self._run_threaded(self.local_alignment)).pack(side=tk.LEFT, padx=5)
        self.result_aln = scrolledtext.ScrolledText(frm, wrap=tk.NONE, font=('Courier', 10))
        self.result_aln.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    def _build_protein_tab(self):
        frm = ttk.Frame(self.tab_prot)
        frm.pack(fill=tk.BOTH, expand=True)
        ttk.Label(frm, text="Protein sequence (or use translated from current DNA)").pack(anchor=tk.W, padx=5, pady=5)
        self.prot_seq_input = scrolledtext.ScrolledText(frm, height=4, wrap=tk.WORD)
        self.prot_seq_input.pack(fill=tk.X, padx=5, pady=5)
        btnf = ttk.Frame(frm)
        btnf.pack(fill=tk.X, padx=5, pady=5)
        ttk.Button(btnf, text="Analyze Protein", command=lambda: self._run_threaded(self.protein_analysis)).pack(side=tk.LEFT, padx=2)
        ttk.Button(btnf, text="Predict 2° Structure", command=lambda: self._run_threaded(self.predict_secondary_structure)).pack(side=tk.LEFT, padx=2)
        self.result_prot = scrolledtext.ScrolledText(frm, wrap=tk.WORD, font=('Consolas', 10))
        self.result_prot.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    def _build_charts_tab(self):
        frm = ttk.Frame(self.tab_charts)
        frm.pack(fill=tk.BOTH, expand=True)
        left = ttk.Frame(frm)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)
        buttons = [
            ("Nucleotide Composition", self.plot_nuc_pie),
            ("GC Window", self.plot_gc_window),
            ("GC Skew", self.plot_gc_skew),
            ("K-mer Frequency", self.plot_kmer_freq),
            ("Codon Frequency", self.plot_codon_freq),
            ("Amino Acid Composition", self.plot_aa_composition),
            ("Hydrophobicity Plot", self.plot_hydrophobicity),
            ("Mutation Distribution", self.plot_mutation_dist),
            ("ORF Length Distribution", self.plot_orf_len_dist),
        ]
        for text, cmd in buttons:
            ttk.Button(left, text=text, command=cmd).pack(pady=2, fill=tk.X)
        self.chart_frame = ttk.Frame(frm)
        self.chart_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

    def _build_log_tab(self):
        frm = ttk.Frame(self.tab_log)
        frm.pack(fill=tk.BOTH, expand=True)
        ttk.Button(frm, text="Generate Full Report", command=self.generate_full_report).pack(pady=5)
        ttk.Label(frm, text="Application Log").pack(anchor=tk.W, padx=5)
        self.log_view = scrolledtext.ScrolledText(frm, height=12, wrap=tk.WORD, font=('Consolas', 9))
        self.log_view.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        handler = LogCapture(self.log_view)
        handler.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
        self.logger.addHandler(handler)

    # ---------- Threading helper ----------
    def _run_threaded(self, target, on_complete=None):
        task = ThreadedTask(self, target, on_complete=on_complete)
        task.start()

    # ---------- Unified Sequence Loading ----------
    def load_sequences(self):
        path = filedialog.askopenfilename(
            filetypes=[("All Supported", "*.fasta *.fa *.fna *.txt *.fastq *.fq *.gb *.gbk"), ("All files", "*.*")]
        )
        if not path:
            return
        ext = os.path.splitext(path)[1].lower()
        if ext in ('.fasta', '.fa', '.fna', '.txt'):
            self.load_fasta(path)
        elif ext in ('.fastq', '.fq'):
            self.load_fastq(path)
        elif ext in ('.gb', '.gbk'):
            self.load_genbank(path)
        else:
            # fallback to FASTA
            self.load_fasta(path)

    def load_fasta(self, filepath: str = None):
        if filepath is None:
            filepath = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa *.fna *.txt"), ("All files", "*.*")])
        if not filepath:
            return
        try:
            with open(filepath, 'r') as f:
                content = f.read()
            self._parse_fasta(content)
            self.set_status(f"Loaded {len(self.sequences)} sequences from FASTA")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load FASTA: {e}")

    def _parse_fasta(self, text: str):
        entries = []
        current_header = "Untitled"
        current_seq = []
        for line in text.splitlines():
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    entries.append((current_header, ''.join(current_seq)))
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_seq:
            entries.append((current_header, ''.join(current_seq)))
        self.sequences = entries
        self._update_seq_listbox()

    def load_fastq(self, filepath: str = None):
        if filepath is None:
            filepath = filedialog.askopenfilename(filetypes=[("FASTQ files", "*.fastq *.fq"), ("All files", "*.*")])
        if not filepath:
            return
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
            entries = []
            i = 0
            while i < len(lines):
                header = lines[i].strip()[1:]
                seq = lines[i+1].strip().upper()
                i += 4
                entries.append((header, seq))
            self.sequences = entries
            self._update_seq_listbox()
            self.set_status(f"Loaded {len(entries)} sequences from FASTQ")
        except Exception as e:
            messagebox.showerror("Error", f"FASTQ parse error: {e}")

    def load_genbank(self, filepath: str = None):
        if filepath is None:
            filepath = filedialog.askopenfilename(filetypes=[("GenBank files", "*.gb *.gbk"), ("All files", "*.*")])
        if not filepath:
            return
        try:
            from Bio import SeqIO
            records = list(SeqIO.parse(filepath, "genbank"))
            entries = [(rec.id, str(rec.seq).upper()) for rec in records]
            self.sequences = entries
            self._update_seq_listbox()
            self.set_status(f"Loaded {len(entries)} sequences from GenBank (via Biopython)")
        except ImportError:
            messagebox.showerror("Error", "Biopython is required for GenBank support. Please install it: pip install biopython")
        except Exception as e:
            messagebox.showerror("Error", f"GenBank parse error: {e}")

    def _update_seq_listbox(self):
        self.seq_listbox.delete(0, tk.END)
        for header, _ in self.sequences:
            self.seq_listbox.insert(tk.END, header[:80])
        if self.sequences:
            self.seq_listbox.selection_set(0)
            self._on_seq_select()

    def _on_seq_select(self, event=None):
        sel = self.seq_listbox.curselection()
        if sel:
            idx = sel[0]
            self.active_seq_index = idx
            header, seq = self.sequences[idx]
            self.seq_display.delete('1.0', tk.END)
            self.seq_display.insert('1.0', seq[:5000])

    def clear_all(self):
        self.sequences.clear()
        self._update_seq_listbox()
        for txt in [self.result_basic, self.result_orf, self.result_mut, self.result_kmer,
                    self.result_codon, self.result_enz, self.result_aln, self.result_prot]:
            txt.delete('1.0', tk.END)
        self.set_status("Cleared")

    def validate_sequences(self):
        if not self.sequences:
            messagebox.showinfo("Validation", "No sequences loaded.")
            return
        report = ""
        for i, (header, seq) in enumerate(self.sequences):
            invalid = [c for c in seq if c not in IUPAC_CODES]
            if invalid:
                report += f"{header}: contains invalid characters {set(invalid)}\n"
            else:
                report += f"{header}: OK (len={len(seq)})\n"
        messagebox.showinfo("Validation Report", report)

    # ---------- Core Analysis Methods ----------
    @staticmethod
    def _clean_seq(seq: str, ambiguous: bool = False) -> str:
        if ambiguous:
            return ''.join(c for c in seq.upper() if c in IUPAC_CODES)
        return ''.join(c for c in seq.upper() if c in 'ATGC')

    def _get_active_seq(self, ambiguous=False):
        if not self.sequences:
            return ""
        _, seq = self.sequences[self.active_seq_index]
        return self._clean_seq(seq, ambiguous)

    def _get_two_seqs(self):
        if len(self.sequences) < 2:
            return None, None
        s1 = self._clean_seq(self.sequences[0][1])
        s2 = self._clean_seq(self.sequences[1][1])
        return s1, s2

    def basic_analysis(self):
        seq = self._get_active_seq()
        if not seq:
            return
        counts = Counter(seq)
        length = len(seq)
        gc = (counts['G'] + counts['C']) / length * 100
        rna = seq.replace('T', 'U')
        protein = self._translate(seq)
        rev = self._reverse_complement(seq)
        report = f"""=== BASIC DNA ANALYSIS ===
Length: {length:,} bp
GC Content: {gc:.2f}%

Nucleotide counts: A={counts['A']:,} T={counts['T']:,} G={counts['G']:,} C={counts['C']:,}

RNA (first 300 bases): {rna[:300]}
Reverse complement (first 500): {rev[:500]}
Protein translation (first 500 aa): {protein[:500]}
"""
        self.result_basic.delete('1.0', tk.END)
        self.result_basic.insert('1.0', report)
        self.analysis_results['basic'] = report

    @staticmethod
    def _translate(seq: str) -> str:
        protein = []
        for i in range(0, len(seq)-2, 3):
            codon = seq[i:i+3]
            protein.append(CODON_TABLE.get(codon, 'X'))
        return ''.join(protein)

    @staticmethod
    def _reverse_complement(seq: str) -> str:
        trans = str.maketrans('ATGC', 'TACG')
        return seq.translate(trans)[::-1]

    def orf_analysis(self):
        seq = self._get_active_seq()
        if not seq:
            return
        min_len = self.min_orf_var.get()
        orfs = []
        for frame in range(3):
            s = seq[frame:]
            orfs.extend(self._find_orfs(s, frame+1, min_len))
        rev = self._reverse_complement(seq)
        for frame in range(3):
            s = rev[frame:]
            orfs_neg = self._find_orfs(s, -(frame+1), min_len)
            for (start, end, length, prot) in orfs_neg:
                orig_start = len(seq) - (frame + end) + 1
                orig_end = len(seq) - (frame + start) + 1
                orfs.append((orig_start, orig_end, length, prot, f"-{frame+1}"))
        self.analysis_results['orfs'] = orfs
        self.result_orf.delete('1.0', tk.END)
        out = f"ORFs found (min {min_len} bp): {len(orfs)}\n\n"
        for idx, (s, e, l, prot, frame) in enumerate(orfs[:100], 1):
            out += f"ORF{idx:4d}: {frame:>3s}  {s:7d}..{e:<7d} len={l:5d}  prot={prot[:50]}\n"
        if len(orfs) > 100:
            out += f"\n... and {len(orfs)-100} more."
        self.result_orf.insert('1.0', out)

    def _find_orfs(self, s: str, frame_label, min_len):
        orfs = []
        i = 0
        n = len(s)
        while i < n-2:
            if s[i:i+3] == 'ATG':
                j = i+3
                while j < n-2:
                    codon = s[j:j+3]
                    if codon in STOP_CODONS:
                        length = j+3 - i
                        if length >= min_len:
                            prot = self._translate(s[i:j+3])
                            orfs.append((i+1, j+3, length, prot))
                        break
                    j += 3
                i = j+1 if j < n-2 else i+1
            else:
                i += 1
        return orfs

    def export_orfs_csv(self):
        orfs = self.analysis_results.get('orfs', [])
        if not orfs:
            messagebox.showinfo("Export", "No ORFs to export.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if not path:
            return
        with open(path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["ORF", "Frame", "Start", "End", "Length(bp)", "Protein"])
            for idx, (s, e, l, prot, frame) in enumerate(orfs, 1):
                writer.writerow([idx, frame, s, e, l, prot])
        messagebox.showinfo("Export", f"Exported {len(orfs)} ORFs to {path}")

    def mutation_analysis(self):
        s1, s2 = self._get_two_seqs()
        if not s1 or not s2:
            messagebox.showwarning("Warning", "Load at least two sequences for comparison.")
            return
        n = min(len(s1), len(s2))
        snps = []
        for i in range(n):
            if s1[i] != s2[i]:
                snps.append((i+1, s1[i], s2[i]))
        insertions = max(0, len(s2) - len(s1))
        deletions = max(0, len(s1) - len(s2))
        report = f"""=== MUTATION ANALYSIS ===
Seq1 length: {len(s1)}, Seq2 length: {len(s2)}
SNPs: {len(snps)} | Insertions: {insertions} | Deletions: {deletions}\n"""
        for pos, a, b in snps[:200]:
            report += f"Pos {pos}: {a} -> {b}\n"
        if len(snps) > 200:
            report += f"... {len(snps)-200} more SNPs"
        self.result_mut.delete('1.0', tk.END)
        self.result_mut.insert('1.0', report)
        self.analysis_results['mutations'] = {'snps': snps, 'insertions': insertions, 'deletions': deletions}

    def kmer_analysis(self):
        seq = self._get_active_seq()
        k = self.k_val_var.get()
        if k < 1:
            return
        kmers = [seq[i:i+k] for i in range(len(seq)-k+1)]
        counts = Counter(kmers)
        total = sum(counts.values())
        top = counts.most_common(50)
        report = f"K-mer analysis (k={k}), distinct: {len(counts)}\n"
        for kmer, cnt in top:
            report += f"{kmer}: {cnt} ({cnt/total*100:.2f}%)\n"
        self.result_kmer.delete('1.0', tk.END)
        self.result_kmer.insert('1.0', report)
        self.analysis_results['kmers'] = {'k': k, 'counts': counts}

    def codon_usage(self):
        seq = self._get_active_seq()
        if len(seq) < 3:
            return
        codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
        counts = Counter(codons)
        total = sum(counts.values())
        out = "Codon Usage Table\n" + "-"*40 + "\n"
        for codon in sorted(CODON_TABLE.keys()):
            cnt = counts.get(codon, 0)
            freq = cnt/total*100 if total else 0
            aa = CODON_TABLE[codon]
            out += f"{codon} {aa} : {cnt:6d} ({freq:5.2f}%)\n"
        self.result_codon.delete('1.0', tk.END)
        self.result_codon.insert('1.0', out)
        self.analysis_results['codons'] = {'counts': counts, 'total': total}

    def restriction_map(self):
        seq = self._get_active_seq()
        out = "Restriction Enzyme Map (exact matches only)\n" + "="*60 + "\n"
        for name, site, cut in ENZYMES:
            if any(c not in 'ATGC' for c in site):
                continue
            positions = [m.start() for m in re.finditer(site, seq)]
            if positions:
                out += f"{name:10s} ({site:15s}): {len(positions)} cut(s) at {positions[:10]}\n"
            else:
                out += f"{name:10s} ({site:15s}): no cuts\n"
        self.result_enz.delete('1.0', tk.END)
        self.result_enz.insert('1.0', out)

    def global_alignment(self):
        s1, s2 = self._get_two_seqs()
        if not s1 or not s2:
            return
        match = self.match_var.get()
        mismatch = self.mismatch_var.get()
        gap = self.gap_var.get()
        score, al1, al2, al3 = self._needleman_wunsch(s1, s2, match, mismatch, gap)
        identity = sum(1 for a,b in zip(al1, al3) if a==b) / max(len(al1),1) * 100
        out = f"Global Alignment (NW)  Score: {score}  Identity: {identity:.2f}%\n\n"
        out += self._format_alignment(al1, al2, al3)
        self.result_aln.delete('1.0', tk.END)
        self.result_aln.insert('1.0', out)

    def local_alignment(self):
        s1, s2 = self._get_two_seqs()
        if not s1 or not s2:
            return
        match = self.match_var.get()
        mismatch = self.mismatch_var.get()
        gap = self.gap_var.get()
        score, al1, al2, al3 = self._smith_waterman(s1, s2, match, mismatch, gap)
        identity = sum(1 for a,b in zip(al1, al3) if a==b) / max(len(al1),1) * 100
        out = f"Local Alignment (SW)  Score: {score}  Identity: {identity:.2f}%\n\n"
        out += self._format_alignment(al1, al2, al3)
        self.result_aln.delete('1.0', tk.END)
        self.result_aln.insert('1.0', out)

    @staticmethod
    def _needleman_wunsch(seq1, seq2, match, mismatch, gap):
        m, n = len(seq1), len(seq2)
        dp = [[0]*(n+1) for _ in range(m+1)]
        for i in range(m+1): dp[i][0] = i*gap
        for j in range(n+1): dp[0][j] = j*gap
        for i in range(1,m+1):
            for j in range(1,n+1):
                diag = dp[i-1][j-1] + (match if seq1[i-1]==seq2[j-1] else mismatch)
                up = dp[i-1][j] + gap
                left = dp[i][j-1] + gap
                dp[i][j] = max(diag, up, left)
        al1, al2, al3 = [], [], []
        i, j = m, n
        while i>0 or j>0:
            if i>0 and j>0 and dp[i][j] == dp[i-1][j-1] + (match if seq1[i-1]==seq2[j-1] else mismatch):
                al1.append(seq1[i-1]); al2.append('|' if seq1[i-1]==seq2[j-1] else ' '); al3.append(seq2[j-1])
                i-=1; j-=1
            elif i>0 and dp[i][j] == dp[i-1][j] + gap:
                al1.append(seq1[i-1]); al2.append(' '); al3.append('-')
                i-=1
            else:
                al1.append('-'); al2.append(' '); al3.append(seq2[j-1])
                j-=1
        return dp[m][n], ''.join(reversed(al1)), ''.join(reversed(al2)), ''.join(reversed(al3))

    @staticmethod
    def _smith_waterman(seq1, seq2, match, mismatch, gap):
        m, n = len(seq1), len(seq2)
        dp = [[0]*(n+1) for _ in range(m+1)]
        max_score = 0
        max_pos = (0,0)
        for i in range(1,m+1):
            for j in range(1,n+1):
                diag = dp[i-1][j-1] + (match if seq1[i-1]==seq2[j-1] else mismatch)
                up = dp[i-1][j] + gap
                left = dp[i][j-1] + gap
                dp[i][j] = max(0, diag, up, left)
                if dp[i][j] > max_score:
                    max_score = dp[i][j]
                    max_pos = (i,j)
        i, j = max_pos
        al1, al2, al3 = [], [], []
        while i>0 and j>0 and dp[i][j] > 0:
            if dp[i][j] == dp[i-1][j-1] + (match if seq1[i-1]==seq2[j-1] else mismatch):
                al1.append(seq1[i-1]); al2.append('|' if seq1[i-1]==seq2[j-1] else ' '); al3.append(seq2[j-1])
                i-=1; j-=1
            elif dp[i][j] == dp[i-1][j] + gap:
                al1.append(seq1[i-1]); al2.append(' '); al3.append('-')
                i-=1
            elif dp[i][j] == dp[i][j-1] + gap:
                al1.append('-'); al2.append(' '); al3.append(seq2[j-1])
                j-=1
            else:
                break
        return max_score, ''.join(reversed(al1)), ''.join(reversed(al2)), ''.join(reversed(al3))

    @staticmethod
    def _format_alignment(al1, al2, al3, line_len=80):
        out = ""
        for i in range(0, max(len(al1), len(al3)), line_len):
            out += al1[i:i+line_len] + "\n"
            out += al2[i:i+line_len] + "\n"
            out += al3[i:i+line_len] + "\n\n"
        return out

    def protein_analysis(self):
        prot = self.prot_seq_input.get('1.0', tk.END).strip().upper()
        prot = ''.join(aa for aa in prot if aa in AA_MW)
        if not prot:
            dna = self._get_active_seq()
            if dna:
                prot = self._translate(dna)[:1000]
                self.prot_seq_input.delete('1.0', tk.END)
                self.prot_seq_input.insert('1.0', prot)
        if not prot:
            messagebox.showerror("Error", "No protein sequence.")
            return
        comp = Counter(prot)
        length = len(prot)
        mw = sum(AA_MW[aa] for aa in prot) + 18.02
        pi = self._calc_pi(prot)
        gravy = sum(KD_HYDROPATHY[aa] for aa in prot) / length
        instab = self._instability_index(prot)
        aro = (comp.get('F',0) + comp.get('W',0) + comp.get('Y',0)) / length
        ala = comp.get('A',0); val = comp.get('V',0); ile = comp.get('I',0); leu = comp.get('L',0)
        aliphatic = ala + 2.9*val + 3.9*ile + 3.9*leu
        report = f"""=== PROTEIN ANALYSIS ===
Length: {length} aa
MW: {mw:.2f} Da
pI: {pi:.2f}
GRAVY: {gravy:.3f}
Instability Index: {instab:.2f} (stable < 40)
Aromaticity: {aro:.3f}
Aliphatic Index: {aliphatic:.2f}

Amino Acid Composition:
"""
        for aa in sorted(AA_MW.keys(), key=lambda x: comp.get(x,0), reverse=True):
            cnt = comp.get(aa,0)
            if cnt>0:
                report += f"{aa}: {cnt:5d} ({cnt/length*100:5.2f}%)\n"
        self.result_prot.delete('1.0', tk.END)
        self.result_prot.insert('1.0', report)
        self.analysis_results['protein'] = {'mw': mw, 'pi': pi, 'gravy': gravy, 'instability': instab,
                                           'aromaticity': aro, 'aliphatic': aliphatic, 'composition': comp}

    def _calc_pi(self, prot):
        best_pi = 7.0
        min_charge = float('inf')
        for pH in [x/100.0 for x in range(0, 1401, 1)]:
            charge = 0
            for aa in prot:
                pka = AA_PKA.get(aa, None)
                if pka:
                    if aa in ('D','E','C','Y'):
                        if pH > pka: charge -= 1
                    elif aa in ('R','K','H'):
                        if pH < pka: charge += 1
            if pH < N_TERM_PKA: charge += 1
            if pH > C_TERM_PKA: charge -= 1
            if abs(charge) < min_charge:
                min_charge = abs(charge)
                best_pi = pH
        return best_pi

    def _instability_index(self, prot):
        if len(prot) < 2:
            return 0
        total = 0.0
        for i in range(len(prot)-1):
            di = prot[i] + prot[i+1]
            total += _dipeptide_weight(di)
        return (10.0 / len(prot)) * total

    def predict_secondary_structure(self):
        prot = self.prot_seq_input.get('1.0', tk.END).strip().upper()
        prot = ''.join(aa for aa in prot if aa in CHOU_FASMAN)
        if len(prot) < 6:
            messagebox.showerror("Error", "Protein too short")
            return
        helix = [CHOU_FASMAN[aa][0] for aa in prot]
        sheet = [CHOU_FASMAN[aa][1] for aa in prot]
        pred = []
        i=0; n=len(prot)
        while i < n:
            if i+5 < n and sum(1 for j in range(i,i+6) if helix[j] > 1.0) >= 4:
                pred.append('H'); i+=1
                while i < n and sum(helix[max(0,i-3):i+1])/4 > 1.0:
                    if i >= len(pred): pred.append('H')
                    else: pred[i] = 'H'
                    i+=1
            elif i+4 < n and sum(1 for j in range(i,i+5) if sheet[j] > 1.0) >= 3:
                pred.append('E'); i+=1
                while i < n and sum(sheet[max(0,i-3):i+1])/4 > 1.0:
                    if i >= len(pred): pred.append('E')
                    else: pred[i] = 'E'
                    i+=1
            else:
                pred.append('C'); i+=1
        pred_str = ''.join(pred[:len(prot)])
        out = f"Secondary Structure (Chou-Fasman):\n{prot}\n{pred_str}"
        self.result_prot.delete('1.0', tk.END)
        self.result_prot.insert('1.0', out)

    # ---------- Chart plotting ----------
    def _clear_chart(self):
        for w in self.chart_frame.winfo_children():
            w.destroy()

    def _draw_chart(self, fig):
        self._clear_chart()
        canvas = FigureCanvasTkAgg(fig, master=self.chart_frame)
        canvas.draw()
        toolbar = NavigationToolbar2Tk(canvas, self.chart_frame)
        toolbar.update()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def plot_nuc_pie(self):
        seq = self._get_active_seq()
        if not seq: return
        cnt = Counter(seq)
        labels = ['A','T','G','C']
        sizes = [cnt[x] for x in labels]
        fig = Figure(figsize=(5,4))
        ax = fig.add_subplot(111)
        ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
        ax.set_title("Nucleotide Composition")
        self._draw_chart(fig)

    def plot_gc_window(self):
        seq = self._get_active_seq()
        if len(seq) < 100: return
        window = 100
        gc = []
        for i in range(0, len(seq)-window+1, 10):
            sub = seq[i:i+window]
            gc.append((sub.count('G')+sub.count('C'))/window*100)
        fig = Figure(figsize=(6,3))
        ax = fig.add_subplot(111)
        ax.plot(range(0, len(seq)-window+1, 10), gc)
        ax.set_xlabel("Position"); ax.set_ylabel("GC %"); ax.set_title("Sliding Window GC Content")
        self._draw_chart(fig)

    def plot_gc_skew(self):
        seq = self._get_active_seq()
        if len(seq) < 100: return
        window = 100
        skew = []
        for i in range(0, len(seq)-window+1, 10):
            sub = seq[i:i+window]
            g = sub.count('G'); c = sub.count('C')
            s = (g-c)/(g+c) if (g+c)>0 else 0
            skew.append(s)
        fig = Figure(figsize=(6,3))
        ax = fig.add_subplot(111)
        ax.plot(range(0, len(seq)-window+1, 10), skew)
        ax.set_xlabel("Position"); ax.set_ylabel("GC Skew"); ax.set_title("GC Skew")
        self._draw_chart(fig)

    def plot_kmer_freq(self):
        kmers_data = self.analysis_results.get('kmers')
        if not kmers_data:
            messagebox.showinfo("Info", "Run K-mer analysis first.")
            return
        counts = kmers_data['counts']
        top = counts.most_common(20)
        fig = Figure(figsize=(6,4))
        ax = fig.add_subplot(111)
        kmers, freqs = zip(*top)
        ax.bar(kmers, freqs)
        ax.set_title("Top K-mer Frequencies")
        ax.tick_params(axis='x', rotation=45)
        fig.tight_layout()
        self._draw_chart(fig)

    def plot_codon_freq(self):
        cdata = self.analysis_results.get('codons')
        if not cdata:
            messagebox.showinfo("Info", "Run Codon Usage first.")
            return
        counts = cdata['counts']
        codons = sorted(CODON_TABLE.keys())
        freqs = [counts.get(c,0) for c in codons]
        fig = Figure(figsize=(8,4))
        ax = fig.add_subplot(111)
        ax.bar(codons, freqs)
        ax.set_title("Codon Frequency")
        ax.tick_params(axis='x', rotation=90, labelsize=7)
        fig.tight_layout()
        self._draw_chart(fig)

    def plot_aa_composition(self):
        prot_data = self.analysis_results.get('protein')
        if not prot_data:
            messagebox.showinfo("Info", "Run Protein Analysis first.")
            return
        comp = prot_data['composition']
        aas = sorted(comp.keys())
        vals = [comp[aa] for aa in aas]
        fig = Figure(figsize=(6,4))
        ax = fig.add_subplot(111)
        ax.bar(aas, vals)
        ax.set_title("Amino Acid Composition")
        self._draw_chart(fig)

    def plot_hydrophobicity(self):
        prot_seq = self._get_active_seq_translated()
        if not prot_seq: return
        hydro = [KD_HYDROPATHY.get(aa, 0) for aa in prot_seq]
        fig = Figure(figsize=(6,3))
        ax = fig.add_subplot(111)
        ax.plot(hydro)
        ax.set_xlabel("Residue"); ax.set_ylabel("Hydropathy"); ax.set_title("Hydrophobicity Plot (Kyte-Doolittle)")
        self._draw_chart(fig)

    def plot_mutation_dist(self):
        mut = self.analysis_results.get('mutations')
        if not mut:
            messagebox.showinfo("Info", "Run Mutation Analysis first.")
            return
        snps = mut['snps']
        if not snps: return
        positions = [p for p,_,_ in snps]
        fig = Figure(figsize=(6,3))
        ax = fig.add_subplot(111)
        ax.hist(positions, bins=50)
        ax.set_title("Mutation Distribution")
        ax.set_xlabel("Position"); ax.set_ylabel("Count")
        self._draw_chart(fig)

    def plot_orf_len_dist(self):
        orfs = self.analysis_results.get('orfs')
        if not orfs:
            messagebox.showinfo("Info", "Run ORF Finder first.")
            return
        lengths = [l for _,_,l,_,_ in orfs]
        fig = Figure(figsize=(6,3))
        ax = fig.add_subplot(111)
        ax.hist(lengths, bins=50)
        ax.set_title("ORF Length Distribution")
        ax.set_xlabel("Length (bp)"); ax.set_ylabel("Count")
        self._draw_chart(fig)

    def _get_active_seq_translated(self):
        seq = self._get_active_seq()
        if seq:
            return self._translate(seq)
        return ""

    # ---------- Export functions ----------
    def export_current_results(self):
        formats = [("PDF", ".pdf"), ("CSV", ".csv"), ("JSON", ".json"), ("HTML", ".html")]
        if HAS_OPENPYXL: formats.append(("Excel", ".xlsx"))
        path = filedialog.asksaveasfilename(defaultextension=".pdf", filetypes=[(f"{n} files", f"*{ext}") for n,ext in formats])
        if not path: return
        ext = os.path.splitext(path)[1].lower()
        if ext == '.pdf':
            self._export_pdf(path)
        elif ext == '.csv':
            self._export_csv(path)
        elif ext == '.json':
            self._export_json(path)
        elif ext == '.html':
            self._export_html(path)
        elif ext == '.xlsx':
            self._export_excel(path)

    def _export_pdf(self, path):
        text = self._get_current_result_text()
        doc = SimpleDocTemplate(path)
        styles = getSampleStyleSheet()
        safe = xml_escape(text).replace('\n', '<br/>')
        doc.build([Paragraph("BioPyLab Report", styles['Title']), Spacer(1,12), Paragraph(safe, styles['BodyText'])])
        messagebox.showinfo("Export", "PDF exported.")

    def _export_csv(self, path):
        orfs = self.analysis_results.get('orfs')
        if orfs:
            with open(path, 'w', newline='') as f:
                w = csv.writer(f)
                w.writerow(["ORF","Frame","Start","End","Length","Protein"])
                for idx, (s,e,l,prot,frame) in enumerate(orfs,1):
                    w.writerow([idx,frame,s,e,l,prot])
            messagebox.showinfo("Export", "ORFs exported as CSV.")
        else:
            messagebox.showinfo("Export", "No tabular data available.")

    def _export_json(self, path):
        with open(path, 'w') as f:
            json.dump(self.analysis_results, f, indent=2, default=str)
        messagebox.showinfo("Export", "Analysis data exported as JSON.")

    def _export_html(self, path):
        text = self._get_current_result_text()
        html = f"<html><head><title>BioPyLab Report</title></head><body><pre>{text}</pre></body></html>"
        with open(path, 'w') as f:
            f.write(html)
        messagebox.showinfo("Export", "HTML exported.")

    def _export_excel(self, path):
        if not HAS_OPENPYXL:
            messagebox.showerror("Error", "openpyxl not installed.")
            return
        wb = openpyxl.Workbook()
        ws = wb.active
        ws.title = "Results"
        orfs = self.analysis_results.get('orfs')
        if orfs:
            ws.append(["ORF","Frame","Start","End","Length","Protein"])
            for idx, (s,e,l,prot,frame) in enumerate(orfs,1):
                ws.append([idx,frame,s,e,l,prot])
        else:
            ws['A1'] = "No tabular data"
        wb.save(path)
        messagebox.showinfo("Export", "Excel exported.")

    def _get_current_result_text(self):
        for i in range(self.notebook.index('end')):
            tab = self.notebook.nametowidget(self.notebook.tabs()[i])
            for child in tab.winfo_children():
                if isinstance(child, scrolledtext.ScrolledText):
                    return child.get('1.0', tk.END)
        return ""

    def generate_full_report(self):
        path = filedialog.asksaveasfilename(defaultextension=".pdf", filetypes=[("PDF", "*.pdf"), ("HTML", "*.html")])
        if not path: return
        ext = os.path.splitext(path)[1].lower()
        if ext == '.pdf':
            self._full_report_pdf(path)
        elif ext == '.html':
            self._full_report_html(path)

    def _full_report_pdf(self, path):
        doc = SimpleDocTemplate(path)
        styles = getSampleStyleSheet()
        story = [Paragraph("BioPyLab Comprehensive Report", styles['Title']), Spacer(1,12)]
        for section_name, key in [("Basic Analysis","basic"), ("ORFs","orfs"), ("Mutations","mutations")]:
            text = self.analysis_results.get(key, "")
            if isinstance(text, str):
                safe = xml_escape(text).replace('\n','<br/>')
                story.append(Paragraph(f"<b>{section_name}</b>", styles['Heading2']))
                story.append(Paragraph(safe, styles['BodyText']))
        doc.build(story)
        messagebox.showinfo("Report", "Full PDF report generated.")

    def _full_report_html(self, path):
        html = "<html><head><title>BioPyLab Report</title></head><body>"
        html += "<h1>BioPyLab Comprehensive Report</h1>"
        for section_name, key in [("Basic Analysis","basic"), ("ORFs","orfs"), ("Mutations","mutations")]:
            text = self.analysis_results.get(key, "")
            if isinstance(text, str):
                html += f"<h2>{section_name}</h2><pre>{text}</pre>"
        html += "</body></html>"
        with open(path, 'w') as f:
            f.write(html)
        messagebox.showinfo("Report", "Full HTML report generated.")

# ====================== MAIN ======================
if __name__ == "__main__":
    app = BioPyLabApp()
    app.root.mainloop()