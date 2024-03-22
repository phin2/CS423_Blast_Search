import tkinter as tk
from tkinter import filedialog
import tempfile

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

import xml.etree.ElementTree as ET 

import pandas as pd

def convertFileToSequence(filename):
    clean = -1
    with open(filename, 'r') as file:
        # read in first line
        header = file.readline().partition(' ')[0] + '\n'
        print(header)
        if (header[0] == '>'):
            print("in FASTA format")
            clean = 1
            # read rest of file
            sequence = file.read()
            # remove all return and newline characters
            sequence = sequence.replace("\r", "")
            sequence = sequence.replace("\n", "")
        else:
            print("invalid format")
    if(clean == -1):
         return -1  # so other functions can check
    return (header,sequence)

def multifasta(files):
    sequences = [convertFileToSequence(file) for file in files]
    multifasta = ""
    for seq in sequences:
        string = ">" + seq[0] + seq[1] + '\n'
        multifasta += string
    
    return multifasta

def xmlparser(xmlfile):
    records = NCBIXML.parse(xmlfile)
    item = next(records)
    for alginment in item.alginments:
        for hsp in alginment.hsp:
            print(alginment.title)

if __name__ == "__main__":
    root = tk.Tk()
    root.withdraw()
    
    #list of pathogenic files
    pathogenic_file = open("pathogenic.txt","r")
    pathogens = pathogenic_file.read()

    pathogens = pathogens.split("\n")

    #prompts user for files to blast
    files = filedialog.askopenfilenames(parent=root,title ="Select Sequence Files")
    
    #formats sequences into multifasta format
    mf = multifasta(files)
    mf = mf.replace(">>",">")


    #performs blast search on the multifasta
    result_handle = NCBIWWW.qblast("blastn", "nt",mf,hitlist_size=20)
    blast_results = result_handle.read()
    
    with open("myBlast.xml", "w") as saveTo:
        saveTo.write(blast_results)
        print(blast_results)
        result_handle.close()

    # open xml file and print its file handle to be sure it was created
    with open("myBlast.xml") as result:
        print(result)

    # go through xml and print the titles that have a low E-value
    for record in NCBIXML.parse(open("myBlast.xml")):
        if record.alignments:
            print("\n")
            print("query: %s" % record.query[:100])
            for align in record.alignments:
                print("match: %s " % align.title[:100])

    

    



    
