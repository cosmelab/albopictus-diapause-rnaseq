#!/usr/bin/env python3
"""
Script to read and extract text from PDF manuscripts in the metadata directory
"""

import os
import PyPDF2
from pathlib import Path

def extract_text_from_pdf(pdf_path):
    """Extract text from a PDF file"""
    try:
        with open(pdf_path, 'rb') as file:
            pdf_reader = PyPDF2.PdfReader(file)
            text = ""
            for page in pdf_reader.pages:
                text += page.extract_text() + "\n"
        return text
    except Exception as e:
        return f"Error reading {pdf_path}: {e}"

def analyze_manuscripts():
    """Analyze all PDF manuscripts in the metadata directory"""
    metadata_dir = Path("data/metadata")
    
    if not metadata_dir.exists():
        print(f"Metadata directory not found: {metadata_dir}")
        return
    
    # Find all PDF files
    pdf_files = list(metadata_dir.glob("*.pdf"))
    
    if not pdf_files:
        print("No PDF files found in metadata directory")
        return
    
    print(f"Found {len(pdf_files)} PDF files:")
    for pdf_file in pdf_files:
        print(f"  - {pdf_file.name}")
    
    print("\n" + "="*80 + "\n")
    
    # Extract text from each PDF and save to .txt files
    for pdf_file in pdf_files:
        print(f"PROCESSING: {pdf_file.name}")
        print("="*60)
        
        text = extract_text_from_pdf(pdf_file)
        
        # Create output filename with .txt extension
        txt_file = pdf_file.with_suffix('.txt')
        
        # Save full text to .txt file
        try:
            with open(txt_file, 'w', encoding='utf-8') as f:
                f.write(text)
            print(f"✓ Saved full text to: {txt_file.name}")
        except Exception as e:
            print(f"✗ Error saving {txt_file.name}: {e}")
        
        # Print first 2000 characters to get an overview
        print("First 2000 characters:")
        print(text[:2000])
        print("\n" + "-"*60 + "\n")
        
        # Look for key sections
        sections_to_find = [
            "methods", "method", "experimental", "design", 
            "results", "differential", "expression", "embryo",
            "timepoint", "photoperiod", "sample", "replicate"
        ]
        
        print("Key sections found:")
        text_lower = text.lower()
        for section in sections_to_find:
            if section in text_lower:
                # Find the context around this keyword
                idx = text_lower.find(section)
                start = max(0, idx - 100)
                end = min(len(text), idx + 200)
                context = text[start:end].replace('\n', ' ')
                print(f"  {section}: {context}")
        
        print("\n" + "="*80 + "\n")

if __name__ == "__main__":
    analyze_manuscripts() 