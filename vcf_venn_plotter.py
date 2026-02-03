#!/usr/bin/env python3
"""
VCF Venn Diagram Plotter
Generates Venn diagrams from VCF comparison TSV files

Author: Brad Till (with assistance from Claude)
Date: February 2026
License: MIT

Usage:
    python vcf_venn_plotter.py -i summary.tsv -o venn_diagram.png
    python vcf_venn_plotter.py -i summary.tsv -o venn.pdf --color1 "#456f01" --color2 "#00688B" --color3 "#ffac12"
"""

import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3, venn2_circles, venn3_circles
import matplotlib.patches as mpatches

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Generate Venn diagrams from VCF comparison TSV files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with defaults
  python vcf_venn_plotter.py -i summary.tsv -o venn.png
  
  # Custom colors (hex format recommended)
  python vcf_venn_plotter.py -i summary.tsv -o venn.pdf --color1 "#456f01" --color2 "#00688B" --color3 "#ffac12"
  
  # Custom fonts and styles
  python vcf_venn_plotter.py -i summary.tsv -o venn.svg --fontsize 14 --font "Arial" --style "bold"
  
  # No outlines
  python vcf_venn_plotter.py -i summary.tsv -o venn.png --outline "none"
        """
    )
    
    # Required arguments
    parser.add_argument('-i', '--input', required=True,
                        help='Input TSV file from VCF comparison')
    parser.add_argument('-o', '--output', required=True,
                        help='Output file for Venn diagram (supported: png, pdf, svg, jpg, tiff, eps)')
    
    # Color options
    parser.add_argument('--color1', default='#456f01',
                        help='Color for first caller circle (hex format recommended, default: #456f01)')
    parser.add_argument('--color2', default='#00688B',
                        help='Color for second caller circle (hex format recommended, default: #00688B)')
    parser.add_argument('--color3', default='#ffac12',
                        help='Color for third caller circle (hex format recommended, default: #ffac12)')
    
    # Outline options
    parser.add_argument('--outline', choices=['none', 'solid', 'dashed', 'dotted'], 
                        default='solid',
                        help='Style of circle outlines (default: solid)')
    parser.add_argument('--outline-width', type=float, default=2.0,
                        help='Width of circle outlines (default: 2.0)')
    
    # Font options for numbers
    parser.add_argument('--fontsize', type=float, default=12.0,
                        help='Font size for numbers in circles (default: 12)')
    parser.add_argument('--font', default='sans-serif',
                        choices=['sans-serif', 'serif', 'monospace', 'Arial', 'Helvetica', 
                                'Times New Roman', 'Courier', 'Palatino'],
                        help='Font family for numbers (default: sans-serif)')
    parser.add_argument('--style', default='normal',
                        choices=['normal', 'bold', 'italic', 'bold-italic'],
                        help='Font style for numbers (default: normal)')
    
    # Font options for labels
    parser.add_argument('--label-fontsize', type=float, default=14.0,
                        help='Font size for caller labels (default: 14)')
    parser.add_argument('--label-font', default='sans-serif',
                        choices=['sans-serif', 'serif', 'monospace', 'Arial', 'Helvetica',
                                'Times New Roman', 'Courier', 'Palatino'],
                        help='Font family for labels (default: sans-serif)')
    parser.add_argument('--label-style', default='normal',
                        choices=['normal', 'bold', 'italic', 'bold-italic'],
                        help='Font style for labels (default: normal)')
    
    # Figure size
    parser.add_argument('--figsize', type=float, nargs=2, default=[10, 8],
                        metavar=('WIDTH', 'HEIGHT'),
                        help='Figure size in inches (default: 10 8)')
    
    # DPI for raster formats
    parser.add_argument('--dpi', type=int, default=300,
                        help='DPI for raster image formats (default: 300)')
    
    # Alpha (transparency)
    parser.add_argument('--alpha', type=float, default=0.5,
                        help='Transparency/opacity of circles (0=transparent, 1=opaque, default: 0.5). '
                             'Overlapping regions will show color blending.')
    
    return parser.parse_args()

def convert_r_color_to_matplotlib(color):
    """Convert R color names to matplotlib-compatible colors"""
    # Map common R colors to matplotlib equivalents
    r_to_mpl_colors = {
        'deepskyblue4': '#00688B',
        'deepskyblue3': '#009ACD',
        'deepskyblue2': '#00B2EE',
        'deepskyblue1': '#00BFFF',
        'deepskyblue': '#00BFFF',
        'skyblue4': '#4A708B',
        'skyblue3': '#6CA6CD',
        'skyblue2': '#7EC0EE',
        'skyblue1': '#87CEFF',
        'steelblue4': '#36648B',
        'steelblue3': '#4F94CD',
        'steelblue2': '#5CACEE',
        'steelblue1': '#63B8FF',
        'dodgerblue4': '#104E8B',
        'dodgerblue3': '#1874CD',
        'dodgerblue2': '#1C86EE',
        'dodgerblue1': '#1E90FF',
    }
    
    # If it's an R color name, convert it
    if color in r_to_mpl_colors:
        return r_to_mpl_colors[color]
    
    # Otherwise return as-is (could be hex code or matplotlib color name)
    return color

def convert_font_style(style):
    """Convert style string to matplotlib font properties"""
    style_map = {
        'normal': ('normal', 'normal'),
        'bold': ('normal', 'bold'),
        'italic': ('italic', 'normal'),
        'bold-italic': ('italic', 'bold')
    }
    return style_map.get(style, ('normal', 'normal'))

def convert_line_style(outline):
    """Convert outline style to matplotlib linestyle"""
    style_map = {
        'none': '',
        'solid': '-',
        'dashed': '--',
        'dotted': ':'
    }
    return style_map.get(outline, '-')

def read_tsv_and_count_variants(tsv_file):
    """
    Read TSV file and count variants for each caller combination
    Returns: dict with caller names and variant counts
    """
    try:
        df = pd.read_csv(tsv_file, sep='\t')
    except Exception as e:
        print(f"Error reading TSV file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Identify which caller columns exist
    caller_columns = []
    caller_names = []
    
    # Look for caller columns (they should be binary 0/1 columns after TotalCallersIdentified)
    possible_callers = ['Bcftools', 'FreeBayes', 'GATK', 'bcftools', 'freebayes', 'gatk']
    
    for caller in possible_callers:
        if caller in df.columns:
            caller_columns.append(caller)
            caller_names.append(caller)
    
    if len(caller_columns) < 2:
        print("Error: Could not find at least 2 caller columns in TSV file", file=sys.stderr)
        print(f"Available columns: {df.columns.tolist()}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(caller_columns)} callers: {', '.join(caller_names)}")
    
    # Count variants for each combination
    counts = {}
    
    if len(caller_columns) == 2:
        # Two callers
        caller1, caller2 = caller_columns[0], caller_columns[1]
        
        # Variants only in caller 1
        counts['10'] = len(df[(df[caller1] == 1) & (df[caller2] == 0)])
        
        # Variants only in caller 2
        counts['01'] = len(df[(df[caller1] == 0) & (df[caller2] == 1)])
        
        # Variants in both
        counts['11'] = len(df[(df[caller1] == 1) & (df[caller2] == 1)])
        
    elif len(caller_columns) == 3:
        # Three callers
        caller1, caller2, caller3 = caller_columns[0], caller_columns[1], caller_columns[2]
        
        # Variants only in caller 1
        counts['100'] = len(df[(df[caller1] == 1) & (df[caller2] == 0) & (df[caller3] == 0)])
        
        # Variants only in caller 2
        counts['010'] = len(df[(df[caller1] == 0) & (df[caller2] == 1) & (df[caller3] == 0)])
        
        # Variants only in caller 3
        counts['001'] = len(df[(df[caller1] == 0) & (df[caller2] == 0) & (df[caller3] == 1)])
        
        # Variants in caller 1 and 2 only
        counts['110'] = len(df[(df[caller1] == 1) & (df[caller2] == 1) & (df[caller3] == 0)])
        
        # Variants in caller 1 and 3 only
        counts['101'] = len(df[(df[caller1] == 1) & (df[caller2] == 0) & (df[caller3] == 1)])
        
        # Variants in caller 2 and 3 only
        counts['011'] = len(df[(df[caller1] == 0) & (df[caller2] == 1) & (df[caller3] == 1)])
        
        # Variants in all three
        counts['111'] = len(df[(df[caller1] == 1) & (df[caller2] == 1) & (df[caller3] == 1)])
    
    return {
        'caller_names': caller_names,
        'counts': counts,
        'num_callers': len(caller_columns)
    }

def create_venn_diagram(data, args):
    """
    Create and save Venn diagram
    """
    caller_names = data['caller_names']
    counts = data['counts']
    num_callers = data['num_callers']
    
    # Convert font styles
    num_style, num_weight = convert_font_style(args.style)
    label_style, label_weight = convert_font_style(args.label_style)
    
    # Convert line style
    linestyle = convert_line_style(args.outline)
    
    # Convert R colors to matplotlib colors
    color1 = convert_r_color_to_matplotlib(args.color1)
    color2 = convert_r_color_to_matplotlib(args.color2)
    color3 = convert_r_color_to_matplotlib(args.color3)
    
    # Create figure
    fig, ax = plt.subplots(figsize=tuple(args.figsize))
    
    # Check if data is very unbalanced (one set much larger than others)
    # This helps with better visual representation
    normalize_to = None
    if num_callers == 3:
        all_counts = [counts['100'], counts['010'], counts['001'], 
                     counts['110'], counts['101'], counts['011'], counts['111']]
        max_count = max(all_counts) if all_counts else 1
        min_nonzero = min([c for c in all_counts if c > 0]) if any(c > 0 for c in all_counts) else 1
        
        # If one region is more than 10x larger than the smallest, normalize
        if max_count > 10 * min_nonzero:
            normalize_to = 'region'
    
    if num_callers == 2:
        # Two-way Venn diagram
        subsets = (counts['10'], counts['01'], counts['11'])
        colors = [color1, color2]
        
        # Use normalize_to parameter for better visual balance
        v = venn2(subsets=subsets, set_labels=caller_names, ax=ax, 
                 alpha=args.alpha, normalize_to=normalize_to)
        
        # Color the patches
        if v.patches[0]:
            v.patches[0].set_facecolor(colors[0])
        if v.patches[1]:
            v.patches[1].set_facecolor(colors[1])
        
        # Set number font properties
        for text in v.subset_labels:
            if text:
                text.set_fontsize(args.fontsize)
                text.set_family(args.font)
                text.set_style(num_style)
                text.set_weight(num_weight)
        
        # Set label font properties
        for text in v.set_labels:
            if text:
                text.set_fontsize(args.label_fontsize)
                text.set_family(args.label_font)
                text.set_style(label_style)
                text.set_weight(label_weight)
        
        # Draw circles with specified outline
        if linestyle:
            c = venn2_circles(subsets=subsets, ax=ax, linestyle=linestyle, 
                            linewidth=args.outline_width)
    
    elif num_callers == 3:
        # Three-way Venn diagram
        subsets = (counts['100'], counts['010'], counts['110'], 
                  counts['001'], counts['101'], counts['011'], counts['111'])
        
        # Set colors for the three main circles
        # matplotlib-venn will automatically blend colors in overlap regions based on alpha
        v = venn3(subsets=subsets, set_labels=caller_names, 
                 set_colors=(color1, color2, color3), 
                 alpha=args.alpha, ax=ax)
        
        # Set number font properties
        for text in v.subset_labels:
            if text:
                text.set_fontsize(args.fontsize)
                text.set_family(args.font)
                text.set_style(num_style)
                text.set_weight(num_weight)
        
        # Set label font properties
        for text in v.set_labels:
            if text:
                text.set_fontsize(args.label_fontsize)
                text.set_family(args.label_font)
                text.set_style(label_style)
                text.set_weight(label_weight)
        
        # Draw circles with specified outline
        if linestyle:
            c = venn3_circles(subsets=subsets, ax=ax, linestyle=linestyle,
                            linewidth=args.outline_width)
    
    # Note: Rotation is not currently supported due to text positioning issues
    # with matplotlib-venn. The rotation parameter is kept for future compatibility
    # but does not affect the output.
    
    # Remove axis
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Save figure
    plt.tight_layout()
    plt.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
    print(f"Venn diagram saved to: {args.output}")
    
    # Print summary
    print("\nVariant counts:")
    if num_callers == 2:
        print(f"  Unique to {caller_names[0]}: {counts['10']}")
        print(f"  Unique to {caller_names[1]}: {counts['01']}")
        print(f"  Common to both: {counts['11']}")
    elif num_callers == 3:
        print(f"  Unique to {caller_names[0]}: {counts['100']}")
        print(f"  Unique to {caller_names[1]}: {counts['010']}")
        print(f"  Unique to {caller_names[2]}: {counts['001']}")
        print(f"  {caller_names[0]} & {caller_names[1]} only: {counts['110']}")
        print(f"  {caller_names[0]} & {caller_names[2]} only: {counts['101']}")
        print(f"  {caller_names[1]} & {caller_names[2]} only: {counts['011']}")
        print(f"  Common to all three: {counts['111']}")

def main():
    args = parse_arguments()
    
    print("VCF Venn Diagram Plotter")
    print("=" * 50)
    
    # Check if matplotlib-venn is available
    try:
        import matplotlib_venn
    except ImportError:
        print("Error: matplotlib-venn is not installed", file=sys.stderr)
        print("Install it with: pip install matplotlib-venn", file=sys.stderr)
        sys.exit(1)
    
    # Read TSV and count variants
    print(f"Reading input file: {args.input}")
    data = read_tsv_and_count_variants(args.input)
    
    # Create Venn diagram
    print(f"Creating {data['num_callers']}-way Venn diagram...")
    create_venn_diagram(data, args)
    
    print("\nDone!")

if __name__ == "__main__":
    main()
