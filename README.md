# VCF-Venn-Diagram-Plotter
A Python command-line tool for generating customizable Venn diagrams from VCF comparison TSV files.
______________________________________________________________________________________________________________________________


Use at your own risk. I cannot provide support. All information obtained/inferred with this script is without any implied warranty of fitness for any purpose or use whatsoever.

This program was built to work with the output files from the Annotated-VCF-Compare tool (https://github.com/bjtill/Annotated-VCF-Compare).  It has not been tested with other tables containing variant counts. Also check the Three VCF analysis tool (https://github.com/bjtill/ThreeVCF_Analysis_TVA-GUI) that also creates Venn diagrams. That program requires R and the R package VennDiagram. The issue is that VennDiagram is not avalable in newer versions of R and so the current tool was built. 

## Features

- Works with 2 or 3 variant callers
- Automatically reads caller names from TSV header
- Highly customizable appearance:
  - Circle colors
  - Outline styles (none, solid, dashed, dotted)
  - Diagram rotation
  - Font sizes and families
  - Font styles (normal, bold, italic, bold-italic)
  - Transparency control
- Multiple output formats (PNG, PDF, SVG, JPEG, TIFF, EPS)
- Command-line interface with sensible defaults

## Installation

### Requirements

- Python 3.6 or higher
- Required packages:
  ```bash
  pip install matplotlib matplotlib-venn pandas
  ```

### Quick Install

```bash
# Clone or download the script
chmod +x vcf_venn_plotter.py

# Install dependencies
pip install matplotlib matplotlib-venn pandas
```

## Usage

### Basic Usage

```bash
# Simplest usage with defaults
python vcf_venn_plotter.py -i summary.tsv -o venn_diagram.png

# Specify output format
python vcf_venn_plotter.py -i summary.tsv -o venn_diagram.pdf
```

### Custom Colors

```bash
# Use hex color codes (most reliable)
python vcf_venn_plotter.py -i summary.tsv -o venn.png \
  --color1 "#456f01" \
  --color2 "#00688B" \
  --color3 "#ffac12"
```

### Custom Styling

```bash
# Bold numbers with larger labels
python vcf_venn_plotter.py -i summary.tsv -o venn.png \
  --fontsize 14 \
  --style bold \
  --label-fontsize 16 \
  --label-style bold

# Italic numbers with serif font
python vcf_venn_plotter.py -i summary.tsv -o venn.pdf \
  --font serif \
  --style italic \
  --label-font "Times New Roman"
```

### Outline Styles

```bash
# No outlines
python vcf_venn_plotter.py -i summary.tsv -o venn.png --outline none

# Dashed outlines
python vcf_venn_plotter.py -i summary.tsv -o venn.png --outline dashed

# Dotted outlines with custom width
python vcf_venn_plotter.py -i summary.tsv -o venn.png \
  --outline dotted \
  --outline-width 3.0
```

### Transparency/Opacity

The `--alpha` parameter controls the transparency of the circles. Lower values make circles more transparent, allowing overlap regions to show color blending naturally.

```bash
# More transparent (better color blending in overlaps)
python vcf_venn_plotter.py -i summary.tsv -o venn.png --alpha 0.4

# More opaque
python vcf_venn_plotter.py -i summary.tsv -o venn.png --alpha 0.7

# Fully opaque (no transparency)
python vcf_venn_plotter.py -i summary.tsv -o venn.png --alpha 1.0
```

### High-Resolution Output

```bash
# High DPI for publication
python vcf_venn_plotter.py -i summary.tsv -o venn.png --dpi 600

# Large figure size
python vcf_venn_plotter.py -i summary.tsv -o venn.png --figsize 12 10
```

### Complete Example (Recommended Settings)

```bash
python vcf_venn_plotter.py -i summary.tsv -o venn.png \
  --color1 "#456f01" \
  --color2 "#00688B" \
  --color3 "#ffac12" \
  --outline solid \
  --fontsize 12 \
  --font sans-serif \
  --style normal \
  --label-fontsize 14 \
  --label-font sans-serif \
  --label-style normal \
  --dpi 300
```

## Command-Line Options

### Required Arguments

| Option | Description |
|--------|-------------|
| `-i`, `--input` | Input TSV file from VCF comparison |
| `-o`, `--output` | Output file (formats: png, pdf, svg, jpg, tiff, eps) |

### Color Options

| Option | Default | Description |
|--------|---------|-------------|
| `--color1` | #456f01 | Color for first caller (hex format recommended) |
| `--color2` | #00688B | Color for second caller (hex format recommended) |
| `--color3` | #ffac12 | Color for third caller (hex format recommended) |

**Note:** Use hex color codes (e.g., `#456f01`) for best compatibility. Some R color names are supported but may not always work.

### Outline Options

| Option | Default | Description |
|--------|---------|-------------|
| `--outline` | solid | Style: none, solid, dashed, dotted |
| `--outline-width` | 2.0 | Width of circle outlines |

### Font Options (Numbers)

| Option | Default | Description |
|--------|---------|-------------|
| `--fontsize` | 12 | Font size for numbers |
| `--font` | sans-serif | Font family |
| `--style` | normal | Style: normal, bold, italic, bold-italic |

### Font Options (Labels)

| Option | Default | Description |
|--------|---------|-------------|
| `--label-fontsize` | 14 | Font size for caller labels |
| `--label-font` | sans-serif | Font family for labels |
| `--label-style` | normal | Style: normal, bold, italic, bold-italic |

### Other Options

| Option | Default | Description |
|--------|---------|-------------|
| `--figsize` | 10 8 | Figure size (width height) |
| `--dpi` | 300 | DPI for raster formats |
| `--alpha` | 0.5 | Transparency (0=transparent, 1=opaque). Controls color blending in overlaps |

**Note:** The `--alpha` parameter controls how transparent the circles are. Lower values (e.g., 0.4) create more transparent circles and allow better color blending in overlap regions. Higher values (e.g., 0.7-1.0) make circles more opaque.

## Input File Format

The script expects a TSV file generated by the VCF comparison C++ program with columns including:
- Caller columns (e.g., `Bcftools`, `FreeBayes`, `GATK`) with binary values (0/1)
- The script automatically detects which callers are present
- Works with 2 or 3 callers

## Output

The script generates:
1. A Venn diagram image in your specified format
2. Console output showing variant counts for each region

Example output:
```
Found 3 callers: Bcftools, FreeBayes, GATK
Creating 3-way Venn diagram...
Venn diagram saved to: venn.png

Variant counts:
  Unique to Bcftools: 1523
  Unique to FreeBayes: 892
  Unique to GATK: 1247
  Bcftools & FreeBayes only: 234
  Bcftools & GATK only: 456
  FreeBayes & GATK only: 189
  Common to all three: 3421
```

## Available Fonts

The following fonts are available through the `--font` and `--label-font` options:
- sans-serif (default)
- serif
- monospace
- Arial
- Helvetica
- Times New Roman
- Courier
- Palatino

## Troubleshooting

### matplotlib-venn not found

```bash
pip install matplotlib-venn
```

### Colors not displaying correctly

Use quotes around color names with special characters:
```bash
--color2 "deepskyblue4"
```

### TSV file not found

Ensure you provide the correct path to your TSV file:
```bash
python vcf_venn_plotter.py -i /path/to/summary.tsv -o output.png
```

## Integration with C++ VCF Comparison Tool

This Python script is designed to work seamlessly with the C++ VCF comparison program (see https://github.com/bjtill/Annotated-VCF-Compare):

```bash
# Step 1: Run C++ comparison
./vcf_compare -b bcftools.vcf -f freebayes.vcf -g gatk.vcf -o summary.tsv

# Step 2: Generate Venn diagram
python vcf_venn_plotter.py -i summary.tsv -o venn_diagram.png
```

## Author

Bradley John Till (with assistance from Claude AI)

## Version History

- v1.0 (February 2026) - Initial release
  - Support for 2-3 callers
  - Customizable colors, fonts, and styles
  - Multiple output formats
