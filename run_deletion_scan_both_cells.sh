#!/bin/bash
#
# Run deletion scan for both mouse cell types
#

set -e  # Exit on error

echo "======================================================================"
echo "Multi-Size Deletion Scan — Both Cell Types"
echo "======================================================================"
echo ""

# Run for Mouse ESC (original)
echo "Running for Mouse ESC (EFO:0004038)..."
python deletion_scan.py jingyun --del-sizes 10 40 80

echo ""
echo "======================================================================"
echo ""

# Run for Olfactory receptor cells
echo "Running for Olfactory receptor cells (CL:0000207)..."
# Temporarily modify the cell type in deletion_scan.py
sed -i.bak "s/CELL_TYPE = 'EFO:0004038'/CELL_TYPE = 'CL:0000207'/" deletion_scan.py
sed -i.bak "s/CELL_NAME = 'Mouse ESC'/CELL_NAME = 'Olfactory receptor cell'/" deletion_scan.py

python deletion_scan.py jingyun --del-sizes 10 40 80

# Restore original
mv deletion_scan.py.bak deletion_scan.py

echo ""
echo "======================================================================"
echo "Both cell types complete!"
echo "Check media/ directory for all plots"
echo "======================================================================"
