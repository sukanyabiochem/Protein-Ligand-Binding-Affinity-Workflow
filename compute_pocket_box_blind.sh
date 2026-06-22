#!/bin/bash
# compute_pocket_box_blind.sh
#
# Defines a Vina search box (docking.center_x/y/z, docking.size_x/y/z)
# when NO co-crystallized ligand is available, using one of three modes:
#
#   whole     Bounding box of the ENTIRE receptor (blind docking).
#             Usage: bash compute_pocket_box_blind.sh whole receptor.pdb [padding]
#
#   residues  Centroid + bounding box of a user-specified list of
#             active-site residue numbers (from literature/UniProt
#             active-site annotation), optionally restricted to one chain.
#             Usage: bash compute_pocket_box_blind.sh residues receptor.pdb "57,102,195" [padding] [chain]
#
#   fpocket   Run fpocket to detect cavities, then report the centroid +
#             bounding box of the TOP-ranked pocket's alpha spheres.
#             Requires: conda install -c bioconda fpocket
#             Usage: bash compute_pocket_box_blind.sh fpocket receptor.pdb [padding]
#
# All modes print a ready-to-paste config.yaml docking: block.

set -euo pipefail

MODE="${1:-}"
PDB="${2:-}"

if [ -z "$MODE" ] || [ -z "$PDB" ] || [ ! -f "$PDB" ]; then
    echo "Usage:"
    echo "  bash $0 whole receptor.pdb [padding]"
    echo "  bash $0 residues receptor.pdb \"57,102,195\" [padding] [chain]"
    echo "  bash $0 fpocket receptor.pdb [padding]"
    exit 1
fi

# Shared awk routine: reads x,y,z triples from stdin (one per line),
# prints centroid + bounding-box (with padding) + a config.yaml snippet.
print_box() {
    local pad="$1"
    awk -v pad="$pad" '
    {
        x=$1+0; y=$2+0; z=$3+0;
        sx+=x; sy+=y; sz+=z; n++;
        if (n==1){minx=maxx=x; miny=maxy=y; minz=maxz=z}
        if (x<minx)minx=x; if (x>maxx)maxx=x;
        if (y<miny)miny=y; if (y>maxy)maxy=y;
        if (z<minz)minz=z; if (z>maxz)maxz=z;
    }
    END {
        if (n==0) { print "Error: no coordinates collected" > "/dev/stderr"; exit 1 }
        cx=sx/n; cy=sy/n; cz=sz/n;
        sizex=(maxx-minx)+2*pad; sizey=(maxy-miny)+2*pad; sizez=(maxz-minz)+2*pad;
        if (sizex<20) sizex=20; if (sizey<20) sizey=20; if (sizez<20) sizez=20;
        printf "Atoms used: %d  (padding %.1f A)\n\n", n, pad;
        printf "Centroid:\n  center_x: %.3f\n  center_y: %.3f\n  center_z: %.3f\n", cx, cy, cz;
        printf "\nBounding box + padding:\n  size_x: %.3f\n  size_y: %.3f\n  size_z: %.3f\n", sizex, sizey, sizez;
        printf "\n--- paste into config.yaml under docking: ---\n";
        printf "  center_x: %.2f\n  center_y: %.2f\n  center_z: %.2f\n", cx, cy, cz;
        printf "  size_x: %.1f\n  size_y: %.1f\n  size_z: %.1f\n", sizex, sizey, sizez;
    }'
}

case "$MODE" in
    whole)
        PADDING="${3:-5}"
        echo "Mode: whole-protein bounding box (blind docking)"
        echo "WARNING: this box covers the entire receptor. Docking accuracy"
        echo "is much lower than a focused box, and runtime increases a lot."
        echo "Use only as a last resort, or to do a first low-exhaustiveness"
        echo "pass to find candidate regions before refining the box."
        echo
        awk '$1=="ATOM" {print substr($0,31,8), substr($0,39,8), substr($0,47,8)}' "$PDB" \
            | print_box "$PADDING"
        ;;

    residues)
        RESLIST="${3:-}"
        PADDING="${4:-8}"
        CHAIN="${5:-}"
        if [ -z "$RESLIST" ]; then
            echo "Error: residues mode requires a comma-separated residue list, e.g. \"57,102,195\"" >&2
            exit 1
        fi
        echo "Mode: residue-based pocket (residues: $RESLIST${CHAIN:+, chain $CHAIN})"
        echo
        # Build an awk regex-friendly set from the comma list.
        RES_AWK=$(echo "$RESLIST" | tr ',' ' ')
        awk -v reslist="$RES_AWK" -v chain="$CHAIN" '
        BEGIN {
            n = split(reslist, arr, " ");
            for (i=1;i<=n;i++) want[arr[i]] = 1;
        }
        $1=="ATOM" {
            resnum = substr($0,23,4) + 0;
            ch = substr($0,22,1);
            if ((resnum in want) && (chain=="" || ch==chain)) {
                print substr($0,31,8), substr($0,39,8), substr($0,47,8);
            }
        }' "$PDB" | print_box "$PADDING"
        ;;

    fpocket)
        PADDING="${3:-4}"
        if ! command -v fpocket >/dev/null 2>&1; then
            echo "Error: fpocket not found. Install with:" >&2
            echo "  conda install -c bioconda -c conda-forge fpocket" >&2
            exit 1
        fi
        echo "Mode: fpocket cavity detection"
        echo "Running fpocket on $PDB ..."
        fpocket -f "$PDB"

        BASE="$(basename "$PDB" .pdb)"
        OUTDIR="$(dirname "$PDB")/${BASE}_out"
        POCKET1="${OUTDIR}/pockets/pocket1_atm.pdb"

        if [ ! -f "$POCKET1" ]; then
            echo "Error: expected $POCKET1 not found. Check $OUTDIR for fpocket output." >&2
            exit 1
        fi

        echo
        echo "Top-ranked pocket: $POCKET1"
        echo "(see ${OUTDIR}/${BASE}_info.txt for druggability scores of all pockets)"
        echo
        awk '$1=="ATOM"||$1=="HETATM" {print substr($0,31,8), substr($0,39,8), substr($0,47,8)}' "$POCKET1" \
            | print_box "$PADDING"
        ;;

    *)
        echo "Unknown mode: $MODE (expected: whole | residues | fpocket)" >&2
        exit 1
        ;;
esac
