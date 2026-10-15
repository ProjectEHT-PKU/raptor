#!/bin/bash
set -eu

echo "RAPTOR setup script"


METRIC="cks"
INT="rk4"
CODE="bhac"
GRID="smr"
RAD="pol"
SF="nosfc"
BFLIP="plus"
SIGMA_CUT=""
UNIFORM=""
THETAE_MIN=""
THETAE_MAX=""

for arg in "$@"
do
    case $arg in
        -c=*|--code=*)
        CODE="${arg#*=}"
        shift # Remove --cache= from processing
        ;;
	-m=*|--metric=*)
	METRIC="${arg#*=}"
	shift # Remove --cache= from processing
	;;
	-i=*|--int=*)
	INT="${arg#*=}"
        shift # Remove --cache= from processing
        ;;
	-g=*|--grid=*)
	GRID="${arg#*=}"
        shift # Remove --cache= from processing
        ;;
	-r=*|--rad=*)
        RAD="${arg#*=}"
        shift # Remove --cache= from processing
        ;;
       	-s=*|--sf=*)
        SF="${arg#*=}"
        shift # Remove --cache= from processing
        ;;
        -b=*|--bflip=*)
        BFLIP="${arg#*=}"
        shift # Remove --cache= from processing
        ;;
        --sigma-cut=*) SIGMA_CUT="${arg#*=}" ;;
        --uniform=*) UNIFORM="${arg#*=}" ;;
        --thetae-min=*) THETAE_MIN="${arg#*=}" ;;
        --thetae-max=*) THETAE_MAX="${arg#*=}" ;;
        *) echo "Unknown setup option: $arg" >&2; exit 1 ;;
    esac
done

# Validate before copying or changing any files. Keep the legacy short options.
case "$CODE" in bhac|harm3d|hamr) ;; *) echo "Invalid code: $CODE" >&2; exit 1 ;; esac
case "$METRIC" in cks|mks) ;; *) echo "Invalid metric: $METRIC" >&2; exit 1 ;; esac
case "$CODE:$METRIC" in harm3d:cks|hamr:cks) echo "CKS requires BHAC" >&2; exit 1 ;; esac
case "$INT" in ver|rk2|rk4|rk45) ;; *) echo "Invalid integrator: $INT" >&2; exit 1 ;; esac
case "$GRID" in smr|SMR) GRID=SMR ;; amr|AMR) GRID=AMR ;; *) echo "Invalid grid: $GRID" >&2; exit 1 ;; esac
case "$RAD" in pol|unpol) ;; *) echo "Invalid radiation: $RAD" >&2; exit 1 ;; esac
case "$SF" in sfc|nosfc) ;; *) echo "Invalid SFC: $SF" >&2; exit 1 ;; esac
case "$BFLIP" in plus|minus) ;; *) echo "Invalid polarity: $BFLIP" >&2; exit 1 ;; esac
case "$UNIFORM" in ''|0|1) ;; *) echo "Invalid uniform flag: $UNIFORM" >&2; exit 1 ;; esac
for value in "$SIGMA_CUT" "$THETAE_MIN" "$THETAE_MAX"; do
    if [ -n "$value" ]; then
        awk -v x="$value" 'BEGIN {
            if (x !~ /^[+]?[0-9]*[.]?[0-9]+([eE][+-]?[0-9]+)?$/ || x+0 <= 0) exit 1
        }' || { echo "Invalid positive numeric option: $value" >&2; exit 1; }
    fi
done
if [ -n "$THETAE_MIN$THETAE_MAX" ]; then
    awk -v lo="${THETAE_MIN:-1e-3}" -v hi="${THETAE_MAX:-100}" \
        'BEGIN { if (lo+0 >= hi+0) exit 1 }' || { echo "Invalid theta_e bounds" >&2; exit 1; }
fi

if [ "${RAPTOR:-}" == "" ] ;
then
	echo "Code dir is empty, set RAPTOR enviroment variable first"
	echo "export RAPTOR=/path/to/code"
	echo "Aborting..."
	exit 1
else
	echo "Code located at $RAPTOR"
	echo " "
fi

rsync -v "$RAPTOR"/run/* .

echo "Copying files"
if [ "$CODE" == "bhac" ] ;
then
        rsync -v "$RAPTOR"/model/bhac/* .

	       echo "done copying"

elif  [ "$CODE" == "harm3d" ] ;
then
        rsync -v "$RAPTOR"/model/harm3d/* .
elif  [ "$CODE" == "hamr" ] ;
then
        rsync -v "$RAPTOR"/model/HAMR/* .
else
 echo "no grmhd model picked. Aborting."
 exit 1
fi

if [ "$CODE" == "bhac" ] ;
then
	if [ "$METRIC" == "mks" ] ;
	then
		sed -i '/#define metric (/s/.*/#define metric (MKSBHAC)/' definitions.h
	fi
    	if [ "$METRIC" == "cks" ] ;
    	then
	    	sed -i  '/#define metric (/s/.*/#define metric (CKS)/' definitions.h
        fi
fi

if [ "$CODE" == "harm3d" ] ;
then
    	if [ "$METRIC" == "mks" ] ;
    	then
	    	sed -i  '/#define metric (/s/.*/#define metric (MKSHARM)/' definitions.h
        fi
fi

if [ "$INT" == "ver" ] ;
then
      	sed -i '/#define int_method (/s/.*/#define int_method (VER)/' definitions.h
fi


if [ "$INT" == "rk2" ] ;
then
      	sed -i '/#define int_method (/s/.*/#define int_method (RK2)/' definitions.h
fi

if [ "$INT" == "rk4" ] ;
then
      	sed -i  '/#define int_method (/s/.*/#define int_method (RK4)/' definitions.h
fi

if [ "$INT" == "rk45" ] ;
then
      	sed -i  '/#define int_method (/s/.*/#define int_method (RK45)/' definitions.h
fi

if [ "$GRID" == "SMR" ] ;
then

      	sed -i  '/#define AMR /s/.*/#define AMR 0/' definitions.h
      	sed -i  '/#define SMR /s/.*/#define SMR 1/' definitions.h
fi

if [ "$GRID" == "AMR" ] ;
then
    	sed -i  '/#define AMR /s/.*/#define AMR 1/' definitions.h
        sed -i  '/#define SMR /s/.*/#define SMR 0/' definitions.h
fi

if [ "$RAD" == "pol" ] ;
then
      	sed -i  '/#define POL (/s/.*/#define POL (1)/' definitions.h
fi

if [ "$RAD" == "unpol" ] ;
then
    	sed -i  '/#define POL (/s/.*/#define POL (0)/' definitions.h
fi

if [ "$SF" == "sfc" ] ;
then
    	sed -i  '/#define SFC /s/.*/#define SFC 1/' model_definitions.h
else
    	sed -i  '/#define SFC /s/.*/#define SFC 0/' model_definitions.h
fi

if [ "$BFLIP" == "minus" ] ;
then
       	sed -i  '/#define BPOL /s/.*/#define BPOL (MINUS)/' model_definitions.h
else
    	sed -i  '/#define BPOL /s/.*/#define BPOL (PLUS)/' model_definitions.h
fi


# Change only macros explicitly requested; the official model headers remain templates.
for setting in "SIGMA_CUT:$SIGMA_CUT" "UNIF:$UNIFORM" "THETAE_MIN:$THETAE_MIN" "THETAE_MAX:$THETAE_MAX"; do
    macro="${setting%%:*}"
    value="${setting#*:}"
    if [ -n "$value" ]; then
        grep -q "^#define $macro " definitions.h || { echo "Missing macro: $macro" >&2; exit 1; }
        sed -i "s/^#define $macro .*/#define $macro ($value)/" definitions.h
    fi
done

echo " "
echo "Succesfully initialized."
echo " "

echo "Start compilation"
echo " "

if ! { make clean && make all; };
then
    echo "Compilations unsuccesful"
    exit 1
else
    echo "Compilations succesful"
fi

echo " "
echo "Done initialization, use make all to recompile."

echo "Setup used code $CODE (-c/--code) metric $METRIC (-m/--metric) int $INT (-i/--int) grid $GRID (-g/--grid) rad $RAD (-r/--rad) sfc $SF (-s/--sf)"
