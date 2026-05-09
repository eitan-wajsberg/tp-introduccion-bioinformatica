#!/bin/bash

echo "================================================"
echo "  Verificación de entorno - TP Bioinformática"
echo "================================================"
echo ""

OK="\e[32m✔\e[0m"
FAIL="\e[31m✘\e[0m"
WARN="\e[33m⚠\e[0m"
ERRORES=0

# --- Perl ---
echo "[ Perl ]"
if command -v perl &> /dev/null; then
    VERSION=$(perl -e 'print $^V')
    echo -e "  $OK Perl instalado: $VERSION"
else
    echo -e "  $FAIL Perl no encontrado"
    ((ERRORES++))
fi
echo ""

# --- Módulos BioPerl ---
echo "[ Módulos BioPerl ]"
MODULOS=(
    "Bio::SeqIO"
    "Bio::Tools::Run::RemoteBlast"
    "Bio::SearchIO"
    "Bio::Seq"
    "Bio::Tools::CodonTable"
)

for MOD in "${MODULOS[@]}"; do
    if perl -e "use $MOD;" &> /dev/null; then
        echo -e "  $OK $MOD"
    else
        echo -e "  $FAIL $MOD  <-- no instalado"
        ((ERRORES++))
    fi
done
echo ""

# --- BLAST ---
echo "[ BLAST+ ]"
if command -v blastp &> /dev/null; then
    BVERSION=$(blastp -version 2>&1 | head -n 1)
    echo -e "  $OK blastp encontrado: $BVERSION"
else
    echo -e "  $WARN blastp no encontrado en el PATH"
    echo -e "        (necesario solo para BLAST local, no para remoto)"
fi

if command -v makeblastdb &> /dev/null; then
    echo -e "  $OK makeblastdb encontrado"
else
    echo -e "  $WARN makeblastdb no encontrado en el PATH"
    echo -e "        (necesario solo para BLAST local)"
fi
echo ""

# --- Resumen ---
echo "================================================"
if [ $ERRORES -eq 0 ]; then
    echo -e "  $OK Todo listo para arrancar el TP"
else
    echo -e "  $FAIL Hay $ERRORES problema(s) que resolver antes de arrancar"
fi
echo "================================================"
