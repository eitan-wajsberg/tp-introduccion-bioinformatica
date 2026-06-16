#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

# ---------------------------------------------------------------------------
# Ex5.pl
#
# Analiza la secuencia del gen HTT usando dos herramientas de EMBOSS:
#
#   1. getorf: encuentra todos los ORFs en la secuencia de nucleótidos
#              del mRNA (data/sequence.gb), los escribe en FASTA y guarda
#              el reporte en results/ex5_orfs.out.
#
#   2. patmatmotifs: busca motivos y dominios funcionales en la secuencia
#                   de aminoácidos (HTT_correcto.fas) contra la base de
#                   datos PROSITE. Guarda el reporte en results/ex5_dominios.out.
#
# Uso:
#   perl scripts/Ex5.pl <secuencia_aa.fas> [prosite.dat]
#
# Ejemplos:
#   perl scripts/Ex5.pl data/HTT_correcto.fas
#   perl scripts/Ex5.pl data/HTT_correcto.fas data/prosite.dat
#
# Outputs:
#   results/ex5_orfs.out       -> ORFs encontrados en el mRNA
#   results/ex5_dominios.out   -> dominios PROSITE encontrados en la proteína
# ---------------------------------------------------------------------------

# --- Argumentos ---
if (@ARGV < 1 || @ARGV > 2) {
    die "Uso: perl scripts/Ex5.pl <secuencia_aa.fas> [prosite.dat]\n" .
        "Ejemplo: perl scripts/Ex5.pl data/HTT_correcto.fas\n";
}

my $archivo_input  = $ARGV[0];
my $prosite_dat    = $ARGV[1] // "data/prosite.dat";
my $archivo_output = "results/ex5_dominios.out";
my $archivo_orfs   = "results/ex5_orfs.out";
my $archivo_gb     = "data/sequence.gb";
my $archivo_mrna   = "data/HTT_mrna.fasta";

system("mkdir -p results");

# ===========================================================================
# 1. Verificar que EMBOSS esté instalado
# ===========================================================================
print "=" x 60 . "\n";
print "Verificando instalación de EMBOSS...\n";
print "=" x 60 . "\n\n";

my $emboss_ok = (system("which patmatmotifs > /dev/null 2>&1") == 0);

unless ($emboss_ok) {
    print "EMBOSS no está instalado o no se encuentra en el PATH.\n\n";
    print "Para instalarlo:\n";
    print "  - Arch Linux:     wget -m 'ftp://emboss.open-bio.org/pub/EMBOSS/'\n";
    print "  - Ubuntu/Debian:  sudo apt-get install emboss\n";
    print "  - macOS (Homebrew): brew install emboss\n";
    print "  - Conda:          conda install -c bioconda emboss\n\n";
    die "Instalá EMBOSS y volvé a correr el script.\n";
}

my $version = `embossversion 2>/dev/null` || "desconocida";
chomp $version;
print "EMBOSS encontrado. Versión: $version\n\n";

# ===========================================================================
# 2. Convertir sequence.gb a FASTA de nucleótidos con seqret (si no existe)
# ===========================================================================
print "=" x 60 . "\n";
print "Preparando secuencia de nucleótidos del mRNA...\n";
print "=" x 60 . "\n\n";

unless (-e $archivo_mrna) {
    unless (-e $archivo_gb) {
        die "No se encontró el archivo GenBank: $archivo_gb\n" .
            "Asegurate de tener data/sequence.gb en el repositorio.\n";
    }

    print "Convirtiendo $archivo_gb a FASTA con seqret...\n";
    my $cmd_seqret = "seqret -sequence $archivo_gb -outseq $archivo_mrna -osformat fasta -auto";
    print "Comando: $cmd_seqret\n\n";

    my $ret_seqret = system($cmd_seqret);
    if ($ret_seqret != 0 || !-e $archivo_mrna || -z $archivo_mrna) {
        die "Error al convertir el archivo GenBank a FASTA.\n";
    }
    print "Archivo FASTA generado: $archivo_mrna\n\n";
} else {
    print "Archivo FASTA ya existe: $archivo_mrna. OK.\n\n";
}

# ===========================================================================
# 3. Correr getorf sobre la secuencia de nucleótidos
# ===========================================================================
print "=" x 60 . "\n";
print "Buscando ORFs con getorf...\n";
print "=" x 60 . "\n\n";

# getorf: encuentra todos los ORFs en los 6 marcos de lectura
# Flags usados:
#   -sequence  : archivo FASTA de nucleótidos de entrada
#   -outseq    : archivo FASTA de salida con las secuencias de los ORFs
#   -find      1: reportar ORFs que empiezan con Met y terminan en stop codon
#   -minsize 300: tamaño mínimo del ORF en nucleótidos (100 aminoácidos)
#   -auto      : no pedir confirmación interactiva
my $cmd_getorf = "getorf " .
                 "-sequence $archivo_mrna " .
                 "-outseq $archivo_orfs " .
                 "-find 1 " .
                 "-minsize 300 " .
                 "-auto";

print "Comando: $cmd_getorf\n\n";

my $ret_getorf = system($cmd_getorf);
if ($ret_getorf != 0) {
    die "Error al ejecutar getorf.\n";
}

unless (-e $archivo_orfs && -s $archivo_orfs) {
    die "No se generó el archivo de ORFs. " .
        "Verificá que la secuencia de entrada sea correcta.\n";
}

# Contar cuántos ORFs encontró
my $n_orfs = 0;
open(my $fh_orfs, '<', $archivo_orfs) or die "No puedo leer $archivo_orfs: $!\n";
while (<$fh_orfs>) { $n_orfs++ if /^>/ }
close $fh_orfs;

print "Se encontraron $n_orfs ORF(s) (con Met inicial, mínimo 100 aa).\n";
print "Reporte escrito en: $archivo_orfs\n\n";

# ===========================================================================
# 4. Verificar/descargar PROSITE
# ===========================================================================
print "=" x 60 . "\n";
print "Verificando base de datos PROSITE...\n";
print "=" x 60 . "\n\n";

my $prosite_doc = $prosite_dat;
$prosite_doc =~ s/prosite\.dat$/prosite.doc/;

unless (-e $prosite_dat && -e $prosite_doc) {
    print "No se encontraron los archivos de PROSITE.\n";
    print "Descargando prosite.dat y prosite.doc desde ExPASy FTP...\n";
    print "(Esto puede tardar unos minutos dependiendo de la conexión)\n\n";

    my $dir = dirname($prosite_dat);
    system("mkdir -p $dir") if $dir && $dir ne ".";

    my $use_wget = (system("which wget > /dev/null 2>&1") == 0);
    my $use_curl = (system("which curl > /dev/null 2>&1") == 0);

    unless ($use_wget || $use_curl) {
        die "No se encontró wget ni curl. Descargá manualmente desde:\n" .
            "  https://ftp.expasy.org/databases/prosite/prosite.dat\n" .
            "  https://ftp.expasy.org/databases/prosite/prosite.doc\n" .
            "y guardalos en: $dir/\n";
    }

    for my $archivo (["prosite.dat", $prosite_dat], ["prosite.doc", $prosite_doc]) {
        my ($nombre, $destino) = @$archivo;
        next if -e $destino;
        my $url = "https://ftp.expasy.org/databases/prosite/$nombre";
        print "Descargando $nombre...\n";
        my $ret;
        if ($use_wget) {
            $ret = system("wget -q --show-progress -O $destino $url");
        } else {
            $ret = system("curl -# -L -o $destino $url");
        }
        if ($ret != 0 || !-e $destino || -z $destino) {
            unlink $destino if -e $destino;
            die "Error al descargar $nombre. Chequeá tu conexión a internet.\n";
        }
        my $size_mb = sprintf("%.1f", (-s $destino) / 1_048_576);
        print "$nombre descargado correctamente (${size_mb} MB).\n";
    }
    print "\n";
} else {
    my $bytes = -s $prosite_dat;
    my $size_mb = defined($bytes) ? sprintf("%.1f", $bytes / 1_048_576) : "?";
    print "prosite.dat encontrado en '$prosite_dat' (${size_mb} MB). OK.\n";
    print "prosite.doc encontrado en '$prosite_doc'. OK.\n\n";
}

# ===========================================================================
# 5. Preprocesar PROSITE con prosextract
# ===========================================================================
print "=" x 60 . "\n";
print "Preprocesando PROSITE con prosextract...\n";
print "=" x 60 . "\n\n";

my $prosite_dir = $prosite_dat;
$prosite_dir =~ s|/[^/]+$||;
$prosite_dir = "." if $prosite_dir eq $prosite_dat;

$ENV{EMBOSS_DATA} = $prosite_dir;

my $prosite_out_dir = "$prosite_dir/PROSITE";
system("mkdir -p $prosite_out_dir");
print "Directorio de salida: $prosite_out_dir\n";

my $prosextract_cmd = "prosextract -prositedir $prosite_dir -auto";
print "Comando: $prosextract_cmd\n";
print "EMBOSS_DATA=$prosite_dir\n\n";

my $ret_prosextract = system($prosextract_cmd);
if ($ret_prosextract != 0) {
    die "Error al ejecutar prosextract.\n" .
        "Verificá que EMBOSS esté correctamente instalado y que prosite.dat sea válido.\n";
}
print "PROSITE preprocesado correctamente.\n\n";

# ===========================================================================
# 6. Correr patmatmotifs
# ===========================================================================
print "=" x 60 . "\n";
print "Corriendo análisis de dominios con patmatmotifs...\n";
print "=" x 60 . "\n\n";

unless (-e $archivo_input) {
    die "No se encontró el archivo de secuencias: $archivo_input\n";
}

my $cmd = "patmatmotifs " .
          "-sequence $archivo_input " .
          "-outfile $archivo_output " .
          "-full " .
          "-auto";

print "Comando: $cmd\n\n";

my $ret = system($cmd);

if ($ret != 0) {
    die "Error al ejecutar patmatmotifs.\n";
}

unless (-e $archivo_output && -s $archivo_output) {
    die "No se generó el archivo de salida.\n";
}

# ===========================================================================
# 7. Parsear y mostrar resumen en pantalla
# ===========================================================================
print "=" x 60 . "\n";
print "Resultados - Dominios PROSITE\n";
print "=" x 60 . "\n\n";

open(my $fh, '<', $archivo_output) or die "No puedo leer $archivo_output: $!\n";
my $contenido = do { local $/; <$fh> };
close $fh;

my @motivos;
while ($contenido =~ /Sequence:\s+(\S+).*?Motif = (\S+).*?Start = (\d+).*?End = (\d+)/gs) {
    push @motivos, {
        secuencia => $1,
        motivo    => $2,
        inicio    => $3,
        fin       => $4,
    };
}

if (@motivos == 0) {
    print "No se encontraron motivos PROSITE en la secuencia.\n\n";
} else {
    printf "Se encontraron %d motivo(s) PROSITE en la secuencia:\n\n", scalar @motivos;
    printf "%-5s  %-30s  %-10s  %-10s\n", "N°", "Motivo", "Inicio", "Fin";
    print "-" x 60 . "\n";
    my $i = 1;
    for my $m (@motivos) {
        printf "%-5d  %-30s  %-10d  %-10d\n",
            $i++, $m->{motivo}, $m->{inicio}, $m->{fin};
    }
    print "\n";
}

print "Reporte de ORFs escrito en:     $archivo_orfs\n";
print "Reporte de dominios escrito en: $archivo_output\n\n";
