#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Bio::SearchIO;

# --- Verificar argumentos ---
if (@ARGV != 2) {
    die "Uso: perl scripts/Ex2_local.pl <archivo.fas> <ruta/blast_db/swissprot>\n" .
        "Ejemplo: perl scripts/Ex2_local.pl data/HTT_correcto.fas blast_db/swissprot\n";
}

my $archivo_input  = $ARGV[0];
my $db_path        = $ARGV[1];
my $archivo_output = "results/blast_local.out";

# --- Verificar que blastp esté disponible ---
unless (system("which blastp > /dev/null 2>&1") == 0) {
    die "Error: blastp no encontrado en el PATH.\n";
}

# --- Verificar que la base de datos exista ---
unless (-e "$db_path.pin" || -e "$db_path.phr") {
    die "Error: No se encontró la base de datos en $db_path\n" .
        "Asegurate de haber corrido makeblastdb primero.\n";
}

print "Corriendo BLAST local contra: $db_path\n";
print "Secuencias input: $archivo_input\n\n";

# --- Ejecutar blastp desde la línea de comando ---
my $cmd = "blastp -query $archivo_input -db $db_path -out $archivo_output " .
          "-evalue 1e-5 -outfmt 0 -num_alignments 10 -num_descriptions 10";

my $resultado = system($cmd);

if ($resultado == 0) {
    print "BLAST local finalizado.\n";
    print "Resultados escritos en: $archivo_output\n\n";

    # --- Parsear y mostrar un resumen en pantalla ---
    my $searchio = Bio::SearchIO->new(
        -file   => $archivo_output,
        -format => 'blast'
    );

    print "=== Resumen de resultados ===\n\n";
    while (my $result = $searchio->next_result()) {
        print "Query: " . $result->query_name() . "\n";
        print "Base de datos: " . $result->database_name() . "\n\n";

        my $contador = 0;
        while (my $hit = $result->next_hit()) {
            $contador++;
            printf "  Hit %-3d: %s\n", $contador, $hit->name();
            printf "           Descripcion: %s\n", $hit->description();
            printf "           Score:       %s\n", $hit->score();
            printf "           E-value:     %s\n", $hit->significance();
            printf "           Identidad:   %.1f%%\n", $hit->frac_identical() * 100;
            print "\n";
        }
        print "Total hits encontrados: $contador\n";
    }
} else {
    die "Error al ejecutar BLAST local.\n";
}

