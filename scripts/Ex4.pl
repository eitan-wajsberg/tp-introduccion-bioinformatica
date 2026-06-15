#!/usr/bin/perl
use strict;
use warnings;
use Bio::SearchIO;
use Bio::DB::GenBank;
use Bio::DB::SwissProt;
use Bio::SeqIO;

# ---------------------------------------------------------------------------
# Ex4.pl
#
# Parsea un reporte de salida de BLAST y busca, dentro de la descripcion
# de cada hit, una coincidencia con un "pattern" dado por parametro
# (por ejemplo: "Mus musculus").
#
# Para los hits que coincidan:
#   - Escribe un listado en texto plano con Accession, Descripcion,
#     Score y E-value.
#   - (Punto extra) Descarga la secuencia completa de cada hit desde
#     GenBank usando Bio::DB::GenBank y la escribe en un archivo FASTA.
#
# Uso:
#   perl Ex4.pl <archivo_blast.out> <pattern> [prefijo_salida]
#
# Ejemplo:
#   perl Ex4.pl results/blast_local.out "Mus musculus" results/ex4
#
# Salidas generadas:
#   <prefijo_salida>_hits.out       -> listado de hits que matchean
#   <prefijo_salida>_sequences.fas  -> secuencias completas (FASTA)
# ---------------------------------------------------------------------------

my ($blast_file, $pattern, $prefix) = @ARGV;

if (!defined $blast_file || !defined $pattern) {
    die "Uso: perl Ex4.pl <archivo_blast.out> <pattern> [prefijo_salida]\n";
}
$prefix ||= "results/ex4";

# ---------------------------------------------------------------------------
# 1. Parsear el reporte BLAST y buscar el pattern en la descripcion de
#    cada hit
# ---------------------------------------------------------------------------
my $in = Bio::SearchIO->new(-format => 'blast', -file => $blast_file);

my @matches;

while (my $result = $in->next_result) {
    while (my $hit = $result->next_hit) {
        my $description = $hit->description;

        if ($description =~ /\Q$pattern\E/i) {

            # Obtener el accession de forma robusta: si Bio::SearchIO no
            # lo pudo separar correctamente, se intenta extraer del
            # campo "name" del hit.
            my $accession = $hit->accession;

            if (!defined $accession || $accession eq '' || $accession eq 'unknown') {
                my $name = $hit->name;
                if ($name =~ /([A-Z0-9]+\.\d+)/) {
                    $accession = $1;
                } else {
                    $accession = $name;
                }
            }

            push @matches, {
                accession   => $accession,
                description => $description,
                score       => $hit->score,
                evalue      => $hit->significance,
            };
        }
    }
}

# ---------------------------------------------------------------------------
# 2. Escribir el listado de hits que matchean (texto plano)
# ---------------------------------------------------------------------------
my $list_file = "${prefix}_hits.out";
open(my $out, '>', $list_file) or die "No se puede escribir $list_file: $!\n";

print $out "Pattern buscado: $pattern\n";
print $out "Archivo analizado: $blast_file\n";
print $out "Hits encontrados: " . scalar(@matches) . "\n\n";

foreach my $m (@matches) {
    print $out "Accession: $m->{accession}\n";
    print $out "Descripcion: $m->{description}\n";
    print $out "Score: $m->{score}\n";
    print $out "E-value: $m->{evalue}\n";
    print $out "---\n";
}
close($out);

print "Listado de hits escrito en: $list_file\n";

if (@matches == 0) {
    print "No se encontraron hits que coincidan con el pattern '$pattern'.\n";
    exit 0;
}

# ---------------------------------------------------------------------------
# 3. Punto extra: descargar la secuencia completa de cada hit desde
#    GenBank y escribirla en formato FASTA
# ---------------------------------------------------------------------------
my $fasta_file = "${prefix}_sequences.fas";
my $seqio_out = Bio::SeqIO->new(-file => ">$fasta_file", -format => 'fasta');

my $sp = Bio::DB::SwissProt->new();
my $gb = Bio::DB::GenBank->new();

foreach my $m (@matches) {
    my $acc = $m->{accession};
    my $acc_query = $acc;

    # Si el accession tiene formato UniProt (ej. P42859.2, Q76P24.1),
    # Bio::DB::GenBank no lo reconoce directamente. En ese caso, se
    # consulta primero a UniProt (Bio::DB::SwissProt) para obtener,
    # a partir de sus cross-references, un accession de RefSeq/GenBank
    # equivalente, que es el que finalmente se le pasa a Bio::DB::GenBank.
    if ($acc =~ /^[OPQ][0-9][A-Z0-9]{3}[0-9](\.\d+)?$/
        || $acc =~ /^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}(\.\d+)?$/) {

        print "Accession $acc tiene formato UniProt. Buscando cross-reference RefSeq/EMBL...\n";

        (my $sp_acc = $acc) =~ s/\.\d+$//;

        my $sp_seq;
        eval {
            $sp_seq = $sp->get_Seq_by_acc($sp_acc);
        };

        if ($@ || !defined $sp_seq) {
            warn "  -> No se pudo consultar UniProt para $sp_acc: $@\n";
            next;
        }

        my $refseq_acc;
        my $embl_acc;
        for my $link ($sp_seq->annotation->get_Annotations('dblink')) {
            my $db = $link->database;
            if ($db eq 'RefSeq' && !$refseq_acc) {
                $refseq_acc = $link->primary_id;
            }
            elsif ($db eq 'EMBL' && !$embl_acc) {
                $embl_acc = $link->primary_id;
            }
        }

        if ($refseq_acc) {
            $acc_query = $refseq_acc;
        }
        elsif ($embl_acc) {
            $acc_query = $embl_acc;
        }
        else {
            warn "  -> No se encontro cross-reference a RefSeq/EMBL para $sp_acc. Se omite.\n";
            next;
        }

        print "  -> Cross-reference encontrada: $acc_query\n";
    }

    print "Descargando secuencia completa para accession: $acc (consulta GenBank: $acc_query) ...\n";

    my $seq;
    eval {
        $seq = $gb->get_Seq_by_acc($acc_query);
    };

    if ($@ || !defined $seq) {
        warn "  -> No se pudo descargar la secuencia de $acc_query: $@\n";
        next;
    }

    $seqio_out->write_seq($seq);
    print "  -> OK (" . $seq->length . " residuos)\n";

    # Buena practica: esperar un segundo entre requests a NCBI para no
    # saturar el servidor.
    sleep(1);
}

print "Secuencias completas escritas en: $fasta_file\n";
