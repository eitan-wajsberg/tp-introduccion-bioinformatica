#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

# --- Verificar que se pase el archivo como argumento ---
if (@ARGV != 1) {
    die "Uso: perl scripts/Ex1.pl <archivo.gb>\n" .
        "Ejemplo: perl scripts/Ex1.pl data/sequence.gb\n";
}

my $archivo_input  = $ARGV[0];
my $archivo_output = "data/ORFs_HTT.fas";

# --- Leer el archivo GenBank ---
my $seqio_in = Bio::SeqIO->new(
    -file   => $archivo_input,
    -format => 'genbank'
);

my $seqio_out = Bio::SeqIO->new(
    -file   => ">$archivo_output",
    -format => 'fasta'
);

print "Leyendo archivo: $archivo_input\n";

while (my $seq = $seqio_in->next_seq()) {

    my $id     = $seq->id();
    my $largo  = $seq->length();
    print "Secuencia: $id | Largo: $largo bp\n\n";

    # --- Detectar el CDS anotado (marco de lectura correcto) ---
    my $cds_start  = undef;
    my $cds_end    = undef;
    my $cds_strand = undef;

    for my $feature ($seq->get_SeqFeatures()) {
        if ($feature->primary_tag eq 'CDS') {
            $cds_start  = $feature->start();
            $cds_end    = $feature->end();
            $cds_strand = $feature->strand();
        }
    }

    # --- Generar los 6 marcos de lectura ---
    # 3 en la hebra directa (+1, +2, +3)
    # 3 en la hebra complementaria inversa (-1, -2, -3)

    my @frames = (
        { offset => 0, strand => 1,  nombre => "+1" },
        { offset => 1, strand => 1,  nombre => "+2" },
        { offset => 2, strand => 1,  nombre => "+3" },
        { offset => 0, strand => -1, nombre => "-1" },
        { offset => 1, strand => -1, nombre => "-2" },
        { offset => 2, strand => -1, nombre => "-3" },
    );

    for my $frame (@frames) {
        my $hebra  = $frame->{strand};
        my $offset = $frame->{offset};
        my $nombre = $frame->{nombre};

        # Obtener la secuencia según la hebra
        my $subseq;
        if ($hebra == 1) {
            $subseq = $seq->subseq($offset + 1, $largo);
        } else {
            my $rev = $seq->revcom();
            $subseq = $rev->subseq($offset + 1, $largo);
        }

        # Traducir a proteína
        my $secuencia_obj = Bio::Seq->new(
            -seq      => $subseq,
            -alphabet => 'dna'
        );
        my $proteina = $secuencia_obj->translate()->seq();
        $proteina =~ s/\*.*//;

        # Determinar si este es el marco correcto según la anotación CDS
        my $es_correcto = "";
        if (defined $cds_start) {
            if ($hebra == 1 && $cds_strand == 1) {
                my $frame_num = ($cds_start - 1) % 3;
                $es_correcto = "_MARCO_CORRECTO" if $frame_num == $offset;
            } elsif ($hebra == -1 && $cds_strand == -1) {
                my $frame_num = ($largo - $cds_end) % 3;
                $es_correcto = "_MARCO_CORRECTO" if $frame_num == $offset;
            }
        }

        # Etiqueta para mostrar en terminal
        my $etiqueta = $es_correcto ? " <-- MARCO CORRECTO (CDS anotado)" : "";
        print "  Marco $nombre$etiqueta\n";

        # Escribir en el archivo FASTA
        my $seq_out = Bio::Seq->new(
            -id   => "${id}_frame_${nombre}${es_correcto}",
            -seq  => $proteina
        );
        $seqio_out->write_seq($seq_out);
    }
}

print "\nResultados escritos en: $archivo_output\n";
