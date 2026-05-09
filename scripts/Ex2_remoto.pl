#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Tools::Run::RemoteBlast;
use Bio::SearchIO;

# --- Verificar argumento ---
if (@ARGV != 1) {
    die "Uso: perl scripts/Ex2_remoto.pl <archivo.fas>\n" .
        "Ejemplo: perl scripts/Ex2_remoto.pl data/HTT_correcto.fas\n";
}

my $archivo_input  = $ARGV[0];
my $archivo_output = "results/blast_remoto.out";

# --- Configurar el BLAST remoto ---
my $blast = Bio::Tools::Run::RemoteBlast->new(
    -prog       => 'blastp',
    -data       => 'swissprot',
    -expect     => 1e-5,
    -readmethod => 'SearchIO',
);

# --- Leer las secuencias del archivo FASTA ---
my $seqio = Bio::SeqIO->new(
    -file   => $archivo_input,
    -format => 'fasta'
);

print "Enviando secuencias al servidor BLAST del NCBI...\n";
print "Esto puede tardar varios minutos. Por favor esperá.\n\n";

open(my $fh, '>', $archivo_output) or die "No puedo abrir $archivo_output: $!\n";

while (my $seq = $seqio->next_seq()) {
    print "Corriendo BLAST para: " . $seq->id() . "\n";

    my $submit = $blast->submit_blast($seq);
    print "Solicitud enviada. Esperando respuesta del servidor NCBI...\n";

    my $intentos = 0;
    WAIT: while (1) {
        sleep(10);
        $intentos++;
        print "  Intento $intentos de recuperar resultado...\n";

        for my $rid ($blast->each_rid()) {
            my $result = $blast->retrieve_blast($rid);

            if (!ref($result)) {
                next;
            }

            print $fh "=== Resultado BLAST Remoto ===\n";
            print $fh "Secuencia query: " . $result->query_description() . "\n";
            print $fh "Base de datos:   " . $result->database_name() . "\n\n";

            my $contador = 0;
            while (my $hit = $result->next_hit()) {
                $contador++;
                printf $fh "Hit %-3d: %s\n", $contador, $hit->name();
                printf $fh "  Descripcion: %s\n", $hit->description();
                printf $fh "  Score:       %s\n", $hit->score();
                printf $fh "  E-value:     %s\n", $hit->significance();
                printf $fh "  Identidad:   %.1f%%\n\n", $hit->frac_identical() * 100;
            }

            print $fh "Total hits: $contador\n\n";
            $blast->remove_rid($rid);
            print "Resultado recibido. $contador hits encontrados.\n";
            last WAIT;
        }

        if ($intentos >= 60) {
            print "Tiempo de espera agotado (10 minutos). El servidor no respondio.\n";
            last WAIT;
        }
    }
}

close($fh);
print "\nResultados escritos en: $archivo_output\n";

