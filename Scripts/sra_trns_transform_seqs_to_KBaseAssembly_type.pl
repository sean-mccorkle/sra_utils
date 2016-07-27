#! /usr/bin/env perl

use strict;
use Getopt::Long;
use Data::Dumper;
use JSON;

use Bio::KBase::HandleService;

my $usage = <<"End_of_Usage";

Usage: $0 [ options ]

Upload one or two FASTA/FASTQ files to the shock server.

Options:

  -h, --help                          - print this help message and exit
  -s, --shock_service_url  URL        - shock service URL (D = https://kbase.us/services/shock-api)
  -n, --handle_service_url URL        - handle service URL (D = https://kbase.us/services/handle_service)
  -o, --output_file_name   json       - output JSON file of KBaseAssembly.PairedEndLibrary type
  --token                  string     - token string
  -t, --type               string     - output KBaseAssembly type (PairedEndLibrary, SingleEndLibrary, ReferenceAssembly)
  --sra                    int        - 1 means convert sra file into fastq files first (default 0, no sra conversion)

Options for PairedEndLibrary:

  -f, --input_file_name    path       - one or two read files (FASTA, FASTQ, or compressed forms)
  --insert                 float      - insert size mean
  --stdev                  float      - insert size standard deviation
  --outward                bool       - this flag is set to 1 if reads in the pair point outward

Options for SingleEndLibrary:

  -f, --input_file_name    path       - one read file (FASTA, FASTQ, or compressed forms)

Options for ReferenceAssembly:

  -f, --input_file_name    path       - one FASTA file containing a reference set of contigs
  --refname                text       - genome name of the reference contig set

Examples:

  $0 -t PairedEndLibrary -d stage/path -f read1.fq -f read2.fq -insert 300 -stdev 60 -o pe.reads.json
  $0 -t SingleEndLibrary -d stage/path -f read.fasta -o se.reads.json
  $0 -t ReferenceAssembly -d stage/path -f genome.fa -o ref.json

End_of_Usage

my ($help, $shock_url, $handle_url);
my ($type, $token);
my (@inputs, $output);
my ($insert, $stdev, $outward, $refname);
my $convert_sra = 0;
my $sra_convert_program = "fastq-dump";

my $rc = GetOptions("h|help"                 => \$help,
                    "s|shock_service_url=s"  => \$shock_url,
                    "n|handle_service_url=s" => \$handle_url,
                    "o|output_file_name=s"   => \$output,
                    "f|input_file_name=s"    => \@inputs,
                    "t|type=s"               => \$type,
                    "sra=i"                  => \$convert_sra,
                    "insert=f"               => \$insert,
                    "stdev=f"                => \$stdev,
                    "outward=i"              => \$outward,
                    "refname=s"              => \$refname);

$token      ||= $ENV{KB_AUTH_TOKEN};
$shock_url  ||= 'https://kbase.us/services/shock-api';
$handle_url ||= 'https://kbase.us/services/handle_service';

$help and die $usage;
$type && $output && @inputs >= 1 && @inputs <= 2 or die $usage;

my $shock = { url => $shock_url, token => $token };
my $handle_service = Bio::KBase::HandleService->new($handle_url);

my $obj;

if ($type eq 'PairedEndLibrary') {
    $obj = upload_pe_lib($shock, \@inputs, $insert, $stdev, $outward, $convert_sra);
} elsif ($type eq 'SingleEndLibrary') {
    $obj = upload_se_lib($shock, \@inputs, $convert_sra);
} elsif ($type eq 'ReferenceAssembly') {
    $obj = upload_ref($shock, \@inputs, $refname);
} else {
    die "Unrecognized output type: $type\n";
}

print_output($output, encode_json($obj));


# system_and_check( $cmd)
#    issue system( $cmd ) then check for error retur

sub  system_and_check
   {
    my $cmd = shift;
    print "### convert sra command is [$cmd]\n";
    system( $cmd );
    if ( $? == -1 )
       {  die "$cmd failed to execute: $!\n";  }
    elsif ( $? & 127 )
       {  printf STDERR "$cmd died with signal %d, %s coredump\n", 
               ($? & 127), ($? & 128 ) ? "with" : "without";
          die "$0 terminating\n";
       }
    print "### looks like successful return\n";
   }


# convert_sra( $filename, $mode) 
#    $filename is SRA filename (string)
#    $mode is 'se' for single end, 'pe' for paired end
# runs SRA conversion, converted files (fastq) will be left
# in working directory, returns converted file name(s) 
#  - with .fastq extention, in a list

sub  convert_sra
   {
    my ( $file, $mode ) = @_;

    my $fileroot = $file;
    $fileroot =~ s/.sra$//;  # remove any .sra extension

    if ( $mode eq 'se' )
       {
        my $outfile = $fileroot . ".fastq";
        system_and_check( "$sra_convert_program $file" );
        if ( -e $outfile )
           { return( ( $outfile ) ); }
        else
           { # maybe try a few other tricks here before bailing
             die "did not file expected $outfile from $sra_convert_program\n";
           }
       }
    elsif ( $mode eq 'pe' )
       {
        my @outfiles = map( $fileroot . "_" . $_ . ".fastq", (1,2) );
        system_and_check( "$sra_convert_program --split-files $file" );
        if ( -e $outfiles[0] && -e $outfiles[1] )
           { return( @outfiles ); }
        else
           { # maybe try a few other tricks here before bailing
             die "did not file expected @outfiles from $sra_convert_program\n";
           }
       }
   }


sub upload_pe_lib {
    my ($shock, $inputs, $insert, $stdev, $outward, $convert_sra) = @_;

    my $obj;

    if ( $convert_sra )
       {
        if ( @$inputs == 1 )
           {  @inputs  = convert_sra( $inputs[0], "pe" );  }
        else
           { die "$0: no handling for multiple SRA files\n"; }
        }
    $obj->{interleaved} = 1 if @$inputs == 1;
    $obj->{insert_size_mean} = $insert if $insert;
    $obj->{insert_size_std_dev} = $stdev if $stdev;
    $obj->{read_orientation_outward} = 1 if $outward;

    my $i;
    for my $file (@$inputs) {
        $obj->{'handle_'.++$i} = validate_seq_file($file);
    }

    return upload_files_in_obj($obj, $shock);
}

sub upload_se_lib {
    my ($shock, $inputs, $convert_sra) = @_;

    my $obj;
    my $file = $inputs->[0];

    if ( $convert_sra )
       {  ($file) = convert_sra( $file, "se" ); }

    $obj->{handle} = validate_seq_file($file);
    return upload_files_in_obj($obj, $shock);
}

sub upload_ref {
    my ($shock, $inputs, $refname) = @_;

    my $obj;
    my $file = $inputs->[0];
    $obj->{handle} = validate_seq_file($file);
    $obj->{reference_name} = $refname if $refname;

    return upload_files_in_obj($obj, $shock);
}

sub validate_seq_file {
    my ($file) = @_;
    -s $file or die "Invalid file: $file\n";
    my $name = $file;
    my $abs_file = $file;
    $name =~ s/\.(gz|gzip|bzip|bzip2|bz|bz2|zip)$//;
    $name =~ s/\.tar$//;
    $name =~ /\.(fasta|fastq|fas|fa|fq|fna|bas\.h5|bax\.h5)$/ or die "Unrecognized file type: $file\n";
    $file =~ s|.*/||;
    #return { file_name => $file, type => '_handle_to_be_' };
    return { full_path_file => $abs_file, file_name => $file, type => '_handle_to_be_' };
}

sub upload_files_in_obj {
    my ($obj, $shock) = @_;
    while (my ($k, $v) = each %$obj) {
        $obj->{$k} = update_handle($v, $shock) if is_handle($k, $v);
    }
    return $obj;
}

sub is_handle {
    my ($k, $v) = @_;
    return 1 if $k =~ /handle/ && $v && ref $v eq 'HASH' && $v->{file_name};
}

sub update_handle {
    my ($handle, $shock) = @_;

    #my $file = $handle->{file_name};
    my $file = $handle->{full_path_file};
    my $id = curl_post_file($file, $shock);


    undef $handle->{full_path_file};

    $handle->{type} = 'shock';
    $handle->{url}  = $shock->{url};
    $handle->{id}   = $id;

    if ($handle_service) {
        my $hid = $handle_service->persist_handle($handle);
        $handle->{hid} = $hid;
    }

    return $handle;
}

sub curl_post_file {
    my ($file, $shock) = @_;
    my $token = $shock->{token};
    my $url   = $shock->{url};
    my $attr  = q('{"filetype":"reads"}'); # should reference have a different type?
    my $cmd   = 'curl --connect-timeout 10 -s -X POST -F attributes=@- -F upload=@'.$file." $url/node ";
    $cmd     .= " -H 'Authorization: OAuth $token'";
    print "curl command is [$cmd]\n";
    my $out   = `echo $attr | $cmd` or die "Connection timeout uploading file to Shock: $file\n";
    print "curl output is [$cmd]\n";
    my $json  = decode_json($out);
    $json->{status} == 200 or die "Error uploading file: $file\n".$json->{status}." ".$json->{error}->[0]."\n";
    return $json->{data}->{id};
}

sub print_output {
    my ($fname, $text) = @_;
    open(F, ">$fname") or die "Could not open $fname";
    print F $text;
    close(F);
}
