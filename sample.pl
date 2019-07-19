#!/opt/perl/5.16.3/bin/perl -w

use strict;
use warnings;
use Text::ParseWords;

sub  EuclideanDist($$);

sub  Read_PCA_Coordinates($$);
sub  Compute_Distances($$$);
sub  Sample_Thru_Pairing($$$$$);
sub  Sample_Thru_DBCS($$$$$);
sub  DBCS_Initial_Seed($$$);
sub  DBCS_MaxMin($$$$); 
sub  DBCS_MaxSum($$$$);

###==============================###
###  Check the input arguments.  ###
###==============================###

my  $infile1;
my  $infile2;
my  $sampling_method;
my  $ratio;

if ( $#ARGV == 3 ) {

    $infile1          =  $ARGV[0];    #-- A file with activity outcomes and molecular fingerprints
    $infile2          =  $ARGV[1];    #-- A file with molecular descriptors.
    $sampling_method  =  $ARGV[2];    #-- A sampling method.
    $ratio            =  $ARGV[3];

} else {

    print STDERR "Usage: $0  FILENAME1  FILENAME2  METHOD  RATIO\n\n";
    print STDERR "  - FILENAME1 : for a file with activity outcomes and molecular properties.\n\n";
    print STDERR "  - FILENAME2 : for a file with molecular fingerprints.\n\n";
    print STDERR "  - METHOD    : a sampling method.\n";
    print STDERR "                1 => active-inactive matching.\n";
    print STDERR "                2 => dissimilarity-based (MaxMin)\n";
    print STDERR "                3 => dissimilarity-based (MaxSum)\n";
    print STDERR "  - RATIO     : Active-to-Inactive ratio (an integer of >=1).\n";
    print STDERR "\n";
    exit;
}

#$infile1          =  "AID_743139_molprop.tab";
#$infile2          =  "molprop_transformed.txt";
#$sampling_method  =  1;

my  @actives   = ();
my  @inactives = ();
my  @inactives_sampled =();


###===========================###
###  Read the property file.  ###
###===========================###

open(INFILE1, $infile1) || die "Can't open file $infile1 for reading:$!\n";

my $firstline = readline(INFILE1);

while(<INFILE1>) {

    # 0  CID
    # 1  Activity
    # 2  CovalentUnitCount
    # 3  MolecularWeight
    # 4  HeavyAtomCount
    # 5  RotatableBondCount
    # 6  Complexity
    # 7  XLogP
    # 8  TPSA
    # 9  HBondDonorCount
    # 10  HBondAcceptorCount
    chomp;
    my @data = parse_line("\t", 0, $_);

    if ( $data[7] ne "" ) {

        if ( ( $data[1] eq "active agonist" ) || ( $data[1] eq "active antagonist" ) ) {

            push( @actives, $data[0] );

        } elsif ( $data[1] eq "inactive" ) {

            push( @inactives, $data[0] );

        }

    }

}

close(INFILE1) || die "Can't close file $infile1: $!\n";

print STDERR "#-- Num. Actives   = ", $#actives   + 1, "\n";
print STDERR "#-- Num. Inactives = ", $#inactives + 1, "\n";



if ( $sampling_method == 1 ) {

    Sample_Thru_Pairing( \@inactives_sampled, \@inactives, \@actives, $infile2, $ratio );

} elsif ( $sampling_method == 2 ) {

    my $num_samples = $ratio * ( $#actives + 1 );
    Sample_Thru_DBCS( \@inactives_sampled, \@inactives, $infile2, $num_samples, "MaxMin" );

} elsif ( $sampling_method == 3 ) {

    my $num_samples = $ratio * ( $#actives + 1 );
    Sample_Thru_DBCS( \@inactives_sampled, \@inactives, $infile2, $num_samples, "MaxSum" );

} else {

    print STDERR "\n\nError: Sampling Method should be a value from 1 to 3\n\n";
    exit;

}


###===================================###
###  Generate hash for sampled cids.  ###
###===================================###

my %samples = ();

for ( my $i=0; $i<=$#actives; $i++) {

    $samples{ $actives[ $i ] } = 1;
}

for ( my $i=0; $i<=$#inactives_sampled; $i++ ) {

    $samples{ $inactives_sampled[ $i ] } = 1

}

my  @all_samples = keys %samples;

print STDERR "# Actives   = ", $#actives + 1, " CIDs\n";
print STDERR "# Inactvies = ", $#inactives_sampled + 1, " CIDs\n";
print STDERR "# Total     = ", $#all_samples + 1, "CIDs\n";

###=================================###
###  Read the property file again.  ###
###=================================###

open(INFILE1, $infile1) || die "Can't open file $infile1 for reading: $!\n";

$firstline = readline(INFILE1);

print $firstline;

while(<INFILE1>) {

    chomp;
    my @data = parse_line("\t", 0, $_); 

    if ( defined $samples{ $data[0] } ) {
        print $_, "\n";
    }

}

close(INFILE1) || die "Can't close file $infile1: $!\n";


#################################
###  Subroutine Definitions.  ###
#################################


###===================================================================###
###  Compute the Euclidean distance between two sets of coordinates.  ###
###===================================================================###

sub  EuclideanDist($$) {

    my ( $aref1, $aref2 ) = @_;

    my  $sum2 = 0;

    if ( $#{ $aref1 } == $#{ $aref2 } ) {

        for ( my $i=0; $i<=$#{ $aref1 }; $i++ ) {
            $sum2 += ( $aref1->[ $i ] - $aref2->[ $i ] )**2;
        } 

    } else {

        print "Error: the dimensions of the two points are not equal ( ", $#{$aref1}, " vs ", $#{$aref2}, " )", "\n";
        exit;

    }

    my $dist = sqrt( $sum2 );

    return $dist;

}



###==============================================================###
###  Find the compound most similar to the remaining compounds.  ###
###==============================================================###


sub  DBCS_Initial_Seed($$$) {

    my ( $aref_seed, $aref_cids, $href_dist ) = @_;

    my  %cid_sum_dist = ();

    for ( my $i=0; $i<=$#{ $aref_cids }; $i++ ) {

        my  $sum_dist  =  0;
        my  $cid1  =  $aref_cids->[ $i ];

        for ( my $j=0; $j<=$#{ $aref_cids }; $j++ ) {

            my  $cid2  =  $aref_cids->[ $j ];

            if ( $i != $j ) {

                if ( $cid1 < $cid2 ) {
                    $sum_dist += $href_dist->{ $cid1 }->{ $cid2 };
                } else {
                    $sum_dist += $href_dist->{ $cid2 }->{ $cid1 };
                }

            }

        }

        $cid_sum_dist{ $cid1 } = $sum_dist;

    }

    @{ $aref_cids } = sort { $cid_sum_dist{ $a } <=> $cid_sum_dist{ $b } } keys %cid_sum_dist;

    push( @{ $aref_seed }, shift( $aref_cids ) );

    print STDERR "The initial seed is CID ", $aref_seed->[0], "\n";

}





###=============================###
###  Read the descriptor file.  ###
###=============================###

sub  Read_PCA_Coordinates($$) {

    my  ( $href_coords, $file ) = @_;

    open(FILE, $file) || die "Can't open file $file for reading: $!\n";

    while(<FILE>) {

        chomp;
        my  @data = parse_line("\t", 0, $_);

        for ( my $i=1; $i<=$#data; $i++ ) {
            push( @{ $href_coords->{ $data[0] } }, $data[$i] );
        }

    }

    close(FILE) || die "Can't close file $file: $!\n";

}


###========================================================================###
###  Compute the distances among the inactives to find the seed compound.  ###
###========================================================================###

sub  Compute_Distances($$$) {

    my ( $href_dist, $aref_cids, $href_coords ) = @_;

    for ( my $i=0; $i<=$#{ $aref_cids }; $i++ ) {

        my  $cid1 = $aref_cids->[$i];

        for ( my $j=$i+1; $j<=$#{ $aref_cids }; $j++ ) {

            my  $cid2 = $aref_cids->[$j];

            my  $dist = EuclideanDist( $href_coords->{ $cid1 }, $href_coords->{ $cid2 } );

            if ( $cid1 < $cid2 ) {
                $href_dist->{ $cid1 }->{ $cid2 } = $dist;
            } else {
                $href_dist->{ $cid2 }->{ $cid1 } = $dist;
            }

        }
    }
}


sub  Compute_Sum_Distances($$$) {

    my ( $href_sum_dist, $aref_cids, $href_coords ) = @_;


    for ( my $i=0; $i<=$#{ $aref_cids }; $i++ ) {

        my  $sum_dist = 0;
        my  $cid1 = $aref_cids->[$i];

        for ( my $j=0; $j<=$#{ $aref_cids }; $j++ ) {

            my  $cid2 = $aref_cids->[$j];

            if ( $i != $j ) {  #-- Compute the distance between two "different" CIDs.

                my $dist = EuclideanDist( $href_coords->{ $cid1 }, $href_coords->{ $cid2 } );

                $sum_dist += $dist;

            }

        }

        $href_sum_dist->{ $cid1 } = $sum_dist;

    }

}



sub  Sample_Thru_Pairing($$$$$) {

    my ( $aref_selected, $aref_all, $aref_ref, $coord_file, $ratio ) = @_;


    my  %Coordinates   = ();
    my  %Sum_Distances = ();

    Read_PCA_Coordinates(\%Coordinates, $coord_file );

    Compute_Sum_Distances( \%Sum_Distances, $aref_ref, \%Coordinates );

    my  @refs_sorted = sort { $Sum_Distances{ $a } <=> $Sum_Distances{ $b } } keys %Sum_Distances;



    ###========================================================###
    ###  Compute the distances between actives and inactives.  ###
    ###========================================================###

    for ( my $k=1; $k<=$ratio; $k++ ) {

        for ( my $i=0; $i<=$#refs_sorted; $i++ ) {

            my $min_dist;
            my $idx_at_min_dist;

            my $cid1 = $refs_sorted[ $i ];

            for ( my $j=0; $j<=$#{ $aref_all }; $j++ ) {

                my $cid2 = $aref_all->[$j];

                my $dist = EuclideanDist( $Coordinates{ $cid1 }, $Coordinates{ $cid2 } );

                if ( ( ! defined $min_dist ) || ( $min_dist > $dist ) ) {

                    $min_dist = $dist;
                    $idx_at_min_dist = $j;

                }

            }

            my $seed = splice( @{ $aref_all }, $idx_at_min_dist, 1 );

            push( @{ $aref_selected }, $seed );

            print STDERR $#{ $aref_selected } + 1, "\t", $seed, "\t", $cid1, "\t", $min_dist, "\n";

        }

    }

}


sub  Sample_Thru_DBCS($$$$$) {

    my ( $aref_selected, $aref_all, $coord_file, $num_samples, $dist_measure ) = @_;

    my  %Coordinates = ();
    my  %Distances   = ();

    Read_PCA_Coordinates(\%Coordinates, $coord_file );

    Compute_Distances( \%Distances, $aref_all, \%Coordinates );

    DBCS_Initial_Seed( $aref_selected, $aref_all, \%Distances);

    if ( lc( $dist_measure ) eq "maxmin" ) {

        DBCS_MaxMin( $aref_selected, $aref_all, \%Distances, $num_samples );

    } elsif ( lc( $dist_measure ) eq "maxsum" ) {

        DBCS_MaxSum( $aref_selected, $aref_all, \%Distances, $num_samples );

    } else {

        print STDERR "Error: \$dist_measure = ", $dist_measure, "\n";

    }

}



###============================###
###  Loop over the inactives.  ###
###============================###

sub  DBCS_MaxMin($$$$) {

    my ( $aref_seed, $aref_cids, $href_dist, $num_seed ) = @_;

    while( @{ $aref_seed } < $num_seed ) {    # If the loop runs when @{ $aref_seed } equals $num_seed
                                              # one more sample will be added at the end, and the array has
        my  $max_min_dist;                    # one more element than $num_seed.
        my  $index_at_max;


        #-- Loop over the remaining compounds in the original set.

        for ( my $i=0; $i<=$#{ $aref_cids }; $i++ ) {   #-- Loop over the remaining compounds
                                                        #-- in the original set.
            my  $min_dist;
            my  $cid1  =  $aref_cids->[ $i ];

            for ( my $j=0; $j<=$#{ $aref_seed }; $j++ ) {  #-- Loop over the sampled compounds
                                                           #-- in the seed set.
                my  $dist;
                my  $cid2  =  $aref_seed->[ $j ];

                if ( $cid1 < $cid2 ) {
                    $dist = $href_dist->{ $cid1 }->{ $cid2 };
                } elsif ( $cid1 > $cid2 ) {
                    $dist = $href_dist->{ $cid2 }->{ $cid1 };
                } else {
                    print "Error: the distance value does not exist for CIDs $cid1 and $cid2.\n";
                    exit;
                }

                if ( ( ! defined $min_dist ) || ( $min_dist > $dist ) ) {
                    $min_dist = $dist;
                }

            }

            if ( ( ! defined $max_min_dist ) || ( $max_min_dist < $min_dist )  ) {

                $max_min_dist = $min_dist;
                $index_at_max = $i;

            }

        }

        my $seed = splice( @{ $aref_cids }, $index_at_max, 1 );

        push( @{ $aref_seed }, $seed );

        print STDERR $#{ $aref_seed }+ 1, "\t", $seed, "\t", $max_min_dist, "\n";

    }

}



###============================###
###  Loop over the inactives.  ###
###============================###

sub  DBCS_MaxSum($$$$) {

    my ( $aref_seed, $aref_cids, $href_dist, $num_seed ) = @_;

    while( @{ $aref_seed } < $num_seed ) {

        my  $max_sum_dist;
        my  $index_at_max;

        for ( my $i=0; $i<=$#{ $aref_cids }; $i++ ) {   # Loop over each remaining inactive.

            my  $sum_dist = 0;
            my  $cid1  =  $aref_cids->[ $i ];

            for ( my $j=0; $j<=$#{ $aref_seed }; $j++ ) {

                my  $cid2  =  $aref_seed->[ $j ];

                if ( $cid1 < $cid2 ) {
                    $sum_dist += $href_dist->{ $cid1 }->{ $cid2 };
                } else {
                    $sum_dist += $href_dist->{ $cid2 }->{ $cid1 };
                }

            }

            if ( ( ! defined $max_sum_dist ) || ( $max_sum_dist < $sum_dist )  ) {

                $max_sum_dist = $sum_dist;
                $index_at_max = $i;

            }

        }

        my $seed = splice( @{ $aref_cids }, $index_at_max, 1 );

        push( @{ $aref_seed }, $seed );

        print STDERR $#{ $aref_seed } + 1, "\t", $seed, "\t", $max_sum_dist / $#{ $aref_seed }, "\n";

    }

}



