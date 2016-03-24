package XGR_app::Controller::Utils;
use POSIX qw(ceil);
use DBI;
use Mojo::DOM;
use Mojo::Base 'Mojolicious::Controller';

##################################
# Define and create files/folders
##################################
our $rootFolder = '/var/www/apps/XGR';
our $dataFolder = "$rootFolder/data";
mkdir $dataFolder unless -e $dataFolder;
our $tmpFolder = 'public/tmp';
mkdir $tmpFolder unless -e $tmpFolder;

##################################
# Connect to SQL
##################################
sub DBConnect {
    my $database = shift || 'ultraDDR';
    my $host = shift || '127.0.0.1';
    my $port = shift || '';
    my $username = shift || 'hfang';
    my $password = shift || '';

    my $dsn = "DBI:mysql:$database:$host:$port";
    my $dbh = DBI->connect( $dsn, $username, $password, { RaiseError => 1 } );

    return $dbh;
}

sub DBDisconnect {
    my $dbh = shift || return;
    $dbh->disconnect();
}

##################################
# Read/Write files
##################################
sub read_from_file {
	my $filename = shift;
	
	my $content='';
	open(IN,"$filename") or die "$!";
	while(<IN>){
		$content.=$_;
	}
	close IN;
	
	return $content;
}

sub export_to_file { 
    my $filename = shift;
    my $content = shift;
    
    open(OUT, ">>$filename") or die "$!";
    print OUT $content;
    close OUT;
}

1;
