package XGR_app;
use Mojo::Base 'Mojolicious';

# This method will run once at server start
sub startup {
	my $self = shift;
	
	$ENV{MOJO_REVERSE_PROXY} = 1;
	$self->config(
		hypnotoad => {
			listen  => ['http://*:3020'],
			workers => 8,
			keep_alive_timeout => 300,
			websocket_timeout => 600,
			proxy => 1
		}
	);
	
	# Documentation browser under "/perldoc"
  	$self->plugin('PODRenderer');

  	# Router
  	my $r = $self->routes;
	
	# Template names are expected to follow the template.format.handler scheme, with template defaulting to controller/action or the route name, format defaulting to html and handler to ep
	
  	# Normal route to controller
  	## Home
  	$r->get('/')->to(template=>'index', controller=>'action', action=>'index');
  	
  	## demo
  	$r->get('/demo')->to(template=>'demo', format=>'html', handler=>'ep', controller=>'action', action=>'index');
  	## about
  	$r->get('/about')->to(template=>'about', format=>'html', handler=>'ep', controller=>'action', action=>'index');

  	
  	## Enricher
  	### Genes
  	$r->get('/enricher/genes')->to(template=>'enricherGenes', format=>'html', handler=>'ep', controller=>'action', action=>'enricher_genes');
  	$r->post('/enricher/genes')->to(template=>'enricherGenes', format=>'html', handler=>'ep', controller=>'action', action=>'enricher_genes');
  	### SNPs
  	$r->get('/enricher/snps')->to(template=>'enricherSNPs', format=>'html', handler=>'ep', controller=>'action', action=>'enricher_snps');
  	$r->post('/enricher/snps')->to(template=>'enricherSNPs', format=>'html', handler=>'ep', controller=>'action', action=>'enricher_snps');
  	### Yours
  	$r->get('/enricher/yours')->to(template=>'enricherYours', format=>'html', handler=>'ep', controller=>'action', action=>'enricher_yours');
  	$r->post('/enricher/yours')->to(template=>'enricherYours', format=>'html', handler=>'ep', controller=>'action', action=>'enricher_yours');
  	
  	## Socialiser
  	### Genes
  	$r->get('/socialiser/genes')->to(template=>'socialiserGenes', format=>'html', handler=>'ep', controller=>'action', action=>'socialiser_genes');
  	$r->post('/socialiser/genes')->to(template=>'socialiserGenes', format=>'html', handler=>'ep', controller=>'action', action=>'socialiser_genes');
  	### SNPs
  	$r->get('/socialiser/snps')->to(template=>'socialiserSNPs', format=>'html', handler=>'ep', controller=>'action', action=>'socialiser_snps');
  	$r->post('/socialiser/snps')->to(template=>'socialiserSNPs', format=>'html', handler=>'ep', controller=>'action', action=>'socialiser_snps');
  	
  	## Subneter
  	### Genes
  	$r->get('/subneter/genes')->to(template=>'subneterGenes', format=>'html', handler=>'ep', controller=>'action', action=>'subneter_genes');
  	$r->post('/subneter/genes')->to(template=>'subneterGenes', format=>'html', handler=>'ep', controller=>'action', action=>'subneter_genes');
  	### SNPs
  	$r->get('/subneter/snps')->to(template=>'subneterSNPs', format=>'html', handler=>'ep', controller=>'action', action=>'subneter_snps');
  	$r->post('/subneter/snps')->to(template=>'subneterSNPs', format=>'html', handler=>'ep', controller=>'action', action=>'subneter_snps');
  	
  	
}

1;
