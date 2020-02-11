default:
	R CMD check norMmix
	R CMD build norMmix --resave-data
	R CMD check norMmix_0.0-1.tar.gz
