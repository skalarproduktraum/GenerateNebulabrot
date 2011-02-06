all: GenerateNebulabrot_GCG GenerateNebulabrot 

GenerateNebulabrot_GCD: GenerateNebulabrot.cxx
	$CXX -DUSE_GCD -O3 -fopenmp `freetype-config --cflags --libs` `libpng-config --cflags --libs` -L./ -lpngwriter -o $@ $^

GenerateNebulabrot: GenerateNebulabrot.cxx
	$CXX -O3 -fopenmp `freetype-config --cflags --libs` `libpng-config --cflags --libs` -L./ -lpngwriter -o $@ $^
