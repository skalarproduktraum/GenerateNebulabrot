all: GenerateNebulabrot_GCD GenerateNebulabrot GMPNebulabrot 

GenerateNebulabrot_GCD: GenerateNebulabrot.cxx
	$(CXX) -DUSE_GCD -O3 -fopenmp `freetype-config --cflags --libs` `libpng-config --cflags --libs` -lpngwriter -o build/$@ $^

GenerateNebulabrot: GenerateNebulabrot.cxx
	$(CXX) -O3 -fopenmp `freetype-config --cflags --libs` `libpng-config --cflags --libs` -lpngwriter -o build/$@ $^

GMPNebulabrot: GMPGenerateNebulabrot.cxx
	$(CXX) -g `freetype-config --cflags --libs` `libpng-config --cflags --libs` -lgmp -lmpc -lmpfr -lpngwriter -o build/$@ $^

clean: 
	rm -rf build/*
