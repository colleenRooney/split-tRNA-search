# build vfp
GPP = g++

vfp:	
	$(GPP) -o vfp main.cpp fasta.cpp matrix.cpp footprint.cpp

# write install dir to sts script
# install
install:
# clean 
clean: 
