
.PHONY: help clean veryclean

help :
	cat ./makefiles/README

clean :
	rm -f *.o *.mod

veryclean : clean
	make -C ext clean
	make -C aenet-1.0.0_src clean
