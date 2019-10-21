clean:
	rm -f *.log *.last *.txt
	rm -f *.png

realclean: clean
	rm -rf *flux *image *psf *residual *model 
	rm -rf *_fixvis.ms
	rm -rf *_fixvis.ms.*
	rm -rf *.ms.*
	rm -rf *.G0 *.G0.*
	rm -rf *.K *.K.*
	rm -rf *.G *.G.*
	rm -rf *.B *.B.*

clobber: realclean
	rm -rf *.ms
	rm -rf *.clean.*
	rm -rf *.dirty.*
	rm -rf TempLattice*

# -fin-
