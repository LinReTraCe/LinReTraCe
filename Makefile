all help default:
	@echo ""
	@echo "---------------------------------------------------------------"
	@echo "LinReTraCe Configuration, Validation, Compilation and Testsuite"
	@echo "---------------------------------------------------------------"
	@echo ""
	@echo "  Configuration:"
	@echo "    Before compilation please configure your setup by"
	@echo "    creating a custom make_include file in the root folder."
	@echo "    Examples are included in ./make_include_examples"
	@echo ""
	@echo "    FC : Fortran Compiler (ifort, mpiifort, gfortran)"
	@echo "    FCDG : Corresponding Single Core Fortran Compiler (ifort, gfortran)"
	@echo "    FFLAGS : Fortran Flags"
	@echo "    FPPFLAGS: Fortran Preprocessor Flags"
	@echo "    HDF5: HDF5 library includes ind links ( HDF5 > 1.12 )"
	@echo ""
	@echo "  - make validate:     validation of make_include file"
	@echo "  - make linretrace:   compilation of LinReTraCe"
	@echo "  - make install:      copy the binary into your ${HOME}/bin folder."
	@echo "  - make test:         run testsuite"
	@echo "  - make clean:        clean repository from temporary compilation objects."

linretrace: binR compile make_include

binR:
	if [ ! -d bin ] ; then mkdir bin ; fi

compile: make_include
	cd src_linretrace/; make

clean:
	cd src_linretrace/; make clean
	cd src_linretrace/digamma; make clean

pristine:
	cd src_linretrace/; make clean
	cd src_linretrace/digamma; make clean
	rmdir bin

ctags:
	ctags -R --exclude=src_linretrace/digamma --exclude=src_pp --exclude=src_linretrace/deprecated

validate: make_include
	cd testsuite/validation; make

validate-clean:
	cd testsuite/validation; make clean

test: bin/linretrace
	./ldft testsuite/tests/Si --interp 3 --output testsuite/tests/Si_input.hdf5
	./bin/linretrace testsuite/tests/config.lrtc
	./lprint -p testsuite/tests/Si_output.hdf5 c-intra xx yy zz
	./lprint -p testsuite/tests/Si_output.hdf5 s-intra xx yy zz

test-clean:
	rm -f testsuite/tests/Si_input.hdf5
	rm -f testsuite/tests/Si_output.hdf5

install: bin/linretrace
	@echo "Creating bin folder in \$$HOME:"
	mkdir -p ${HOME}/bin
	cp bin/linretrace ${HOME}/bin
	@echo "Installation complete."
	@echo
	@echo "In order to make linretrace available execute the following command:"
	@echo "  echo 'export PATH=\"${HOME}/bin:\$$PATH\"' >> ~/.bashrc"
	@echo
	@echo "and reload the bashrc:"
	@echo "  source \$$HOME/.bashrc"
	@echo
