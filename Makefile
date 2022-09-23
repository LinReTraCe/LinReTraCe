all: binR linretrace make_include

binR:
	if [ ! -d bin ] ; then mkdir bin ; fi

linretrace: make_include
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

test: bin/linretrace
	./ldft testsuite/Si --interp 3 --output testsuite/Si_input.hdf5
	./bin/linretrace testsuite/config.lrtc
	./lprint -p testsuite/Si_output.hdf5 c-intra xx yy zz
	./lprint -p testsuite/Si_output.hdf5 s-intra xx yy zz
	rm testsuite/Si_input.hdf5
	rm testsuite/Si_output.hdf5

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
