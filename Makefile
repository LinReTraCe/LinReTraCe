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
