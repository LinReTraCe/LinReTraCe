all: binR linretrace

binR:
	if [ ! -d bin ] ; then mkdir bin ; fi

linretrace:
	cd src/; make

clean:
	cd src/; make clean

pristine:
	cd src/; make clean
	rmdir bin

ctags:
	ctags -R --exclude=src/digamma
